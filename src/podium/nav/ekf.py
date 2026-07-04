"""Relative-navigation Kalman filter (Joseph form, static-subset style).

Kernels are pure functions with fixed shapes and closed-form models —
the shape the C emitter expects. The covariance update uses the Joseph
form P = (I-KH) P (I-KH)' + K R K', which preserves symmetry and
positive semidefiniteness under floating-point roundoff; the naive
(I-KH) P form does not, and a filter that slowly loses PD-ness fails in
the field, not in the unit test. The innovation and its covariance are
returned so callers can run consistency monitors (NEES/NIS) — the
receipts treat statistical consistency as a hard requirement, not a
tuning nicety.

Process noise: discrete white-noise-acceleration model. q_accel is the
acceleration PSD [ (m/s^2)^2 * s ]; it must absorb both real actuation
noise and the CW-vs-truth model error at the chosen rate (quantified in
the tests).
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from podium.core import cw
from podium.verify import Interval, contract, shapes

F64 = NDArray[np.float64]

# position-only and position+velocity measurement matrices
H_POS = np.hstack([np.eye(3), np.zeros((3, 3))])
H_POSVEL = np.eye(6)


@contract(dt=Interval(1e-3, 600.0), q_accel=Interval(0.0, 1.0))
def process_noise_wna(dt: float, q_accel: float) -> F64:
    """Discrete white-noise-acceleration Q (6x6, per-axis blocks)."""
    q = np.zeros((6, 6))
    q_pp = q_accel * dt * dt * dt / 3.0
    q_pv = q_accel * dt * dt / 2.0
    q_vv = q_accel * dt
    for i in range(3):
        q[i, i] = q_pp
        q[i, i + 3] = q_pv
        q[i + 3, i] = q_pv
        q[i + 3, i + 3] = q_vv
    return q


@shapes(x=(6,), p=(6, 6), phi=(6, 6), q=(6, 6))
def predict(x: F64, p: F64, phi: F64, q: F64) -> tuple[F64, F64]:
    """KF time update: x+ = Phi x, P+ = Phi P Phi' + Q (symmetrized).

    Emits to C as written (matmul lowering); NumPy's BLAS accumulation
    order differs from the emitted naive loops, so golden vectors for
    this kernel are a relative-tolerance class (see tests)."""
    x_out: F64 = phi @ x
    p_prop = phi @ p @ phi.T + q
    p_out: F64 = 0.5 * (p_prop + p_prop.T)
    return x_out, p_out


@contract(r_var=Interval(1e-2, 1e6))
@shapes(x=(6,), p=(6, 6), z=(3,))
def update_sequential(x: F64, p: F64, z: F64, r_var: float) -> tuple[F64, F64]:
    """Sequential scalar Joseph updates for position measurements
    (H = [I3 0], R = r_var * I3).

    The flight-side form: one measurement component at a time, so the
    innovation covariance is a SCALAR (division, not linalg.solve — no
    matrix factorization in flight code). For diagonal R this is
    algebraically equivalent to the batch Joseph update, which the
    receipts verify to near machine precision. Static subset: fixed
    shapes, bounded loops, one division guarded by r_var > 0.
    """
    xs = np.empty(6)
    ps = np.empty((6, 6))
    for i in range(6):
        xs[i] = x[i]
        for j in range(6):
            ps[i, j] = p[i, j]
    pn = np.empty((6, 6))
    kk = np.empty(6)
    for m in range(3):
        s = ps[m, m] + r_var
        # covariance-repair clamp: a valid covariance has ps[m,m] >= 0
        # so s >= r_var already; the clamp makes the division PROVABLY
        # safe for any input (EVA discharges it from the r_var contract)
        if s < r_var:  # noqa: PLR1730 — max() is outside the C subset
            s = r_var
        for i in range(6):
            kk[i] = ps[i, m] / s
        nu = z[m] - xs[m]
        for i in range(6):
            xs[i] = xs[i] + kk[i] * nu
        # Joseph form for H = e_m: P' = (I - k e_m') P (I - k e_m')'
        #                              + r k k'
        for i in range(6):
            for j in range(6):
                pn[i, j] = (ps[i, j] - kk[i] * ps[m, j]
                            - ps[i, m] * kk[j]
                            + kk[i] * ps[m, m] * kk[j]
                            + r_var * kk[i] * kk[j])
        for i in range(6):
            for j in range(6):
                ps[i, j] = 0.5 * (pn[i, j] + pn[j, i])
    return xs, ps


def _joseph_core(
    x: F64, p: F64, nu: F64, h: F64, r: F64
) -> tuple[F64, F64, F64]:
    s = h @ p @ h.T + r
    k = np.linalg.solve(s.T, (p @ h.T).T).T  # K = P H' S^-1
    x_out: F64 = x + k @ nu
    ikh = np.eye(len(x)) - k @ h
    p_j = ikh @ p @ ikh.T + k @ r @ k.T
    p_out: F64 = 0.5 * (p_j + p_j.T)
    return x_out, p_out, s


def update_joseph(
    x: F64, p: F64, z: F64, h: F64, r: F64
) -> tuple[F64, F64, F64, F64]:
    """Joseph-form measurement update (linear H).

    Returns (x, P, innovation, innovation_covariance). The gain solve is
    against S = H P H' + R (never inverted explicitly).
    """
    nu = z - h @ x
    x_out, p_out, s = _joseph_core(x, p, nu, h, r)
    return x_out, p_out, nu, s


def update_joseph_nonlinear(
    x: F64,
    p: F64,
    z: F64,
    h_fn: Callable[[F64], F64],
    h_jac: Callable[[F64], F64],
    r: F64,
    angle_rows: tuple[int, ...] = (),
) -> tuple[F64, F64, F64, F64]:
    """EKF update: nonlinear measurement prediction h(x), Jacobian H(x),
    Joseph-form covariance. angle_rows lists innovation components that
    are angles — they are wrapped to (-pi, pi] before the update, which
    matters for bearing sensors near the +/-pi seam."""
    h = np.asarray(h_jac(x), dtype=np.float64)
    pred = np.asarray(h_fn(x), dtype=np.float64)
    nu = z - pred
    for i in angle_rows:
        nu[i] = math.atan2(math.sin(nu[i]), math.cos(nu[i]))
    x_out, p_out, s = _joseph_core(x, p, nu, h, r)
    return x_out, p_out, nu, s


class RelNavEkf:
    """CW-model relative-navigation filter (sandbox convenience wrapper).

    Impulsive burns are known inputs: pass the commanded dv to step() and
    it feeds through the velocity estimate at the tick, so maneuvers are
    not mistaken for state error.
    """

    def __init__(
        self,
        n: float,
        dt: float,
        q_accel: float,
        r_pos: float,
        r_vel: float | None = None,
        x0: F64 | None = None,
        p0: F64 | None = None,
    ) -> None:
        self.phi = cw.stm(n, dt)
        self.q = process_noise_wna(dt, q_accel)
        if r_vel is None:
            self.h = H_POS
            self.r = np.eye(3) * r_pos * r_pos
        else:
            self.h = H_POSVEL
            self.r = np.diag([r_pos**2] * 3 + [r_vel**2] * 3)
        self.x: F64 = np.zeros(6) if x0 is None else np.asarray(x0, dtype=np.float64).copy()
        self.p: F64 = (
            np.diag([100.0**2] * 3 + [1.0**2] * 3) if p0 is None else p0.copy()
        )
        self.last_nu: F64 = np.zeros(self.h.shape[0])
        self.last_s: F64 = np.eye(self.h.shape[0])

    def step(self, z: F64, dv: F64 | None = None) -> F64:
        """One measurement-then-predict cycle; returns the post-update
        estimate valid at the measurement time. dv (if any) is the
        commanded impulse applied at this tick, fed through prediction."""
        self.x, self.p, self.last_nu, self.last_s = update_joseph(
            self.x, self.p, z, self.h, self.r
        )
        est = self.x.copy()
        x_pred = self.x.copy()
        if dv is not None:
            x_pred[3:6] += dv
        self.x, self.p = predict(x_pred, self.p, self.phi, self.q)
        return est
