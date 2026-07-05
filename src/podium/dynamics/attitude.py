"""Rigid-body rotational dynamics (truth model).

State: quaternion q (scalar-first, body -> reference) + body angular
velocity w [rad/s]. Euler's equations with a full inertia tensor and an
external body-frame torque; RK4 with quaternion renormalization per step
(the receipt suite pins energy and angular-momentum conservation, which
is where unnormalized quaternions and sloppy integrators show up).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from podium.core import quat

F64 = NDArray[np.float64]


def euler_deriv(w: F64, inertia: F64, torque_body: F64) -> F64:
    """Body angular acceleration: I w' = tau - w x (I w)."""
    out: F64 = np.linalg.solve(inertia, torque_body - np.cross(w, inertia @ w))
    return out


def step(q: F64, w: F64, inertia: F64, torque_body: F64, dt: float) -> tuple[F64, F64]:
    """One RK4 step of the coupled kinematics + Euler equations."""
    def f(qq: F64, ww: F64) -> tuple[F64, F64]:
        return quat.deriv(qq, ww), euler_deriv(ww, inertia, torque_body)

    k1q, k1w = f(q, w)
    k2q, k2w = f(q + 0.5 * dt * k1q, w + 0.5 * dt * k1w)
    k3q, k3w = f(q + 0.5 * dt * k2q, w + 0.5 * dt * k2w)
    k4q, k4w = f(q + dt * k3q, w + dt * k3w)
    q_out = quat.normalize(q + (dt / 6.0) * (k1q + 2 * k2q + 2 * k3q + k4q))
    w_out = w + (dt / 6.0) * (k1w + 2 * k2w + 2 * k3w + k4w)
    return q_out, w_out


def gravity_gradient_torque(nadir_body: F64, inertia: F64, n: float) -> F64:
    """Gravity-gradient torque in the body frame for a circular orbit,

        tau = 3 n^2 (o_hat x I o_hat),

    where o_hat is the unit NADIR direction (toward the primary's
    centre) expressed in body coordinates and n is the orbital mean
    motion. Vanishes when a principal axis points at nadir (the
    gravity-gradient equilibrium); to first order it restores the body
    toward the local-vertical/local-horizontal orientation, the basis
    of passive gravity-gradient stabilization."""
    o = np.asarray(nadir_body, dtype=np.float64)
    o = o / np.linalg.norm(o)
    tau: F64 = 3.0 * n * n * np.cross(o, inertia @ o)
    return tau


def kinetic_energy(w: F64, inertia: F64) -> float:
    return 0.5 * float(w @ (inertia @ w))


def momentum_inertial(q: F64, w: F64, inertia: F64) -> F64:
    """Angular momentum rotated into the reference frame (conserved when
    the external torque is zero)."""
    return quat.rotate(q, inertia @ w)
