"""Analytic torque-free attitude oracle (#43): the Euler+quaternion
integrator is cross-validated against EXACT closed-form solutions of
rigid-body rotation — a mathematical ground truth, not another
numerical stack.

- Asymmetric body: Euler's equations have a Jacobi-elliptic-function
  solution; omega(t) must match to ~1e-11.
- Axisymmetric body: closed-form regular precession — the symmetry
  axis traces a constant-angle nutation cone about the conserved
  angular momentum H and precesses about it at |H|/It; omega_perp
  rotates at the body rate n(Ia-It)/It.
"""

import numpy as np
import pytest

pytest.importorskip("scipy")
from scipy.special import ellipj  # noqa: E402

from podium.core import quat  # noqa: E402
from podium.dynamics import attitude as att  # noqa: E402


def test_asymmetric_omega_matches_jacobi_elliptic():
    """General torque-free tumble (I1<I2<I3): the integrated angular
    velocity equals the exact Jacobi-elliptic solution to ~1e-11."""
    i1, i2, i3 = 80.0, 100.0, 120.0
    inertia = np.diag([i1, i2, i3])
    iv = np.array([i1, i2, i3])
    # start on the analytic t=0 state (cn=dn=1, sn=0): w=(a1,0,a3)
    w0_seed = np.array([0.2, 0.15, 0.6])
    two_t = float(np.sum(iv * w0_seed**2))
    l2 = float(np.sum(iv**2 * w0_seed**2))
    assert l2 > two_t * i2   # motion about the major axis (this case)
    a1 = np.sqrt((two_t * i3 - l2) / (i1 * (i3 - i1)))
    a2 = np.sqrt((two_t * i3 - l2) / (i2 * (i3 - i2)))
    a3 = np.sqrt((l2 - two_t * i1) / (i3 * (i3 - i1)))
    rate = np.sqrt((i3 - i2) * (l2 - two_t * i1) / (i1 * i2 * i3))
    m = ((i2 - i1) * (two_t * i3 - l2)) / ((i3 - i2) * (l2 - two_t * i1))

    q = np.array([1.0, 0.0, 0.0, 0.0])
    w = np.array([a1, 0.0, a3])
    dt = 0.01
    err = 0.0
    for k in range(1, 4001):
        q, w = att.step(q, w, inertia, np.zeros(3), dt)
        sn, cn, dn, _ = ellipj(rate * (k * dt), m)
        wa = np.array([a1 * cn, a2 * sn, a3 * dn])
        err = max(err, float(np.max(np.abs(w - wa))))
    assert err < 1e-11, err


def test_axisymmetric_omega_matches_regular_precession():
    """Axisymmetric torque-free body: omega_3 is constant and the
    transverse omega rotates at the body-precession rate lam =
    n(Ia-It)/It — the closed form."""
    it, ia = 100.0, 60.0
    inertia = np.diag([it, it, ia])
    n, wp = 0.5, 0.2
    lam = n * (ia - it) / it
    q = np.array([1.0, 0.0, 0.0, 0.0])
    w = np.array([wp, 0.0, n])
    dt = 0.005
    err = 0.0
    for k in range(1, 4001):
        q, w = att.step(q, w, inertia, np.zeros(3), dt)
        t = k * dt
        wa = np.array([wp * np.cos(lam * t), wp * np.sin(lam * t), n])
        err = max(err, float(np.max(np.abs(w - wa))))
    assert err < 1e-9, err


def test_axisymmetric_nutation_cone_and_precession_rate():
    """The attitude predictions: the symmetry axis keeps a constant
    angle to the conserved H (nutation cone) and precesses about H at
    the inertial rate |H|/It. H itself is conserved (torque-free)."""
    it, ia = 100.0, 60.0
    inertia = np.diag([it, it, ia])
    w = np.array([0.2, 0.0, 0.5])
    q = np.array([1.0, 0.0, 0.0, 0.0])
    h0 = inertia @ w                 # inertial ang. momentum (q0 = I)
    hmag = float(np.linalg.norm(h0))
    hhat = h0 / hmag
    ref = np.cross(hhat, [1.0, 0.0, 0.0])
    ref /= np.linalg.norm(ref)
    ref2 = np.cross(hhat, ref)
    dt = 0.005
    az, az_prev, cos_lo, cos_hi, hdev = 0.0, None, 1e9, -1e9, 0.0
    for _ in range(4000):
        q, w = att.step(q, w, inertia, np.zeros(3), dt)
        h = quat.rotate(q, inertia @ w)
        hdev = max(hdev, float(np.linalg.norm(h - h0)))
        ax = quat.rotate(q, np.array([0.0, 0.0, 1.0]))
        cosn = float(np.dot(ax, hhat))
        cos_lo, cos_hi = min(cos_lo, cosn), max(cos_hi, cosn)
        perp = ax - cosn * hhat
        a = float(np.arctan2(np.dot(perp, ref2), np.dot(perp, ref)))
        if az_prev is not None:
            d = (a - az_prev + np.pi) % (2 * np.pi) - np.pi
            az += d
        az_prev = a
    assert hdev < 1e-10, hdev                    # momentum conserved
    assert cos_hi - cos_lo < 1e-9                 # nutation cone constant
    measured = az / (4000 * dt)
    assert abs(measured - hmag / it) < 0.01 * (hmag / it)  # precession


def test_torque_free_energy_and_momentum_conserved():
    """Reaffirm the physical invariants the analytic solutions assume:
    kinetic energy and the inertial momentum magnitude are constant."""
    inertia = np.diag([80.0, 100.0, 120.0])
    w = np.array([0.3, 0.2, 0.5])
    q = np.array([1.0, 0.0, 0.0, 0.0])
    e0 = att.kinetic_energy(w, inertia)
    l0 = float(np.linalg.norm(att.momentum_inertial(q, w, inertia)))
    dt = 0.01
    for _ in range(5000):
        q, w = att.step(q, w, inertia, np.zeros(3), dt)
    assert abs(att.kinetic_energy(w, inertia) - e0) < 1e-9
    lf = float(np.linalg.norm(att.momentum_inertial(q, w, inertia)))
    assert abs(lf - l0) < 1e-9
