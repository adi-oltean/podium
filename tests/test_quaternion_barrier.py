"""Exact SOS barrier for the polynomial quaternion-feedback closed loop
(#53) --- the RPOD-core nonlinear system the paper names.

Single-axis attitude with quaternion feedback: state (q_s, q_v, w),
torque tau = -kp q_v - kd w. The closed loop is a genuinely nonlinear
(quadratic vector field) polynomial system on the unit-quaternion
manifold. Its Lyapunov barrier V = 2 kp (1 - q_s) + 1/2 I w^2 has
V-dot = -kd w^2 --- an exact cancellation of the kp cross terms ---
certifying an infinite-horizon attitude-stability invariant, verified
over the rationals, and validated from an untrusted float SDP.
"""

from fractions import Fraction as F

import numpy as np
import pytest

from podium.verify import sos


def _closed_loop(kp: F, kd: F, inv_i: F):
    """Polynomials for the single-axis quaternion-feedback closed loop
    (vars q_s=x0, q_v=x1, w=x2) and its Lyapunov barrier V."""
    f = [
        {(0, 1, 1): F(-1, 2)},                       # q_s dot = -1/2 w q_v
        {(1, 0, 1): F(1, 2)},                         # q_v dot =  1/2 w q_s
        {(0, 1, 0): -inv_i * kp, (0, 0, 1): -inv_i * kd},  # w dot
    ]
    v = {(0, 0, 0): 2 * kp, (1, 0, 0): -2 * kp, (0, 0, 2): F(1, 2) / inv_i}
    return f, v


def test_lie_derivative_cross_terms_cancel_exactly():
    """For any rational gains, V-dot = -kd w^2 exactly: the kp q_v w
    cross terms cancel --- the kinematics/dynamics coupling, verified
    over the rationals."""
    for kp, kd, inv_i in [(F(2), F(1), F(1)), (F(5), F(3), F(1, 4)),
                          (F(7, 2), F(2), F(2, 3))]:
        f, v = _closed_loop(kp, kd, inv_i)
        vdot = sos.lie_derivative(v, f)
        assert vdot == {(0, 0, 2): -kd}, (kp, kd, inv_i, vdot)


def test_barrier_lie_derivative_is_sos():
    """-V-dot = kd w^2 is SOS (basis [w]) -> V is non-increasing along
    the closed loop, so every sub-level set is an infinite-horizon
    invariant of the nonlinear quaternion system."""
    f, v = _closed_loop(F(2), F(1), F(1))
    neg_vdot = sos.pscale(F(-1), sos.lie_derivative(v, f))
    ok, prob = sos.is_sos(neg_vdot, [(0, 0, 1)], [[F(1)]])
    assert ok, prob


@pytest.mark.slow
def test_barrier_certificate_validated_from_untrusted_sdp():
    """The SOS certificate of -V-dot synthesized by an UNTRUSTED float
    SDP over a richer basis [w, q_v, q_s] is validated to an EXACT
    rational certificate --- the validated-SOS pipeline on the real
    RPOD closed loop."""
    cp = pytest.importorskip("cvxpy")
    f, v = _closed_loop(F(2), F(1), F(1))
    target = sos.pscale(F(-1), sos.lie_derivative(v, f))   # kd w^2
    basis = [(0, 0, 1), (0, 1, 0), (1, 0, 0)]              # w, q_v, q_s

    g = cp.Variable((3, 3), symmetric=True)
    cons = [g >> 0]
    # coefficient constraints: -V-dot = kd w^2 = z^T G z
    cons += [g[0, 0] == 1,                # w^2
             g[1, 1] == 0, g[2, 2] == 0,  # q_v^2, q_s^2
             g[0, 1] == 0, g[0, 2] == 0, g[1, 2] == 0]
    cp.Problem(cp.Minimize(0), cons).solve(solver=cp.CLARABEL)

    g_exact = sos.validate_gram(target, basis, g.value.tolist())
    assert g_exact is not None
    ok, prob = sos.is_sos(target, basis, g_exact)
    assert ok, prob


def test_closed_loop_converges_and_preserves_norm():
    """Numerically: the quaternion-feedback loop drives (q_v, w) -> 0,
    V is monotone non-increasing, and the quaternion norm q_s^2 + q_v^2
    is preserved (the certified invariant is physical)."""
    kp, kd, inv_i = 2.0, 1.0, 1.0

    def field(x):
        qs, qv, w = x
        return np.array([-0.5 * w * qv, 0.5 * w * qs,
                         inv_i * (-kp * qv - kd * w)])

    def energy(x):
        return 2 * kp * (1 - x[0]) + 0.5 / inv_i * x[2]**2

    ang = 1.2
    x = np.array([np.cos(ang / 2), np.sin(ang / 2), 0.3])
    v_prev, norm0, dt = energy(x), x[0]**2 + x[1]**2, 0.01
    for _ in range(6000):
        k1 = field(x)
        k2 = field(x + 0.5 * dt * k1)
        k3 = field(x + 0.5 * dt * k2)
        k4 = field(x + dt * k3)
        x = x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        v = energy(x)
        assert v <= v_prev + 1e-9
        assert abs((x[0]**2 + x[1]**2) - norm0) < 1e-6   # ||q|| preserved
        v_prev = v
    assert abs(x[1]) < 1e-2 and abs(x[2]) < 1e-2         # converged
