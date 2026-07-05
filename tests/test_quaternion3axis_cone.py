"""Higher-degree three-axis quaternion-feedback safe-attitude-cone
barrier (#54) --- the last SOS frontier the paper named.

The full three-axis attitude closed loop (7 states q0..q3, w1..w3;
torque tau = -kp q_vec - kd w, asymmetric inertia) is a genuinely
nonlinear polynomial system on the unit-quaternion manifold. Two exact
results, both verified over the rationals:

1. V-dot = -kd ||w||^2 for the Lyapunov barrier V = 2kp(1-q0) +
   1/2 w' I w --- the kp coupling AND the gyroscopic w x I w term both
   cancel (the latter because w . (w x I w) = 0, so rigid-body
   rotation is workless);
2. a POSITIVSTELLENSATZ containment {V <= c} subset {q0 >= q0_min}
   (a safe attitude cone) via the exact identity
   q0 - q0_min = (1/2kp)(c - V) + (1/4kp) w' I w with the remainder SOS.

Together: starting inside {V <= c} the attitude never leaves the cone,
for all time --- an infinite-horizon safe-cone certificate for the
full three-axis nonlinear loop.
"""

from fractions import Fraction as F

import numpy as np
import pytest

from podium.verify import sos

_VARS = "q0 q1 q2 q3 w1 w2 w3".split()


def _m(**kw):
    e = [0] * 7
    for k, v in kw.items():
        e[_VARS.index(k)] = v
    return tuple(e)


def _system(kp, kd, i1, i2, i3):
    f = [
        {_m(q1=1, w1=1): F(-1, 2), _m(q2=1, w2=1): F(-1, 2),
         _m(q3=1, w3=1): F(-1, 2)},
        {_m(q0=1, w1=1): F(1, 2), _m(q2=1, w3=1): F(1, 2),
         _m(q3=1, w2=1): F(-1, 2)},
        {_m(q0=1, w2=1): F(1, 2), _m(q1=1, w3=1): F(-1, 2),
         _m(q3=1, w1=1): F(1, 2)},
        {_m(q0=1, w3=1): F(1, 2), _m(q1=1, w2=1): F(1, 2),
         _m(q2=1, w1=1): F(-1, 2)},
        {_m(q1=1): -kp / i1, _m(w1=1): -kd / i1,
         _m(w2=1, w3=1): -(i3 - i2) / i1},
        {_m(q2=1): -kp / i2, _m(w2=1): -kd / i2,
         _m(w3=1, w1=1): -(i1 - i3) / i2},
        {_m(q3=1): -kp / i3, _m(w3=1): -kd / i3,
         _m(w1=1, w2=1): -(i2 - i1) / i3},
    ]
    v = {_m(): 2 * kp, _m(q0=1): -2 * kp,
         _m(w1=2): i1 / 2, _m(w2=2): i2 / 2, _m(w3=2): i3 / 2}
    return f, v


def test_vdot_kp_and_gyroscopic_terms_cancel_exactly():
    """V-dot = -kd(w1^2 + w2^2 + w3^2) for several asymmetric inertias
    and gains: both the kp coupling and the gyroscopic w x I w term
    cancel exactly."""
    for kp, kd, i1, i2, i3 in [(F(2), F(1), F(2), F(3), F(4)),
                               (F(5), F(3), F(1), F(5, 2), F(7)),
                               (F(3, 2), F(2), F(4), F(4), F(1))]:
        f, v = _system(kp, kd, i1, i2, i3)
        vdot = sos.lie_derivative(v, f)
        assert vdot == {_m(w1=2): -kd, _m(w2=2): -kd, _m(w3=2): -kd}


def test_neg_vdot_is_sos():
    """-V-dot = kd||w||^2 is SOS (basis [w1,w2,w3])."""
    f, v = _system(F(2), F(1), F(2), F(3), F(4))
    neg = sos.pscale(F(-1), sos.lie_derivative(v, f))
    ok, prob = sos.is_sos(neg, [_m(w1=1), _m(w2=1), _m(w3=1)],
                          [[F(1), F(0), F(0)], [F(0), F(1), F(0)],
                           [F(0), F(0), F(1)]])
    assert ok, prob


def _cone_remainder(kp, kd, i1, i2, i3, q0_min):
    _f, v = _system(kp, kd, i1, i2, i3)
    c = 2 * kp * (1 - q0_min)
    cmv = sos.psub({_m(): c}, v)
    return sos.psub(sos.psub({_m(q0=1): F(1)}, {_m(): q0_min}),
                    sos.pscale(F(1, 2) / kp, cmv)), c


def test_safe_cone_positivstellensatz_exact():
    """The containment {V<=c} subset {q0>=q0_min}: the exact remainder
    q0 - q0_min - (1/2kp)(c - V) equals (1/4kp) w' I w, which is SOS ->
    on {V<=c} (c-V>=0), q0 - q0_min = nonneg + nonneg >= 0."""
    kp, kd, i1, i2, i3, q0_min = F(2), F(1), F(2), F(3), F(4), F(4, 5)
    rem, _c = _cone_remainder(kp, kd, i1, i2, i3, q0_min)
    assert rem == {_m(w1=2): i1 / (4 * kp), _m(w2=2): i2 / (4 * kp),
                   _m(w3=2): i3 / (4 * kp)}
    ok, prob = sos.is_sos(rem, [_m(w1=1), _m(w2=1), _m(w3=1)],
                          [[i1 / (4 * kp), F(0), F(0)],
                           [F(0), i2 / (4 * kp), F(0)],
                           [F(0), F(0), i3 / (4 * kp)]])
    assert ok, prob


@pytest.mark.slow
def test_cone_remainder_validated_from_untrusted_sdp():
    """The remainder's SOS certificate synthesized by an UNTRUSTED float
    SDP over [w1,w2,w3] is validated to an exact rational certificate on
    the 7-state system."""
    cp = pytest.importorskip("cvxpy")
    kp, kd, i1, i2, i3, q0_min = F(2), F(1), F(2), F(3), F(4), F(4, 5)
    rem, _c = _cone_remainder(kp, kd, i1, i2, i3, q0_min)
    basis = [_m(w1=1), _m(w2=1), _m(w3=1)]
    diag = [float(i1 / (4 * kp)), float(i2 / (4 * kp)), float(i3 / (4 * kp))]

    g = cp.Variable((3, 3), symmetric=True)
    cons = [g >> 0, g[0, 0] == diag[0], g[1, 1] == diag[1],
            g[2, 2] == diag[2], g[0, 1] == 0, g[0, 2] == 0, g[1, 2] == 0]
    cp.Problem(cp.Minimize(0), cons).solve(solver=cp.CLARABEL)
    g_exact = sos.validate_gram(rem, basis, g.value.tolist())
    assert g_exact is not None
    ok, prob = sos.is_sos(rem, basis, g_exact)
    assert ok, prob


def test_three_axis_loop_stays_in_cone():
    """Numerically: from inside {V<=c} the full three-axis loop keeps
    q0 >= q0_min for all time, preserves ||q||, and converges."""
    kp, kd = 2.0, 1.0
    inrt = np.array([2.0, 3.0, 4.0])
    q0_min = 0.8
    c = 2 * kp * (1 - q0_min)

    def field(x):
        q, w = x[:4], x[4:]
        qd = 0.5 * np.array([
            -(q[1]*w[0] + q[2]*w[1] + q[3]*w[2]),
            q[0]*w[0] + q[2]*w[2] - q[3]*w[1],
            q[0]*w[1] - q[1]*w[2] + q[3]*w[0],
            q[0]*w[2] + q[1]*w[1] - q[2]*w[0]])
        gyro = np.cross(w, inrt * w)
        wd = (-kp * q[1:4] - kd * w - gyro) / inrt
        return np.concatenate([qd, wd])

    def energy(x):
        return 2*kp*(1 - x[0]) + 0.5*float(inrt @ x[4:]**2)

    ang = 0.5   # start well inside the cone
    axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
    q = np.concatenate([[np.cos(ang/2)], np.sin(ang/2)*axis])
    x = np.concatenate([q, np.array([0.05, -0.03, 0.04])])
    assert energy(x) <= c              # start inside {V<=c}
    dt, norm0 = 0.01, float(q @ q)
    for _ in range(8000):
        k1 = field(x)
        k2 = field(x + 0.5*dt*k1)
        k3 = field(x + 0.5*dt*k2)
        k4 = field(x + dt*k3)
        x = x + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)
        assert x[0] >= q0_min - 1e-6            # never leaves the cone
        assert abs(float(x[:4] @ x[:4]) - norm0) < 1e-6   # ||q|| kept
    assert np.linalg.norm(x[4:]) < 1e-2         # converged
