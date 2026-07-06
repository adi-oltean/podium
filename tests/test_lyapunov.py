"""Certificate-carrying Lyapunov ellipsoid receipts (#50): an LQR
controller's value matrix P, re-verified in exact rational arithmetic
as a closed-loop Lyapunov certificate, and shown to bound and contract
the trajectory it controls.
"""

import math
from fractions import Fraction as F

import numpy as np

from podium.control import lqr
from podium.verify import lyapunov as ly

MU = 3.986004418e14
N = math.sqrt(MU / 6_778_137.0**3)


def _cw_lqr():
    a, b = lqr.cw_discrete(N, 5.0)
    q = np.diag([1.0, 1.0, 1.0, 100.0, 100.0, 100.0])
    r = np.eye(3) * 1e3
    k, p = lqr.dlqr_cert(a, b, q, r)
    a_cl = a - b @ k
    return a, b, k, p, a_cl


def test_lqr_value_matrix_certifies_the_closed_loop():
    """The rationalized Riccati P passes the exact Lyapunov check:
    P >= 0 and P - A_cl' P A_cl >= 0 (the decrease carries the full Q
    margin, so exact verification is robust to rationalization)."""
    _a, _b, _k, p, a_cl = _cw_lqr()
    rep = ly.verify_lyapunov(ly.rationalize_matrix(a_cl.tolist()),
                             ly.rationalize_matrix(p.tolist()))
    assert rep.certified(), (rep.p_positive, rep.decrease_psd,
                             rep.problems)


def test_value_decreases_monotonically_along_trajectory():
    """x' P x is strictly non-increasing under the closed loop from a
    random start — the ellipsoid contracts."""
    _a, _b, _k, p, a_cl = _cw_lqr()
    rng = np.random.default_rng(50)
    x = rng.normal(0.0, np.array([50, 50, 50, 0.1, 0.1, 0.1]))
    v_prev = float(x @ p @ x)
    for _ in range(300):
        x = a_cl @ x
        v = float(x @ p @ x)
        assert v <= v_prev + 1e-9 * abs(v_prev), (v, v_prev)
        v_prev = v
    assert v_prev < 1e-3 * float(np.max(np.abs(p)))  # converged toward 0


def test_ellipsoid_sublevel_set_is_invariant():
    """If x0 is inside the ellipsoid {x'Px <= c}, every future state
    stays inside — the defining property of an invariant set."""
    _a, _b, _k, p, a_cl = _cw_lqr()
    rng = np.random.default_rng(7)
    x = rng.normal(0.0, np.array([30, 30, 30, 0.05, 0.05, 0.05]))
    c = float(x @ p @ x)
    for _ in range(500):
        x = a_cl @ x
        assert float(x @ p @ x) <= c + 1e-9 * c


def test_exact_certificate_object_and_negatives():
    """The exact certificate rejects a bogus P (identity is not a
    Lyapunov matrix for this open-loop-unstable plant) and a P that
    breaks the decrease."""
    a, _b, _k, _p, a_cl = _cw_lqr()
    # identity P: the decrease I - A_cl' A_cl need not be PSD
    ident = [[F(1) if i == j else F(0) for j in range(6)] for i in range(6)]
    bogus = ly.verify_lyapunov(ly.rationalize_matrix(a_cl.tolist()), ident)
    assert not bogus.certified()
    # the OPEN-loop CW plant is not contracting: P from the closed loop
    # fails the Lyapunov decrease against the open-loop A
    _a2, _b2, _k2, p, _acl2 = _cw_lqr()
    open_loop = ly.verify_lyapunov(ly.rationalize_matrix(a.tolist()),
                                   ly.rationalize_matrix(p.tolist()))
    assert not open_loop.decrease_psd

    # the EllipsoidInvariant value is exact
    inv = ly.EllipsoidInvariant(ident)
    assert inv.value([F(2), F(0), F(0), F(0), F(0), F(0)]) == F(4)


def test_zero_p_is_rejected_not_positive_definite():
    """P = 0 is PSD and trivially non-increasing, but {x : x'Px <= c} is
    all of R^n, not a bounded ellipsoid. The certificate requires P > 0
    (positive definite), so the checker must reject P = 0."""
    _a, _b, _k, _p, a_cl = _cw_lqr()
    zero = [[F(0)] * 6 for _ in range(6)]
    rep = ly.verify_lyapunov(ly.rationalize_matrix(a_cl.tolist()), zero)
    assert rep.decrease_psd          # 0 - A'0A = 0 is trivially PSD...
    assert not rep.p_positive        # ...but P = 0 is not positive definite
    assert not rep.certified()
