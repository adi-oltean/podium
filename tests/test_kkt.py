"""Verified-KKT-checker receipts (#40): an untrusted solver's solution
is re-verified for optimality in EXACT rational arithmetic.

The headline receipt is exactness — a QP with a known rational optimum
verifies to residuals that are Fraction(0) EXACTLY, so the checker
certifies true optima with zero slack, not merely 'small'. Then an
untrusted-solver solution of a Layer-0 min-energy rendezvous QP
certifies within tolerance, and perturbing the primal or flipping a
dual sign makes the exact residuals blow up so the checker rejects.
"""

import math
from fractions import Fraction as F

import numpy as np
import pytest

from podium import constants as const
from podium.core import cw
from podium.verify import kkt


def test_equality_qp_exact_optimum_zero_residual():
    """min 1/2(x1^2 + x2^2) s.t. x1 + x2 = 2  ->  x=(1,1), nu=-1.
    Every KKT residual is EXACTLY zero."""
    p = [[F(1), F(0)], [F(0), F(1)]]
    q = [F(0), F(0)]
    a = [[F(1), F(1)]]
    b = [F(2)]
    rep = kkt.verify_qp(p, q, [], [], a, b,
                        x=[F(1), F(1)], mu=[], nu=[F(-1)])
    assert rep.stationarity == F(0)
    assert rep.eq_residual == F(0)
    assert rep.duality_gap == F(0)
    assert rep.primal_obj == F(1)   # 1/2(1+1)
    assert rep.certified(tol=F(0))  # exact, not just within tolerance


def test_inequality_qp_exact_optimum_zero_residual():
    """min 1/2 x^2 s.t. x >= 1 (as -x <= -1) -> x=1, mu=1, active.
    Complementary slackness and stationarity hold EXACTLY."""
    rep = kkt.verify_qp(
        p=[[F(1)]], q=[F(0)], g=[[F(-1)]], h=[F(-1)], a=[], b=[],
        x=[F(1)], mu=[F(1)], nu=[])
    assert rep.stationarity == F(0)
    assert rep.ineq_violation == F(0)
    assert rep.dual_violation == F(0)
    assert rep.duality_gap == F(0)     # mu * slack = 1 * 0
    assert rep.certified(tol=F(0))


def test_perturbed_primal_is_rejected():
    """Nudging x off the optimum makes the exact stationarity/eq
    residual large — the checker rejects it."""
    p = [[F(1), F(0)], [F(0), F(1)]]
    a = [[F(1), F(1)]]
    bad = kkt.verify_qp(p, [F(0), F(0)], [], [], a, [F(2)],
                        x=[F(11, 10), F(1)], mu=[], nu=[F(-1)])
    assert not bad.certified()
    # x1+x2 = 2.1 != 2, and stationarity broken
    assert bad.eq_residual == F(1, 10)
    assert bad.stationarity >= F(1, 10)


def test_flipped_dual_sign_is_rejected():
    """A negative multiplier is dual-INFEASIBLE; the checker catches it
    exactly even though stationarity might look plausible."""
    rep = kkt.verify_qp(
        p=[[F(1)]], q=[F(0)], g=[[F(-1)]], h=[F(-1)], a=[], b=[],
        x=[F(1)], mu=[F(-1)], nu=[])
    assert rep.dual_violation == F(1)
    assert not rep.certified()


def test_nonsymmetric_p_is_flagged():
    rep = kkt.verify_qp(
        p=[[F(1), F(2)], [F(0), F(1)]], q=[F(0), F(0)],
        g=[], h=[], a=[], b=[], x=[F(0), F(0)], mu=[], nu=[])
    assert any("symmetric" in m for m in rep.problems)
    assert not rep.certified()


@pytest.mark.slow
def test_untrusted_solver_min_energy_rendezvous_certified():
    """A Layer-0 min-energy rendezvous QP (min 1/2||u||^2 s.t. the
    impulses reach the target state) solved by the UNTRUSTED cvxpy/
    Clarabel path, then verified exactly. The rationalized solution
    certifies within a tight tolerance."""
    cp = pytest.importorskip("cvxpy")
    n = math.sqrt(const.MU_EARTH / 6_778_137.0**3)
    ts = [0.0, 300.0, 600.0]
    tf = 900.0
    bmat = np.vstack([np.zeros((3, 3)), np.eye(3)])  # dv -> state
    # A (6 x 9): final state = sum stm(tf - t_k) @ B @ u_k
    blocks = [cw.stm(n, tf - t) @ bmat for t in ts]
    amat = np.hstack(blocks)
    target = np.array([0.0, -5.0, 0.0, 0.0, 0.0, 0.0])

    u = cp.Variable(9)
    prob = cp.Problem(cp.Minimize(0.5 * cp.sum_squares(u)),
                      [amat @ u == target])
    prob.solve(solver=cp.CLARABEL)
    assert prob.status == "optimal"
    nu_dual = prob.constraints[0].dual_value

    rep = kkt.verify_qp(
        p=kkt.rationalize_mat(np.eye(9).tolist()),
        q=[kkt.Frac(0)] * 9,
        g=[], h=[],
        a=kkt.rationalize_mat(amat.tolist()),
        b=kkt.rationalize_vec(target.tolist()),
        x=kkt.rationalize_vec(u.value.tolist()),
        mu=[],
        nu=kkt.rationalize_vec(nu_dual.tolist()))
    # exact residuals from a float solver are tiny but nonzero
    assert rep.certified(tol=kkt.Frac(1, 10**5)), (
        float(rep.stationarity), float(rep.eq_residual))
    # sanity: the certified objective matches the solver's
    assert abs(float(rep.primal_obj) - prob.value) < 1e-6
