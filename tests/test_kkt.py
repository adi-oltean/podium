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


# --- SOCP (#41): exact cone complementarity ---------------------------

def test_soc_membership_predicate_is_exact():
    """The sqrt-free 'within tol of the second-order cone' test, on
    exact points: interior, boundary, exterior, and the tolerance
    boundary — all decided in exact Fraction arithmetic."""
    ok = kkt._soc_margin_ok
    assert ok([F(2), F(1), F(1)], F(0))         # ||(1,1)||<2 interior
    assert ok([F(1), F(1), F(0)], F(0))         # ||1||=1=s0 boundary
    assert not ok([F(1), F(1), F(1)], F(0))     # ||(1,1)||=sqrt2>1 out
    assert not ok([F(-1), F(0), F(0)], F(0))    # s0<0 out
    # exterior point pulled in by a large enough tolerance:
    # ||(1,1)||^2 = 2 <= (1+tol)^2 requires tol >= sqrt2 - 1 ~ 0.4142
    assert not ok([F(1), F(1), F(1)], F(41, 100))
    assert ok([F(1), F(1), F(1)], F(42, 100))


def test_lp_as_socp_exact_zero_residual():
    """A pure-nonneg-cone standard form (an LP) with a rational
    optimum: min x s.t. x >= 1 -> x=1, s=0, z=1. Every conic-KKT
    residual is EXACTLY zero and both cones hold at tol 0."""
    rep = kkt.verify_socp(
        c=[F(1)], a=[], b=[], g=[[F(-1)]], h=[F(-1)],
        dims={"l": 1, "q": []},
        x=[F(1)], y=[], z=[F(1)], s=[F(0)])
    assert rep.stationarity == F(0)
    assert rep.conic_residual == F(0)
    assert rep.comp_slack == F(0)
    assert rep.s_in_cone and rep.z_in_cone
    assert rep.certified(tol=F(0))


def test_negated_cone_dual_is_rejected():
    """Flipping the cone dual out of K* fails membership exactly."""
    rep = kkt.verify_socp(
        c=[F(1)], a=[], b=[], g=[[F(-1)]], h=[F(-1)],
        dims={"l": 1, "q": []},
        x=[F(1)], y=[], z=[F(-1)], s=[F(0)])
    assert not rep.z_in_cone
    assert not rep.certified()


@pytest.mark.slow
def test_untrusted_ecos_socp_certified():
    """An UNTRUSTED ECOS solve (native standard-form x/y/z/s duals) of
    a min-norm-with-thrust-cone SOCP is verified exactly: min t s.t.
    (t, u) in SOC and A u = b. The cone is active (t = ||u||) and the
    conic KKT certifies within tolerance."""
    ecos = pytest.importorskip("ecos")
    from scipy import sparse  # noqa: PLC0415

    # x = [u0,u1,u2, t]; min t
    c = np.array([0.0, 0.0, 0.0, 1.0])
    # equality A u = b  (underdetermined -> nontrivial min-norm)
    a = sparse.csc_matrix(np.array([[1.0, 0.0, 0.0, 0.0],
                                    [0.0, 1.0, 0.0, 0.0]]))
    b = np.array([0.5, -0.3])
    # conic s = (t, u) in SOC^4:  G x + s = h,  h = 0,  s = -G x
    g = sparse.csc_matrix(np.array([
        [0.0, 0.0, 0.0, -1.0],   # s0 = t
        [-1.0, 0.0, 0.0, 0.0],   # s1 = u0
        [0.0, -1.0, 0.0, 0.0],
        [0.0, 0.0, -1.0, 0.0]]))
    h = np.zeros(4)
    dims = {"l": 0, "q": [4]}
    sol = ecos.solve(c, g, h, dims, a, b, verbose=False)
    assert sol["info"]["exitFlag"] == 0

    rep = kkt.verify_socp(
        c=kkt.rationalize_vec(c.tolist()),
        a=kkt.rationalize_mat(a.toarray().tolist()),
        b=kkt.rationalize_vec(b.tolist()),
        g=kkt.rationalize_mat(g.toarray().tolist()),
        h=kkt.rationalize_vec(h.tolist()),
        dims=dims,
        x=kkt.rationalize_vec(sol["x"].tolist()),
        y=kkt.rationalize_vec(sol["y"].tolist()),
        z=kkt.rationalize_vec(sol["z"].tolist()),
        s=kkt.rationalize_vec(sol["s"].tolist()),
        cone_tol=kkt.Frac(1, 10**6))
    # The exact re-check MEASURES that ECOS's solution satisfies the conic
    # KKT residuals to a tight tolerance. certified() itself requires EXACT
    # conic-dual feasibility (stationarity == 0), which a floating-point
    # solve does not provide, so for approximate output this is an exact
    # residual report, not a suboptimality certificate (see SOCPReport).
    tol = kkt.Frac(1, 10**5)
    assert rep.stationarity <= tol
    assert rep.conic_residual <= tol
    assert abs(rep.comp_slack) <= tol
    assert rep.s_in_cone and rep.z_in_cone
    # the optimum is t = ||u|| with u the min-norm reach -> ~0.5831
    assert abs(float(rep.primal_obj) - np.hypot(0.5, 0.3)) < 1e-6


@pytest.mark.slow
def test_embedded_ecos_layer0_socp_certified():
    """The embedded-solver + verified-KKT loop on a real Layer-0
    guidance problem: a min-fuel rendezvous SOCP (per-step thrust cones
    ||u_k|| <= t_k, min sum t_k, CW-STM reach) solved by the EMBEDDED
    ECOS solver and re-verified exactly by certify_ecos."""
    cp = pytest.importorskip("cvxpy")
    pytest.importorskip("ecos")

    nn = math.sqrt(const.MU_EARTH / 6_778_137.0**3)
    ts = [0.0, 250.0, 500.0, 750.0]
    tf = 1000.0
    bmat = np.vstack([np.zeros((3, 3)), np.eye(3)])
    blocks = [cw.stm(nn, tf - t) @ bmat for t in ts]      # 6x3 each
    target = np.array([0.0, -8.0, 0.0, 0.0, 0.0, 0.0])

    u = [cp.Variable(3) for _ in ts]
    tmag = cp.Variable(len(ts))
    cons = [cp.norm(u[k]) <= tmag[k] for k in range(len(ts))]
    cons.append(sum(blocks[k] @ u[k] for k in range(len(ts))) == target)
    prob = cp.Problem(cp.Minimize(cp.sum(tmag)), cons)

    sol, rep = kkt.certify_ecos(prob)
    assert sol["info"]["exitFlag"] == 0
    # exact residual measurement of the embedded solver's output (see the
    # note in test_untrusted_ecos_socp_certified): the conic-KKT residuals
    # are tight, though certified() requires exact conic-dual feasibility.
    tol = kkt.Frac(1, 10**4)
    assert rep.stationarity <= tol
    assert rep.conic_residual <= tol
    assert abs(rep.comp_slack) <= tol
    assert rep.s_in_cone and rep.z_in_cone
    # the exact-measured objective matches the embedded solver's cost
    assert abs(float(rep.primal_obj) - sol["info"]["pcost"]) < 1e-6
    # at least one thrust cone is active (a real burn happened)
    prob.solve(solver=cp.ECOS)
    assert sum(np.linalg.norm(uk.value) for uk in u) > 1e-3


def test_certify_ecos_rejects_tampered_primal():
    """If the returned solution is corrupted, re-verification catches
    it — the certificate is not a rubber stamp. (Done at the checker
    level, since certify_ecos wraps a trustworthy solve.)"""
    # a feasible standard-form point, then a tampered one
    good = kkt.verify_socp(
        c=[F(1)], a=[], b=[], g=[[F(-1)]], h=[F(-1)],
        dims={"l": 1, "q": []}, x=[F(1)], y=[], z=[F(1)], s=[F(0)])
    assert good.certified(tol=F(0))
    tampered = kkt.verify_socp(
        c=[F(1)], a=[], b=[], g=[[F(-1)]], h=[F(-1)],
        dims={"l": 1, "q": []}, x=[F(3)], y=[], z=[F(1)], s=[F(0)])
    assert not tampered.certified()   # slack now inconsistent


def test_indefinite_qp_is_rejected_as_nonconvex():
    """A stationary point of a non-convex (indefinite-P) QP is not a
    global optimum. x=0 satisfies every KKT residual for min -1/2 x^2,
    but the checker must reject it by verifying P >= 0 -- otherwise it
    would 'certify' a maximum as optimal."""
    rep = kkt.verify_qp(p=[[F(-1)]], q=[F(0)], g=[], h=[], a=[], b=[],
                        x=[F(0)], mu=[], nu=[])
    assert rep.stationarity == F(0)                # KKT residuals vanish...
    assert not rep.certified()                     # ...but P is not PSD
    assert any("positive semidefinite" in p for p in rep.problems)


def test_epsilon_stationary_lp_is_not_certified():
    """SOUNDNESS REGRESSION: a small stationarity residual is NOT an
    objective-gap bound when P is rank-deficient. min -eps*x s.t.
    0 <= x <= M with x=0 is eps-stationary but 100 away from the true
    optimum (x=M). The checker must NOT certify it: the residual is not
    in range(P) = {0}, so the Lagrangian dual is unbounded and there is
    no valid suboptimality bound."""
    M = F(10**12)
    eps = F(1, 10**10)
    rep = kkt.verify_qp(
        p=[[F(0)]], q=[-eps],
        g=[[-F(1)], [F(1)]], h=[F(0), M], a=[], b=[],
        x=[F(0)], mu=[F(0), F(0)], nu=[])
    assert rep.stationarity == eps          # under any tol >= 1e-10
    assert rep.duality_gap == F(0)          # complementarity term alone is 0
    assert rep.suboptimality_bound is None  # no valid dual bound
    assert not rep.certified()              # ...so NOT certified
    assert not rep.certified(tol=F(1))      # not even at a loose tolerance


def test_malformed_qp_shape_is_rejected():
    """A too-short G/h/mu silently drops inequality rows from the
    residuals; the checker must reject the malformed instance rather than
    certify it."""
    rep = kkt.verify_qp(
        p=[[F(1)]], q=[F(0)],
        g=[[F(-1)], [F(1)]], h=[F(-1)], a=[], b=[],   # h/mu shorter than G
        x=[F(1)], mu=[F(1)], nu=[])
    assert rep.problems
    assert not rep.certified()


def test_float_input_is_rejected():
    """A float in the trusted path could round a residual to zero; the
    exact-arithmetic contract is enforced, so a float instance is
    rejected rather than accepted by cancellation."""
    rep = kkt.verify_qp(
        p=[[1.0]], q=[0.0], g=[], h=[], a=[], b=[],
        x=[0.0], mu=[], nu=[])
    assert any("Fraction" in m or "inexact" in m for m in rep.problems)
    assert not rep.certified()


def test_socp_uncovered_cone_tail_is_rejected():
    """If dims do not cover all of s/z, a tail entry (here a negative,
    out-of-cone slack) is never membership-checked; the coverage guard
    rejects it."""
    rep = kkt.verify_socp(
        c=[F(1)], a=[], b=[], g=[[F(-1)], [F(0)]], h=[F(-1), F(0)],
        dims={"l": 1, "q": []},                 # covers only 1 of 2 entries
        x=[F(1)], y=[], z=[F(1), F(1)], s=[F(0), F(-5)])
    assert rep.problems
    assert not rep.certified()


def test_suboptimality_bound_is_exact_for_pd_qp():
    """For a positive-definite P the rigorous bound equals the true gap.
    min 1/2 x^2 (optimum x=0, p*=0); the point x=1 has stationarity
    residual r=1 and true suboptimality 1/2, and the checker reports
    exactly 1/2 (= 1/2 r' P^-1 r) rather than the meaningless
    complementarity gap 0."""
    rep = kkt.verify_qp(p=[[F(1)]], q=[F(0)], g=[], h=[], a=[], b=[],
                        x=[F(1)], mu=[], nu=[])
    assert rep.duality_gap == F(0)                    # no constraints
    assert rep.suboptimality_bound == F(1, 2)         # = true gap, exact
    assert not rep.certified(tol=F(1, 10))            # 1/2 > 1/10
    assert rep.certified(tol=F(1))                    # 1/2 <= 1
