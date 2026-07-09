"""Refusal-branch coverage for the trusted conic-KKT checker
(podium.verify.kkt).

These exercise the SOCP shape/cone-coverage guards, the exact
second-order-cone membership on genuine SOC blocks (not just the
nonneg-cone LP forms of tests/test_kkt.py), and the embedded-ECOS
re-verification wrapper. Each malformed instance must be REFUSED
(the specific guard fires and .certified() is False), never certified.
"""

from fractions import Fraction as F

import pytest

from podium.verify import kkt

# --- helper contract ---------------------------------------------------

def test_mtv_of_empty_matrix_is_empty():
    """M' v for an empty (no-row) M is the empty vector: the transpose
    product has no columns to accumulate. Guards the g/a-absent SOCP/QP
    paths that short-circuit to a zero stationarity contribution."""
    assert kkt._mtv([], [F(1), F(2)]) == []


# --- SOCP shape / cone-coverage refusals -------------------------------

def test_socp_float_input_is_rejected():
    """A single float in the trusted path is refused: an inexact entry
    could cancel a residual to a spurious zero."""
    rep = kkt.verify_socp(
        c=[1.0], a=[], b=[], g=[[F(-1)]], h=[F(-1)],
        dims={"l": 1, "q": []}, x=[F(1)], y=[], z=[F(1)], s=[F(0)])
    assert any("Fraction" in m or "inexact" in m for m in rep.problems)
    assert not rep.certified()


def test_socp_zero_size_second_order_cone_is_rejected():
    """A declared second-order cone must have size >= 1; a size-0 cone is
    malformed and refused."""
    rep = kkt.verify_socp(
        c=[F(1)], a=[], b=[], g=[], h=[],
        dims={"l": 0, "q": [0]}, x=[F(1)], y=[], z=[], s=[])
    assert any("size >= 1" in m for m in rep.problems)
    assert not rep.certified()


def test_socp_wrong_length_c_is_rejected():
    """len(c) must equal the number of primal variables."""
    rep = kkt.verify_socp(
        c=[F(1), F(2)], a=[], b=[], g=[], h=[],
        dims={"l": 0, "q": []}, x=[F(1)], y=[], z=[], s=[])
    assert any("len(c) != n" in m for m in rep.problems)
    assert not rep.certified()


def test_socp_equality_block_shape_mismatch_is_rejected():
    """A row of A whose width != n is refused (would misalign A x - b)."""
    rep = kkt.verify_socp(
        c=[F(1)], a=[[F(1), F(2)]], b=[F(0)], g=[], h=[],
        dims={"l": 0, "q": []}, x=[F(1)], y=[F(0)], z=[], s=[])
    assert any("A, b shape mismatch" in m for m in rep.problems)
    assert not rep.certified()


def test_socp_dual_y_length_mismatch_is_rejected():
    """len(y) must equal rows(A): a short y silently drops an equality
    multiplier from the stationarity check."""
    rep = kkt.verify_socp(
        c=[F(1)], a=[[F(1)]], b=[F(0)], g=[], h=[],
        dims={"l": 0, "q": []}, x=[F(1)], y=[], z=[], s=[])
    assert any("len(y)" in m for m in rep.problems)
    assert not rep.certified()


def test_socp_inequality_block_shape_mismatch_is_rejected():
    """A row of G whose width != n is refused (would misalign G x + s)."""
    rep = kkt.verify_socp(
        c=[F(1)], a=[], b=[], g=[[F(1), F(2)]], h=[F(0)],
        dims={"l": 1, "q": []}, x=[F(1)], y=[], z=[F(1)], s=[F(0)])
    assert any("G, h shape mismatch" in m for m in rep.problems)
    assert not rep.certified()


# --- genuine second-order-cone membership ------------------------------

def test_socp_second_order_cone_out_of_cone_point_is_rejected():
    """A real SOC block (q=[3]) whose slack and dual are OUTSIDE the
    second-order cone is refused. Here stationarity is exactly zero yet
    z=(0,1,0) has ||(1,0)|| = 1 > 0 = z0, so z is not in K*: no
    suboptimality bound is issued and the point is not certified. This
    drives the SOC branches of _cone_blocks, the tolerance membership
    test, and the exact (tol-0) dual-cone test."""
    identity_neg = [[F(-1), F(0), F(0)],
                    [F(0), F(-1), F(0)],
                    [F(0), F(0), F(-1)]]
    # G = -I, h = 0  =>  s = -G x = x,  stationarity = c + G' z = c - z
    z = [F(0), F(1), F(0)]          # out of SOC: z0=0 < ||(1,0)|| = 1
    rep = kkt.verify_socp(
        c=list(z), a=[], b=[], g=identity_neg, h=[F(0), F(0), F(0)],
        dims={"l": 0, "q": [3]},
        x=[F(0), F(1), F(0)],       # s = x = (0,1,0): also out of SOC
        y=[], z=z, s=[F(0), F(1), F(0)])
    assert rep.stationarity == F(0)          # exact conic-dual stationarity
    assert not rep.s_in_cone                 # SOC slack fails membership
    assert not rep.z_in_cone                 # SOC dual fails membership
    assert rep.suboptimality_bound is None   # z not EXACTLY in K* -> no bound
    assert not rep.certified()


def test_socp_second_order_cone_feasible_point_certifies():
    """Positive control: an exactly conic-feasible SOC point (interior
    slack, dual on the cone boundary, exact complementarity) certifies.
    min t s.t. (t, u) in SOC^3, u = (0,0) fixed by A -> optimum t=0."""
    identity_neg = [[F(-1), F(0), F(0)],
                    [F(0), F(-1), F(0)],
                    [F(0), F(0), F(-1)]]
    # x = (t, u0, u1); c = (1,0,0); A pins u = 0; SOC on s = (t, u0, u1)
    rep = kkt.verify_socp(
        c=[F(1), F(0), F(0)],
        a=[[F(0), F(1), F(0)], [F(0), F(0), F(1)]], b=[F(0), F(0)],
        g=identity_neg, h=[F(0), F(0), F(0)],
        dims={"l": 0, "q": [3]},
        x=[F(0), F(0), F(0)],       # t = 0, u = 0
        y=[F(0), F(0)],
        z=[F(1), F(0), F(0)],       # dual on SOC boundary, exactly in K*
        s=[F(0), F(0), F(0)])
    assert rep.stationarity == F(0)
    assert rep.conic_residual == F(0)
    assert rep.comp_slack == F(0)
    assert rep.s_in_cone and rep.z_in_cone
    assert rep.suboptimality_bound == F(0)   # exact dual bound issued
    assert rep.certified(tol=F(0))


# --- embedded-ECOS re-verification wrapper -----------------------------

def test_certify_ecos_embedded_socp_is_reverified():
    """certify_ecos solves a min-norm SOCP with the EMBEDDED ECOS solver
    and re-checks its x/y/z/s exactly. The conic-KKT residuals of the
    floating-point solve are tight, so the exact re-check MEASURES a valid
    solution (certified() itself needs exact conic-dual feasibility)."""
    cp = pytest.importorskip("cvxpy")
    pytest.importorskip("ecos")
    import numpy as np  # noqa: PLC0415

    u = cp.Variable(2)
    t = cp.Variable()
    prob = cp.Problem(cp.Minimize(t),
                      [cp.norm(u) <= t, u[0] == 0.5, u[1] == -0.3])
    sol, rep = kkt.certify_ecos(prob)
    assert sol["info"]["exitFlag"] == 0
    tol = kkt.Frac(1, 10**5)
    assert rep.stationarity <= tol
    assert rep.conic_residual <= tol
    assert abs(rep.comp_slack) <= tol
    assert rep.s_in_cone and rep.z_in_cone
    # optimum t = ||u|| = hypot(0.5, 0.3)
    assert abs(float(rep.primal_obj) - float(np.hypot(0.5, 0.3))) < 1e-6
