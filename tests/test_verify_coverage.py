"""Refusal-branch and helper coverage for the trusted exact-rational
checkers: bracket (nonconvex-QCQP bracket), sos (Positivstellensatz),
barrier (abort-safety), lyapunov (ellipsoid invariant), and contracts.

Every case asserts CORRECT behavior of a verification-critical branch:
malformed shapes raise, dual-infeasible / non-PSD / singular inputs are
REFUSED, and the exact helpers compute the right value. No source or
existing test is modified.
"""

from fractions import Fraction as F

import pytest

from podium.verify import barrier, bracket, contracts, lyapunov, sos

# ======================================================================
# bracket.py
# ======================================================================

def test_certify_upper_bound_wrong_x_dimension_raises():
    """A point whose dimension disagrees with the QCQP is a caller bug,
    not a silent truncation -- certify_upper_bound raises."""
    p0, q0, r0, p1, q1, r1 = bracket.keepout_qcqp((F(0), F(0)), F(1))
    with pytest.raises(ValueError, match="wrong dimension"):
        bracket.certify_upper_bound(p0, q0, r0, p1, q1, r1, x=[F(0)])


def test_multi_objective_dimension_mismatch_raises():
    """_check_multi rejects an objective whose P0 is not n x n rather than
    indexing past the end later."""
    with pytest.raises(ValueError, match="objective dimensions"):
        bracket.certify_lower_bound_multi(
            p0=[[F(1)]], q0=[F(0), F(0)], r0=F(0),
            cons=[], lams=[], t=F(0))


def test_multi_constraint_dimension_mismatch_raises():
    """A constraint block of the wrong width is rejected."""
    con = ([[F(1), F(2)]], [F(0)], F(0))     # 1x2 P_k with n = 1
    with pytest.raises(ValueError, match="constraint dimensions"):
        bracket.certify_lower_bound_multi(
            p0=[[F(1)]], q0=[F(0)], r0=F(0),
            cons=[con], lams=[F(0)], t=F(0))


def test_multi_lambda_count_mismatch_raises():
    """One multiplier per constraint is required; a mismatched count is a
    malformed dual and is rejected."""
    con = ([[F(1)]], [F(0)], F(0))
    with pytest.raises(ValueError, match="len\\(lams\\)"):
        bracket.certify_lower_bound_multi(
            p0=[[F(1)]], q0=[F(0)], r0=F(0),
            cons=[con], lams=[], t=F(0))


def test_multi_negative_lambda_is_not_certified():
    """A negative S-procedure multiplier is dual-INFEASIBLE: the weak-
    duality proof needs lam >= 0, so no lower bound is certified."""
    con = ([[F(1)]], [F(0)], F(0))
    assert not bracket.certify_lower_bound_multi(
        p0=[[F(1)]], q0=[F(0)], r0=F(0),
        cons=[con], lams=[F(-1)], t=F(0))


def test_multi_upper_bound_wrong_x_dimension_raises():
    con = ([[F(1)]], [F(0)], F(0))
    with pytest.raises(ValueError, match="wrong dimension"):
        bracket.certify_upper_bound_multi(
            p0=[[F(1)]], q0=[F(0)], r0=F(0), cons=[con], x=[F(0), F(0)])


def test_exact_solve_refuses_singular_system():
    """The exact linear solve returns None on a singular matrix (here two
    dependent rows) rather than a spurious solution."""
    assert bracket._solve([[F(1), F(2)], [F(2), F(4)]], [F(1), F(1)]) is None


def test_dual_value_nondiagonal_pd_matrix_is_exact():
    """g(lam) for a non-diagonal positive-definite A = P0 - lam*P1 is the
    exact Lagrangian dual value; exercises the elimination path of the
    exact solve. A = [[2,1],[1,2]] (eigenvalues 1, 3 > 0)."""
    p0 = [[F(2), F(1)], [F(1), F(2)]]
    q0 = [F(2), F(0)]
    p1 = [[F(1), F(0)], [F(0), F(1)]]
    q1 = [F(0), F(0)]
    g = bracket.dual_value(p0, q0, F(0), p1, q1, F(0), lam=F(0))
    # min_x x'A x + q0'x = -1/4 q0' A^{-1} q0, A=[[2,1],[1,2]], A^{-1}=
    # [[2,-1],[-1,2]]/3; q0=[2,0] -> q0'A^{-1}q0 = 8/3 -> g = -2/3
    assert g == F(-2, 3)


def test_dual_value_hard_case_indefinite_A_returns_none():
    """When A = P0 - lam P1 is not positive definite (the 'hard case'),
    the closed-form dual value is undefined and None is returned."""
    p0, q0, r0, p1, q1, r1 = bracket.keepout_qcqp((F(0), F(0)), F(1))
    # lam = 2 makes A = (1 - 2) I = -I, indefinite
    assert bracket.dual_value(p0, q0, r0, p1, q1, r1, lam=F(2)) is None


def test_recover_lower_bound_negative_multiplier_returns_none():
    """A rounded multiplier below zero is dual-infeasible -> no bound."""
    p0, q0, r0, p1, q1, r1 = bracket.keepout_qcqp((F(0), F(0)), F(1))
    assert bracket.recover_lower_bound(
        p0, q0, r0, p1, q1, r1, lam_float=-1.0) is None


def test_recover_lower_bound_infeasible_multiplier_returns_none():
    """A multiplier that lands in the hard case (A not PD) yields no exact
    dual value, so recover_lower_bound returns None."""
    p0, q0, r0, p1, q1, r1 = bracket.keepout_qcqp((F(0), F(0)), F(1))
    assert bracket.recover_lower_bound(
        p0, q0, r0, p1, q1, r1, lam_float=2.0) is None


# ======================================================================
# sos.py
# ======================================================================

def test_mono_builds_exponent_tuple():
    assert sos._mono(2, 0, 1) == (2, 0, 1)


def test_is_sos_rejects_nonsquare_gram():
    """A Gram whose size does not match the basis cannot encode z^T G z;
    it is refused before any PSD test."""
    ok, problems = sos.is_sos({}, basis=[(1, 0), (0, 1)], gram=[[F(1)]])
    assert not ok
    assert any("square" in m for m in problems)


def test_is_sos_rejects_nonsymmetric_gram():
    """A non-symmetric Gram is not a valid quadratic form and is refused."""
    ok, problems = sos.is_sos(
        {}, basis=[(1, 0), (0, 1)],
        gram=[[F(1), F(2)], [F(3), F(1)]])
    assert not ok
    assert any("symmetric" in m for m in problems)


def test_validate_gram_wrong_float_shape_returns_none():
    """A float Gram whose shape disagrees with the basis is rejected."""
    assert sos.validate_gram(
        target={(2, 0): F(1)}, basis=[(1, 0), (0, 1)],
        gram_float=[[1.0]]) is None


def test_validate_gram_indefinite_correction_is_refused():
    """Correcting toward an off-diagonal-only target (2*x*y) forces an
    indefinite Gram [[0,1],[1,0]]; validate_gram must return None rather
    than ship a non-PSD (hence non-SOS) certificate. Exercises the
    off-diagonal residual-correction branch."""
    assert sos.validate_gram(
        target={(1, 1): F(2)}, basis=[(1, 0), (0, 1)],
        gram_float=[[0.0, 0.0], [0.0, 0.0]]) is None


# ======================================================================
# barrier.py
# ======================================================================

def _valid_case() -> barrier.AbortSafetyCase:
    return barrier.AbortSafetyCase(
        center=tuple(F(0) for _ in range(6)),
        radii=tuple(F(1) for _ in range(6)),
        koz_radius=F(1))


def test_det_requires_pivot_row_swap():
    """Exact determinant with a zero leading pivot forces a row swap and a
    sign flip: det([[0,1],[1,0]]) = -1."""
    assert barrier._det([[F(0), F(1)], [F(1), F(0)]]) == F(-1)


def test_barrier_nonpositive_epsilon_is_rejected():
    """The S-procedure margins eps0/eps_u must be strictly positive; a zero
    margin gives no separation and is refused."""
    cert = barrier.BarrierCertificate(
        a=tuple(F(0) for _ in range(5)), lam0=F(0), lam_u=F(0),
        eps0=F(0), eps_u=F(1, 10), case=_valid_case())
    problems = barrier.verify_certificate(cert)
    assert any("eps0 must be a positive Fraction" in m for m in problems)


def test_barrier_float_coefficient_is_rejected():
    """A non-Fraction barrier coefficient breaks the exact-arithmetic
    contract and is refused."""
    cert = barrier.BarrierCertificate(
        a=(1.0, F(0), F(0), F(0), F(0)), lam0=F(0), lam_u=F(0),
        eps0=F(1, 10), eps_u=F(1, 10), case=_valid_case())
    problems = barrier.verify_certificate(cert)
    assert any("coefficients must be Fractions" in m for m in problems)


def test_barrier_malformed_case_center_is_rejected():
    """A case whose center is not length-6 is malformed and refused."""
    bad_case = barrier.AbortSafetyCase(
        center=tuple(F(0) for _ in range(5)),
        radii=tuple(F(1) for _ in range(6)), koz_radius=F(1))
    cert = barrier.BarrierCertificate(
        a=tuple(F(0) for _ in range(5)), lam0=F(0), lam_u=F(0),
        eps0=F(1, 10), eps_u=F(1, 10), case=bad_case)
    problems = barrier.verify_certificate(cert)
    assert any("center/radii" in m for m in problems)


def test_barrier_nonpositive_koz_radius_is_rejected():
    """The keep-out radius must be strictly positive."""
    bad_case = barrier.AbortSafetyCase(
        center=tuple(F(0) for _ in range(6)),
        radii=tuple(F(1) for _ in range(6)), koz_radius=F(0))
    cert = barrier.BarrierCertificate(
        a=tuple(F(0) for _ in range(5)), lam0=F(0), lam_u=F(0),
        eps0=F(1, 10), eps_u=F(1, 10), case=bad_case)
    problems = barrier.verify_certificate(cert)
    assert any("koz_radius must be a positive" in m for m in problems)


def test_barrier_trivial_certificate_fails_safety_psd():
    """A well-formed but vacuous certificate (all coefficients zero) is a
    flow invariant (B = 0, so dB/dtau = 0) yet cannot establish the C1
    safety inequality B <= -eps0 on X0: C1 = -eps0*I is not PSD, so the
    checker reports the failure instead of certifying nothing as safe."""
    cert = barrier.BarrierCertificate(
        a=tuple(F(0) for _ in range(5)), lam0=F(0), lam_u=F(0),
        eps0=F(1, 10), eps_u=F(1, 10), case=_valid_case())
    problems = barrier.verify_certificate(cert)
    assert any("C1 not PSD" in m for m in problems)
    # the conservation check itself passes (no flow-invariant complaint)
    assert not any("flow invariant" in m for m in problems)


# ======================================================================
# lyapunov.py
# ======================================================================

def test_lyapunov_nonsymmetric_p_is_rejected():
    """A Lyapunov matrix must be symmetric; an asymmetric P is flagged and
    the certificate is refused."""
    rep = lyapunov.verify_lyapunov(
        a_cl=[[F(0), F(0)], [F(0), F(0)]],
        p=[[F(1), F(2)], [F(0), F(1)]])
    assert any("not symmetric" in m for m in rep.problems)
    assert not rep.certified()


# ======================================================================
# contracts.py
# ======================================================================

def test_interval_contains_and_annotation():
    iv = contracts.Interval(0.0, 1.0)
    assert iv.contains(0.5)
    assert not iv.contains(1.5)
    assert not iv.contains([0.5, 2.0])       # any element outside -> False
    assert iv.to_annotation() == "range(0.0,1.0)"


def test_contract_on_unknown_argument_raises():
    with pytest.raises(TypeError, match="unknown args"):
        @contracts.contract(nonesuch=contracts.Interval(0.0, 1.0))
        def _f(x):
            return x


def test_contract_enforces_argument_range():
    """An in-range call passes through; an out-of-range argument raises
    ContractError (enforcement is on by default in the sandbox)."""
    @contracts.contract(x=contracts.Interval(0.0, 1.0))
    def scaled(x):
        return x * 2

    assert scaled(0.5) == 1.0
    with pytest.raises(contracts.ContractError, match="outside"):
        scaled(5.0)


def test_shapes_decorator_records_metadata():
    @contracts.shapes(a=(3, 3))
    def kernel(a):
        return a

    assert kernel.__podium_shapes__ == {"a": (3, 3)}


def test_prove_passes_and_raises():
    contracts.prove(True)                    # satisfied invariant: no raise
    with pytest.raises(contracts.ContractError, match="must hold"):
        contracts.prove(False, "must hold")
