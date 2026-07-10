"""Trusted-checker input guards and refusal paths.

The exact-rational verifiers must REFUSE inexact (float) input and malformed
(wrong-shape) data rather than certify on it — a silent-wrong-answer here would
defeat the point. These lock down the refusal branches that the well-formed
tests elsewhere do not exercise.
"""
from fractions import Fraction as F

import pytest

from podium.verify import barrier, bracket, kkt, lyapunov, scvx_cut, sos

# --- float-input guards (the trusted path is float-free) --------------------

def test_is_psd_refuses_float_input():
    # exact PSD test cannot certify from inexact entries
    assert not barrier.is_psd([[1.0, 0.0], [0.0, 1.0]])


def test_is_sos_refuses_float_input():
    ok, _ = sos.is_sos({(2,): F(1)}, [(1,)], [[1.0]])
    assert not ok


def test_verify_lyapunov_refuses_float_input():
    rep = lyapunov.verify_lyapunov([[1.0, 0.0], [0.0, 1.0]],
                                   [[1.0, 0.0], [0.0, 1.0]])
    assert not rep.certified()


def test_check_qcqp_refuses_float_input():
    with pytest.raises(ValueError):
        bracket._check_qcqp([[1.0]], [F(0)], [[1.0]], [F(0)])


# --- kkt.verify_qp shape validation (malformed data must not certify) -------

def _wellformed_qp() -> dict:
    # min 1/2||x||^2 s.t. x1 + x2 = 2  ->  x=(1,1), nu=-1; residuals exactly 0
    return dict(p=[[F(1), F(0)], [F(0), F(1)]], q=[F(0), F(0)],
                g=[], h=[], a=[[F(1), F(1)]], b=[F(2)],
                x=[F(1), F(1)], mu=[], nu=[F(-1)])


def test_verify_qp_wellformed_baseline_certifies():
    assert kkt.verify_qp(**_wellformed_qp()).certified(tol=F(0))


def test_verify_qp_refuses_wrong_q_length():
    d = _wellformed_qp()
    d["q"] = [F(0)]
    assert not kkt.verify_qp(**d).certified()


def test_verify_qp_refuses_wrong_g_width():
    d = _wellformed_qp()
    d.update(g=[[F(1)]], h=[F(0)], mu=[F(0)])
    assert not kkt.verify_qp(**d).certified()


def test_verify_qp_refuses_mismatched_a_b_nu():
    d = _wellformed_qp()
    d["b"] = [F(2), F(3)]
    assert not kkt.verify_qp(**d).certified()


def test_verify_qp_refuses_wrong_a_width():
    d = _wellformed_qp()
    d.update(a=[[F(1)]], nu=[F(-1)])
    assert not kkt.verify_qp(**d).certified()


def test_verify_qp_refuses_non_square_p():
    d = _wellformed_qp()
    d["p"] = [[F(1)]]
    assert not kkt.verify_qp(**d).certified()


# --- shape mismatches in the other checkers ---------------------------------

def test_check_qcqp_refuses_dim_mismatch():
    with pytest.raises(ValueError):
        bracket._check_qcqp([[F(1), F(0)], [F(0), F(1)]], [F(0)],
                            [[F(1)]], [F(0)])


def test_verify_lyapunov_refuses_shape_mismatch():
    rep = lyapunov.verify_lyapunov([[F(1), F(0)], [F(0), F(1)]], [[F(1)]])
    assert not rep.certified()


# --- multi-constraint bracket: the upper-bound leg decides feasibility with
#     an exact _quad(...) < 0 test, so (unlike the is_psd-guarded lower leg) a
#     float datum could wrong-accept an eps-infeasible point and return an
#     upper bound BELOW J*, collapsing the bracket beneath the true optimum.

def _multi() -> dict:
    # min ||x||^2 s.t. ||x - (2,0)|| >= 3   (a single keep-out, as a list)
    con = ([[F(1), F(0)], [F(0), F(1)]], [F(-4), F(0)], F(-5))
    return dict(p0=[[F(1), F(0)], [F(0), F(1)]], q0=[F(0), F(0)], r0=F(0),
                cons=[con], x=[F(0), F(3)])


def test_certify_upper_bound_multi_baseline():
    d = _multi()
    assert bracket.certify_upper_bound_multi(**d) == F(9)


def test_certify_upper_bound_multi_refuses_float_data():
    d = _multi()
    d["cons"] = [([[1.0, F(0)], [F(0), F(1)]], [F(-4), F(0)], F(-5))]
    with pytest.raises(ValueError):
        bracket.certify_upper_bound_multi(**d)


def test_certify_upper_bound_multi_refuses_float_x():
    d = _multi()
    d["x"] = [0.0, F(3)]            # exact data, but a float candidate point
    with pytest.raises(ValueError):
        bracket.certify_upper_bound_multi(**d)


def test_certify_upper_bound_refuses_float_x():
    # single-constraint upper leg: r0/r1/x also feed _quad and must be exact
    with pytest.raises(ValueError):
        bracket.certify_upper_bound(
            [[F(1), F(0)], [F(0), F(1)]], [F(0), F(0)], F(0),
            [[F(1), F(0)], [F(0), F(1)]], [F(-4), F(0)], F(-5),
            [0.0, F(3)])


# --- scvx_cut.certify_cut: exact-rational Positivstellensatz witness --------

def _cut() -> tuple:
    q, cut, mult, basis, gram = scvx_cut.superquadric_diagonal_certificate(F(1))
    return q, [(mult, cut)], basis, gram


def test_certify_cut_baseline():
    q, cuts, basis, gram = _cut()
    assert scvx_cut.certify_cut(q, cuts, basis, gram).certified


def test_certify_cut_refuses_shape_mismatch():
    q, cuts, basis, gram = _cut()
    with pytest.raises(ValueError):
        scvx_cut.certify_cut(q, cuts, basis, gram[:-1])


def test_certify_cut_refuses_float_data():
    q, cuts, basis, gram = _cut()
    with pytest.raises(ValueError):
        scvx_cut.certify_cut(q, [(1.0, cuts[0][1])], basis, gram)
