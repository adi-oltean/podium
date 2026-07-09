"""Trusted-checker input guards and refusal paths.

The exact-rational verifiers must REFUSE inexact (float) input and malformed
(wrong-shape) data rather than certify on it — a silent-wrong-answer here would
defeat the point. These lock down the refusal branches that the well-formed
tests elsewhere do not exercise.
"""
from fractions import Fraction as F

import pytest

from podium.verify import barrier, bracket, kkt, lyapunov, sos

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
