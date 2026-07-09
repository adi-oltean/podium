"""Trusted-checker input guards.

The exact-rational verifiers must REFUSE inexact (float) input rather than
certify on it — a CONTRIBUTING non-negotiable (the trusted path is float-free).
`kkt.verify_qp` already enforces this in test_kkt; these lock down the rest.
"""
from fractions import Fraction as F

import pytest

from podium.verify import barrier, bracket, lyapunov, sos


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
