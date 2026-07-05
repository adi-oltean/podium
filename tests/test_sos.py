"""Exact-rational SOS / Positivstellensatz receipts (#51): the
quadratic S-procedure barrier (#20) generalized to higher-degree
polynomial systems, verified exactly.

The flagship is a genuinely nonlinear (cubic) Duffing oscillator with a
QUARTIC Lyapunov barrier whose Lie derivative -dV/dt is certified
sum-of-squares over the rationals, proving the sub-level set is an
infinite-horizon invariant --- the higher-degree analogue of the
quadratic abort-safety barrier.
"""

from fractions import Fraction as F

import numpy as np

from podium.verify import sos


def test_basic_sos_diagonal_and_coupled():
    """x1^4 + x2^4 is SOS (diagonal Gram); (x1^2+x2^2)^2 is SOS with a
    rank-1 non-diagonal Gram --- both certified by exact identity + PSD."""
    basis = [(2, 0), (0, 2)]          # z = [x1^2, x2^2]
    quartic_sum = {(4, 0): F(1), (0, 4): F(1)}
    ok, prob = sos.is_sos(quartic_sum, basis, [[F(1), F(0)], [F(0), F(1)]])
    assert ok, prob
    perfect_square = {(4, 0): F(1), (2, 2): F(2), (0, 4): F(1)}
    ok2, prob2 = sos.is_sos(perfect_square, basis,
                            [[F(1), F(1)], [F(1), F(1)]])
    assert ok2, prob2


def test_indefinite_polynomial_is_rejected():
    """x1^2 - x2^2 is not SOS: the only Gram reproducing it is
    indefinite, and the exact PSD check catches it."""
    basis = [(1, 0), (0, 1)]
    indef = {(2, 0): F(1), (0, 2): F(-1)}
    ok, prob = sos.is_sos(indef, basis, [[F(1), F(0)], [F(0), F(-1)]])
    assert not ok
    assert any("semidefinite" in m for m in prob)
    # a Gram that does not reproduce the polynomial is caught too
    ok2, prob2 = sos.is_sos(indef, basis, [[F(1), F(0)], [F(0), F(1)]])
    assert not ok2
    assert any("identity" in m for m in prob2)


def test_lie_derivative_is_computed_exactly():
    """Duffing: xdot1 = x2, xdot2 = -x1 - x1^3 - x2, with the quartic
    energy V = 1/2 x1^2 + 1/4 x1^4 + 1/2 x2^2 gives dV/dt = -x2^2
    EXACTLY (a nontrivial cancellation of the cubic cross terms)."""
    v = {(2, 0): F(1, 2), (4, 0): F(1, 4), (0, 2): F(1, 2)}
    f = [{(0, 1): F(1)},
         {(1, 0): F(-1), (3, 0): F(-1), (0, 1): F(-1)}]
    vdot = sos.lie_derivative(v, f)
    assert vdot == {(0, 2): F(-1)}       # -x2^2, cross terms cancel


def test_duffing_quartic_barrier_lie_derivative_is_sos():
    """The infinite-horizon nonlinear barrier: -dV/dt = x2^2 is SOS
    (basis [x2], Gram [[1]]), so V is non-increasing along the flow and
    every sub-level set {V <= c} is invariant --- for a genuinely cubic
    system with a quartic certificate, verified exactly."""
    v = {(2, 0): F(1, 2), (4, 0): F(1, 4), (0, 2): F(1, 2)}
    f = [{(0, 1): F(1)},
         {(1, 0): F(-1), (3, 0): F(-1), (0, 1): F(-1)}]
    neg_vdot = sos.pscale(F(-1), sos.lie_derivative(v, f))
    ok, prob = sos.is_sos(neg_vdot, [(0, 1)], [[F(1)]])
    assert ok, prob


def test_duffing_sublevel_set_invariant_in_simulation():
    """Numerically confirm the certified invariant: the quartic V is
    monotone non-increasing along the cubic Duffing flow, so a start
    inside {V <= c} never leaves."""
    def field(x):
        return np.array([x[1], -x[0] - x[0]**3 - x[1]])

    def energy(x):
        return 0.5 * x[0]**2 + 0.25 * x[0]**4 + 0.5 * x[1]**2

    x = np.array([1.2, 0.8])
    c = energy(x)
    dt, v_prev = 0.002, energy(x)
    for _ in range(20000):
        k1 = field(x)
        k2 = field(x + 0.5 * dt * k1)
        k3 = field(x + 0.5 * dt * k2)
        k4 = field(x + dt * k3)
        x = x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        v = energy(x)
        assert v <= v_prev + 1e-9
        assert v <= c + 1e-9
        v_prev = v
