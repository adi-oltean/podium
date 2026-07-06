"""Exact-arithmetic optimality-gap certificates for nonconvex QCQPs.

For a nonconvex quadratically-constrained quadratic program

    J* = min  x'P0 x + q0'x + r0    s.t.  x'P1 x + q1'x + r1 >= 0

(the constraint concave, e.g. a spherical keep-out ||x-c|| >= R, so the
feasible set is nonconvex) this module brackets the true optimum J* by
two certificates, each checked in exact rationals with no floating point
in the trusted path -- the same "check the answer, not the run"
discipline as the barrier / KKT / Lyapunov / SOS checkers:

* Lower bound (S-procedure LMI dual). t <= J* is certified iff there is
  a multiplier lambda >= 0 with the (n+1)x(n+1) matrix

      M(lambda, t) = [[ P0 - lambda P1,        (q0 - lambda q1)/2 ],
                      [ (q0 - lambda q1)'/2,    r0 - lambda r1 - t ]]

  positive semidefinite. Then t <= J* by weak duality (any dual-feasible
  point lower-bounds the primal). `certify_lower_bound` checks
  lambda >= 0 and is_psd(M) exactly -- no matrix inverse.

* Upper bound. Any point feasible for the true nonconvex constraint gives
  J* <= its objective. `podium.verify.scvx_cut` produces such a point as
  the solution of the convex program over a certified SOUND half-space
  inner-approximation of the keep-out, and `podium.verify.kkt` certifies
  that solution exactly.

When the two meet (J_lb == J_ub) the bracket closes: an exact certificate
of the GLOBAL optimum of a nonconvex problem.
"""

from __future__ import annotations

from fractions import Fraction as F

from podium.verify.barrier import is_psd

Vec = list[F]
Mat = list[list[F]]


def keepout_qcqp(center: tuple[F, ...], radius: F
                 ) -> tuple[Mat, Vec, F, Mat, Vec, F]:
    """QCQP data for  min ||x||^2  s.t.  ||x - center|| >= radius.
    Returns (P0, q0, r0, P1, q1, r1) with the keep-out written as
    x'P1 x + q1'x + r1 >= 0."""
    n = len(center)
    ident: Mat = [[F(1) if i == j else F(0) for j in range(n)]
                  for i in range(n)]
    p0, q0, r0 = ident, [F(0)] * n, F(0)
    p1 = ident
    q1 = [-2 * center[i] for i in range(n)]
    r1 = sum((center[i] ** 2 for i in range(n)), F(0)) - radius ** 2
    return p0, q0, r0, p1, q1, r1


def lower_bound_matrix(p0: Mat, q0: Vec, r0: F, p1: Mat, q1: Vec, r1: F,
                       lam: F, t: F) -> Mat:
    """Assemble the S-procedure LMI M(lambda, t) (see module docstring)."""
    n = len(q0)
    m: Mat = [[F(0)] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(n):
            m[i][j] = p0[i][j] - lam * p1[i][j]
        b = (q0[i] - lam * q1[i]) / 2
        m[i][n] = b
        m[n][i] = b
    m[n][n] = r0 - lam * r1 - t
    return m


def certify_lower_bound(p0: Mat, q0: Vec, r0: F, p1: Mat, q1: Vec, r1: F,
                        lam: F, t: F) -> bool:
    """Exact certificate that t <= J* for the nonconvex QCQP: returns True
    iff lam >= 0 and M(lam, t) >= 0. All inputs Fractions (rationalize an
    SDP-dual solution first)."""
    return lam >= 0 and is_psd(lower_bound_matrix(
        p0, q0, r0, p1, q1, r1, lam, t))
