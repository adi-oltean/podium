"""Verified-KKT checker: an exact-rational optimality certificate for a
convex quadratic program.

The online convex solvers (Clarabel, ECOS, an embedded QOCOGEN kernel)
are UNTRUSTED. This checker re-verifies a claimed solution's optimality
in `fractions.Fraction` arithmetic, so the certificate carries no
floating-point uncertainty: every residual reported is the EXACT KKT
violation of the given rationalized point.

For a convex QP

    minimize   1/2 x' P x + q' x
    subject to G x <= h,   A x = b            (P symmetric PSD)

the KKT conditions at a primal x with duals (mu >= 0 for G, nu for A)
are stationarity  P x + q + G' mu + A' nu = 0, primal feasibility
G x <= h and A x = b, dual feasibility mu >= 0, and complementary
slackness mu_i (h - G x)_i = 0. For a CONVEX QP these are necessary AND
sufficient for global optimality, and when x is primal-feasible and
(mu, nu) dual-feasible the duality gap is exactly sum_i mu_i (h-Gx)_i,
which bounds the suboptimality p(x) - p*. `verify_qp` returns all of
these as exact Fractions.

Trusted-checker discipline (as in podium.verify.barrier): only Fraction
arithmetic, no floats in the verification path; `rationalize_*` convert
a solver's float output to exact rationals before checking.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

Frac = Fraction
Vec = list[Fraction]
Mat = list[list[Fraction]]


def _mv(m: Mat, v: Vec) -> Vec:
    """Matrix-vector product (exact)."""
    return [sum((row[j] * v[j] for j in range(len(v))), Frac(0))
            for row in m]


def _mtv(m: Mat, v: Vec) -> Vec:
    """Transpose-matrix-vector product M' v (exact)."""
    if not m:
        return []
    n = len(m[0])
    return [sum((m[i][j] * v[i] for i in range(len(m))), Frac(0))
            for j in range(n)]


def _rat(x: float, max_den: int = 10**12) -> Frac:
    return Frac(x).limit_denominator(max_den)


def rationalize_vec(v: list[float], max_den: int = 10**12) -> Vec:
    return [_rat(float(x), max_den) for x in v]


def rationalize_mat(m: list[list[float]], max_den: int = 10**12) -> Mat:
    return [[_rat(float(x), max_den) for x in row] for row in m]


@dataclass
class KKTReport:
    """Exact KKT residuals of a claimed QP solution. All Fractions."""

    stationarity: Frac      # max |P x + q + G' mu + A' nu|
    eq_residual: Frac       # max |A x - b|
    ineq_violation: Frac    # max(0, (G x - h)_i)  — primal infeasibility
    dual_violation: Frac    # max(0, -mu_i)        — dual infeasibility
    duality_gap: Frac       # sum_i mu_i (h - G x)_i
    primal_obj: Frac        # 1/2 x' P x + q' x
    problems: list[str] = field(default_factory=list)

    def certified(self, tol: Frac = Frac(1, 10**9)) -> bool:
        """True iff every exact residual is within tol (and no
        structural problem was found). For a convex QP this certifies
        x is within `duality_gap` of the global optimum."""
        return (not self.problems
                and self.stationarity <= tol
                and self.eq_residual <= tol
                and self.ineq_violation <= tol
                and self.dual_violation <= tol
                and abs(self.duality_gap) <= tol)


def verify_qp(
    p: Mat, q: Vec, g: Mat, h: Vec, a: Mat, b: Vec,
    x: Vec, mu: Vec, nu: Vec,
) -> KKTReport:
    """Exact KKT verification of (x, mu, nu) for the convex QP. All
    inputs must be Fractions (use rationalize_* on solver output).
    G/h and A/b may be empty ([]) for an unconstrained/equality-free
    or inequality-free problem."""
    problems: list[str] = []
    n = len(x)
    if len(p) != n or any(len(row) != n for row in p):
        problems.append("P must be n x n")
    # symmetry (exact) — a PSD-objective QP must have symmetric P
    for i in range(len(p)):
        for j in range(i):
            if p[i][j] != p[j][i]:
                problems.append(f"P not symmetric at ({i},{j})")
                break

    def _absmax(v: Vec) -> Frac:
        return max((abs(t) for t in v), default=Frac(0))

    # stationarity: P x + q + G' mu + A' nu
    stat = [sum((p[i][j] * x[j] for j in range(n)), Frac(0)) + q[i]
            for i in range(n)]
    gtmu = _mtv(g, mu) if g else [Frac(0)] * n
    atnu = _mtv(a, nu) if a else [Frac(0)] * n
    stat = [stat[i] + gtmu[i] + atnu[i] for i in range(n)]

    # inequality slack s = h - G x  (>= 0 for feasibility)
    if g:
        gx = _mv(g, x)
        slack = [h[i] - gx[i] for i in range(len(h))]
        ineq_viol = max((-s for s in slack), default=Frac(0))
        ineq_viol = max(ineq_viol, Frac(0))
        gap = sum((mu[i] * slack[i] for i in range(len(slack))), Frac(0))
        dual_viol = max((-m for m in mu), default=Frac(0))
        dual_viol = max(dual_viol, Frac(0))
    else:
        ineq_viol = Frac(0)
        gap = Frac(0)
        dual_viol = Frac(0)

    # equality residual A x - b
    if a:
        ax = _mv(a, x)
        eq_res = _absmax([ax[i] - b[i] for i in range(len(b))])
    else:
        eq_res = Frac(0)

    # primal objective 1/2 x' P x + q' x
    px = _mv(p, x) if p else [Frac(0)] * n
    obj = (sum((x[i] * px[i] for i in range(n)), Frac(0)) / 2
           + sum((q[i] * x[i] for i in range(n)), Frac(0)))

    return KKTReport(
        stationarity=_absmax(stat),
        eq_residual=eq_res,
        ineq_violation=ineq_viol,
        dual_violation=dual_viol,
        duality_gap=gap,
        primal_obj=obj,
        problems=problems,
    )
