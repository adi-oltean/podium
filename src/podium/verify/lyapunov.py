"""Certificate-carrying Lyapunov ellipsoid invariants (credible
autocoding, Feron-style).

An LQR controller synthesizes a value-function matrix P; the ellipsoid
{x : x' P x <= c} is invariant and contracting for the closed loop
A_cl = A - B K, because

    P - A_cl' P A_cl = Q + K' R K  >= 0,

so x' P x is non-increasing along every closed-loop trajectory. This
module treats P as a CERTIFICATE and re-verifies that inequality in
`fractions.Fraction` arithmetic with no floating point in the trusted
path, reusing the exact all-principal-minors PSD check from
podium.verify.barrier. The decrease term Q + K'RK carries a large
margin (Q >= 0 by design), so rationalizing the float Riccati solution
leaves the exact inequality comfortably satisfied — the certificate is
robust, not knife-edge.

The C emitter can render P as a quadratic PROVE/ACSL obligation on the
flight controller, closing the credible-autocoding loop; here the
sandbox check and the exact re-verification stand in for that.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from podium.verify.barrier import Frac, Mat, _det, is_psd

Vec = list[Fraction]


def _matmul(x: Mat, y: Mat) -> Mat:
    n, k, m = len(x), len(y), len(y[0])
    return [[sum((x[i][t] * y[t][j] for t in range(k)), Frac(0))
             for j in range(m)] for i in range(n)]


def _transpose(x: Mat) -> Mat:
    return [[x[i][j] for i in range(len(x))] for j in range(len(x[0]))]


def _sub(x: Mat, y: Mat) -> Mat:
    return [[x[i][j] - y[i][j] for j in range(len(x[0]))]
            for i in range(len(x))]


def rationalize_matrix(m: list[list[float]], max_den: int = 10**12) -> Mat:
    return [[Frac(float(v)).limit_denominator(max_den) for v in row]
            for row in m]


@dataclass
class EllipsoidInvariant:
    """A Lyapunov certificate P (exact rationals) for a closed loop."""

    p: Mat

    def value(self, x: Vec) -> Frac:
        """x' P x (exact)."""
        px = [sum((self.p[i][j] * x[j] for j in range(len(x))), Frac(0))
              for i in range(len(x))]
        return sum((x[i] * px[i] for i in range(len(x))), Frac(0))


@dataclass
class LyapunovReport:
    p_positive: bool          # P > 0 (bounded ellipsoid sub-level set)
    decrease_psd: bool        # P - A_cl' P A_cl >= 0 (non-increasing)
    problems: list[str] = field(default_factory=list)

    def certified(self) -> bool:
        return self.p_positive and self.decrease_psd and not self.problems


def verify_lyapunov(a_cl: Mat, p: Mat) -> LyapunovReport:
    """Exactly verify that P certifies the closed loop A_cl: P > 0
    (positive definite, so {x : x'Px <= c} is a bounded ellipsoid) and
    the Lyapunov decrease P - A_cl' P A_cl >= 0, by the exact
    all-principal-minors test. P > 0 is checked as PSD plus nonsingular
    (a PSD matrix with nonzero determinant has all-positive eigenvalues);
    without strictness P = 0 would spuriously certify. Inputs are
    Fraction matrices (rationalize float syntheses first)."""
    problems: list[str] = []
    n = len(p)
    if (any(len(row) != n for row in p) or len(a_cl) != n
            or any(len(row) != n for row in a_cl)):
        problems.append("P and A_cl must be square and same size")
        return LyapunovReport(False, False, problems)
    # symmetry of P (a Lyapunov matrix is symmetric)
    for i in range(n):
        for j in range(i):
            if p[i][j] != p[j][i]:
                problems.append(f"P not symmetric at ({i},{j})")
                break
    p_pos = is_psd(p) and _det(p) != 0        # positive definite (P > 0)
    decrease = _sub(p, _matmul(_matmul(_transpose(a_cl), p), a_cl))
    # symmetrize exactly (A_cl'PA_cl is symmetric in exact arithmetic
    # when P is; guard against asymmetry from a non-symmetric P input)
    dec_sym = [[(decrease[i][j] + decrease[j][i]) / Frac(2) for j in range(n)]
               for i in range(n)]
    dec_psd = is_psd(dec_sym)
    return LyapunovReport(p_positive=p_pos, decrease_psd=dec_psd,
                          problems=problems)
