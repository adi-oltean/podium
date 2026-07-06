"""Certified sound convexification of a spherical keep-out for SCvx.

Successive convexification linearizes the nonconvex keep-out
``||r - c|| >= R`` into a half-space cut ``n . (r - c) >= R`` at the
reference. With a slack this is soft (iterates may violate); applied
HARD with a rational-unit normal it is a *sound convex
inner-approximation* -- every point of the half-space genuinely lies
outside the ball -- so every SCvx iterate is feasible for the true
nonconvex problem, not just the limit.

This module makes that soundness a machine-checked, exact-rational
certificate (the same discipline as barrier/KKT/Lyapunov/SOS): for a
rational unit ``n`` it verifies the Positivstellensatz identity

    ||u||^2 - R^2  =  ( (n_perp . u)^2 + (n . u - R)^2 )  +  2R (n . u - R)
                       \\-------- SOS s0(r) --------/       \\-- s1 h --/

with ``u = r - c``, ``h = n . u - R``, and multiplier ``s1 = 2R >= 0``.
On the cut ``{h >= 0}`` the right side is ``SOS + nonneg >= 0``, so
``||u|| >= R``. The SOS block is carried by an exact rational Gram that
can be re-synthesized from an untrusted float SDP and validated
(``sos.validate_gram``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction as F

from podium.verify import barrier, sos

Mono = tuple[int, ...]      # (deg rx, deg ry)
Poly = dict[Mono, F]


def snap_rational_unit(vx: float, vy: float, max_den: int = 1000
                       ) -> tuple[F, F]:
    """Nearest rational unit vector to the direction (vx, vy), exact:
    parametrize the circle by t = tan(theta/2) in Q, so
    (n_x, n_y) = ((1 - t^2), 2t) / (1 + t^2) with n_x^2 + n_y^2 = 1
    exactly."""
    theta = math.atan2(vy, vx)
    t = F(math.tan(theta / 2.0)).limit_denominator(max_den)
    d = 1 + t * t
    return (1 - t * t) / d, (2 * t) / d


@dataclass
class KOZCutReport:
    """Exact soundness certificate for one half-space keep-out cut."""

    sos_certified: bool
    identity_ok: bool
    gram: list[list[F]]
    lhs: dict[Mono, F]
    h: dict[Mono, F]
    problems: list[str] = field(default_factory=list)

    def sound(self) -> bool:
        return self.sos_certified and self.identity_ok and not self.problems


def _koz_polys(n: tuple[F, F], radius: F, center: tuple[F, F]
               ) -> tuple[Poly, Poly, Poly]:
    nx, ny = n
    kx, ky = center
    u2 = sos.padd(
        {(2, 0): F(1), (1, 0): -2 * kx, (0, 0): kx * kx},
        {(0, 2): F(1), (0, 1): -2 * ky, (0, 0): ky * ky})
    ndu: Poly = {(1, 0): nx, (0, 1): ny, (0, 0): -(nx * kx + ny * ky)}
    h = sos.psub(ndu, {(0, 0): radius})
    lhs = sos.psub(sos.psub(u2, {(0, 0): radius * radius}),
                   sos.pscale(2 * radius, h))
    return u2, h, lhs


def _exact_gram(n: tuple[F, F], radius: F, center: tuple[F, F]
                ) -> list[list[F]]:
    """Gram of lhs = (n_perp.u)^2 + (n.u - R)^2 over basis [1, rx, ry]."""
    nx, ny = n
    kx, ky = center
    # affine forms as coefficient vectors over [1, rx, ry]
    perp = [ny * kx - nx * ky, -ny, nx]                 # n_perp . u
    face = [-(nx * kx + ny * ky + radius), nx, ny]      # n.u - R
    g = [[F(0)] * 3 for _ in range(3)]
    for v in (perp, face):
        for i in range(3):
            for j in range(3):
                g[i][j] += v[i] * v[j]
    return g


def certify_halfspace_koz(n: tuple[F, F], radius: F,
                          center: tuple[F, F]) -> KOZCutReport:
    """Exact certificate that the hard cut n.(r - center) >= R is a
    sound inner-approximation of ||r - center|| >= R. Requires n a
    rational unit vector."""
    problems: list[str] = []
    if n[0] * n[0] + n[1] * n[1] != 1:
        problems.append("n is not an exact rational unit vector")
    if radius < 0:
        # the S-procedure multiplier is 2R, which must be >= 0; a negative
        # radius makes the certificate meaningless (the geometry inverts).
        problems.append("radius must be nonnegative")
    _u2, h, lhs = _koz_polys(n, radius, center)
    gram = _exact_gram(n, radius, center)
    basis: list[Mono] = [(0, 0), (1, 0), (0, 1)]
    sos_ok, sos_probs = sos.is_sos(lhs, basis, gram)
    problems += sos_probs
    # S-procedure identity: lhs + 2R h == ||u||^2 - R^2
    recon = sos.padd(lhs, sos.pscale(2 * radius, h))
    target = sos.psub(_u2, {(0, 0): radius * radius})
    identity_ok = (recon == target)
    if not identity_ok:
        problems.append("S-procedure identity failed")
    return KOZCutReport(sos_ok, identity_ok, gram, lhs, h, problems)


def gram_constraints(n: tuple[F, F], radius: F, center: tuple[F, F]
                     ) -> tuple[list[Mono], Poly]:
    """(basis, lhs) for re-synthesizing the SOS Gram from an untrusted
    float SDP and validating it with sos.validate_gram."""
    _u2, _h, lhs = _koz_polys(n, radius, center)
    basis: list[Mono] = [(0, 0), (1, 0), (0, 1)]
    return basis, lhs


# --- higher-degree (superquadric) certified cut -----------------------

# degree-2 monomial basis for a quartic SOS block
QUARTIC_BASIS: list[Mono] = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]


@dataclass
class CutReport:
    """Exact S-procedure certificate for a polynomial keep-out cut:
    q = z^T G z + sum_i c_i cut_i, with G PSD and each c_i >= 0, so
    q >= 0 wherever all cut_i >= 0."""

    certified: bool
    problems: list[str]
    gram: list[list[F]]


def certify_cut(q: Poly, cuts: list[tuple[F, Poly]],
                basis: list[Mono], gram: list[list[F]]) -> CutReport:
    """Verify q = z^T G z + sum c_i cut_i exactly, with G PSD and every
    multiplier c_i >= 0 -- a Positivstellensatz witness that q >= 0 on
    the intersection of the cut half-/sub-spaces. All exact rational."""
    problems: list[str] = []
    recon = dict(sos.gram_poly(basis, gram))
    for c, cut in cuts:
        if c < 0:
            problems.append("negative S-procedure multiplier")
        recon = sos.padd(recon, sos.pscale(c, cut))
    if sos.psub(q, recon):
        problems.append("identity q = z^T G z + sum c*cut failed")
    if not barrier.is_psd(gram):
        problems.append("SOS Gram is not positive semidefinite")
    return CutReport(not problems, problems, gram)


def superquadric_diagonal_certificate(
        p: F) -> tuple[Poly, Poly, F, list[Mono], list[list[F]]]:
    """Quartic keep-out rx^4 + ry^4 >= 2 p^4 (nonconvex exterior) with
    the tangent half-space cut rx + ry >= 2p at the diagonal boundary
    point (p, p). Returns (q, cut, multiplier, basis, gram) for the
    exact degree-4 Positivstellensatz

        rx^4 + ry^4 - 2 p^4 = sigma0 + 4 p^3 (rx + ry - 2p),

    sigma0 = (rx^2 - p^2)^2 + 2p^2 (rx - p)^2 + (same in ry), SOS.
    """
    q: Poly = {(4, 0): F(1), (0, 4): F(1), (0, 0): -2 * p**4}
    cut: Poly = {(1, 0): F(1), (0, 1): F(1), (0, 0): -2 * p}
    mult = 4 * p**3
    squares: list[tuple[list[F], F]] = [
        ([-p * p, F(0), F(0), F(1), F(0), F(0)], F(1)),   # rx^2 - p^2
        ([-p * p, F(0), F(0), F(0), F(0), F(1)], F(1)),   # ry^2 - p^2
        ([-p, F(1), F(0), F(0), F(0), F(0)], 2 * p * p),  # sqrt(2)p (rx-p)
        ([-p, F(0), F(1), F(0), F(0), F(0)], 2 * p * p),  # sqrt(2)p (ry-p)
    ]
    gram = [[F(0)] * 6 for _ in range(6)]
    for vec, w in squares:
        for i in range(6):
            for j in range(6):
                gram[i][j] += w * vec[i] * vec[j]
    return q, cut, mult, QUARTIC_BASIS, gram
