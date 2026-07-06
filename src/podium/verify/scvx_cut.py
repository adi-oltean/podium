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

from podium.verify import sos

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
