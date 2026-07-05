"""Exact-rational sum-of-squares / Positivstellensatz certificate
checking.

Generalizes the quadratic S-procedure barrier (podium.verify.barrier,
degree two) to genuinely higher-degree polynomial systems. A polynomial
p is certified SOS by exhibiting a monomial basis z and a symmetric
Gram matrix G with

    p(x) = z(x)^T G z(x)   and   G >= 0,

both checked EXACTLY: the polynomial identity is matched coefficient by
coefficient over the rationals, and G's positive semidefiniteness is
decided by the exact all-principal-minors test reused from
podium.verify.barrier. There is no floating point in the trusted path;
a floating-point SOS solver may synthesize (G, z), but the certificate
that ships is the rational Gram matrix, re-verified here.

Applied to a barrier/Lyapunov function V and a polynomial vector field
f, certifying that the Lie derivative -dV/dt is SOS proves dV/dt <= 0
everywhere, hence the sub-level set {V <= c} is an INFINITE-HORIZON
invariant of a nonlinear system --- the higher-degree analogue of the
quadratic abort-safety barrier.

Polynomials are represented as dicts mapping an exponent tuple (one
entry per variable) to a Fraction coefficient; zero coefficients are
dropped so equality is exact dict equality.
"""

from __future__ import annotations

from fractions import Fraction

from podium.verify.barrier import Frac, is_psd

Mono = tuple[int, ...]
Poly = dict[Mono, Fraction]


def _clean(p: Poly) -> Poly:
    return {m: c for m, c in p.items() if c != 0}


def padd(*ps: Poly) -> Poly:
    out: Poly = {}
    for p in ps:
        for m, c in p.items():
            out[m] = out.get(m, Frac(0)) + c
    return _clean(out)


def pscale(s: Fraction, p: Poly) -> Poly:
    return _clean({m: s * c for m, c in p.items()})


def psub(a: Poly, b: Poly) -> Poly:
    return padd(a, pscale(Frac(-1), b))


def pmul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            m = tuple(ma[i] + mb[i] for i in range(len(ma)))
            out[m] = out.get(m, Frac(0)) + ca * cb
    return _clean(out)


def pdiff(p: Poly, var: int) -> Poly:
    """Partial derivative with respect to variable `var`."""
    out: Poly = {}
    for m, c in p.items():
        if m[var] == 0:
            continue
        nm = list(m)
        k = nm[var]
        nm[var] -= 1
        out[tuple(nm)] = out.get(tuple(nm), Frac(0)) + c * k
    return _clean(out)


def lie_derivative(v: Poly, f: list[Poly]) -> Poly:
    """dV/dt = sum_i (dV/dx_i) * f_i along the vector field f (a list of
    polynomials, one per variable)."""
    terms = [pmul(pdiff(v, i), f[i]) for i in range(len(f))]
    return padd(*terms)


def _mono(*exps: int) -> Mono:
    return tuple(exps)


def gram_poly(basis: list[Mono], gram: list[list[Fraction]]) -> Poly:
    """Expand z^T G z into a polynomial (exact), z = basis."""
    out: Poly = {}
    n = len(basis)
    for i in range(n):
        for j in range(n):
            g = gram[i][j]
            if g == 0:
                continue
            m = tuple(basis[i][k] + basis[j][k]
                      for k in range(len(basis[i])))
            out[m] = out.get(m, Frac(0)) + g
    return _clean(out)


def is_sos(p: Poly, basis: list[Mono],
           gram: list[list[Fraction]]) -> tuple[bool, list[str]]:
    """Certify p is SOS via (basis, Gram): exact identity p = z^T G z
    AND G >= 0 (exact PSD). Returns (certified, problems)."""
    problems: list[str] = []
    n = len(basis)
    if any(len(row) != n for row in gram) or len(gram) != n:
        problems.append("Gram must be square, matching the basis")
        return False, problems
    for i in range(n):
        for j in range(i):
            if gram[i][j] != gram[j][i]:
                problems.append(f"Gram not symmetric at ({i},{j})")
                break
    if psub(p, gram_poly(basis, gram)):        # nonempty => not identical
        problems.append("polynomial identity p = z^T G z fails")
    if not is_psd(gram):
        problems.append("Gram is not positive semidefinite")
    return (not problems), problems
