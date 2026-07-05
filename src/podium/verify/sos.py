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


def validate_gram(target: Poly, basis: list[Mono],
                  gram_float: list[list[float]],
                  margin: Fraction = Fraction(1, 10**6),
                  max_den: int = 10**9) -> list[list[Fraction]] | None:
    """Round-and-correct an UNTRUSTED float SOS Gram into an EXACT
    rational Gram that reproduces `target` identically and stays PSD
    (validated SOS). The synthesis (an SDP) may be floating point; the
    shipped certificate is exact.

    Method: rationalize the float Gram, inflate the diagonal by a small
    rational `margin` (the float interior-point Gram is strictly PD, so
    this preserves PSD with slack), then absorb the exact coefficient
    residual monomial by monomial. Each Gram entry contributes to
    exactly one product monomial, so the correction decouples: for each
    residual monomial pick one entry with that product and adjust it.
    Returns the exact Gram, or None if `target` has a monomial no basis
    pair can produce (the basis is too small for an SOS form).
    """
    n = len(basis)
    g = [[Frac(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            g[i][j] = Frac(float(gram_float[i][j])).limit_denominator(max_den)
    # symmetrize exactly, then add the PSD-slack margin on the diagonal
    for i in range(n):
        for j in range(i):
            s = (g[i][j] + g[j][i]) / 2
            g[i][j] = g[j][i] = s
        g[i][i] += margin

    # map each product monomial -> a preferred entry (diagonal first)
    entry_for: dict[Mono, tuple[int, int]] = {}
    for i in range(n):
        for j in range(i, n):
            m = tuple(basis[i][k] + basis[j][k] for k in range(len(basis[i])))
            if m not in entry_for or i == j:  # prefer a diagonal entry
                entry_for.setdefault(m, (i, j))
                if i == j:
                    entry_for[m] = (i, j)

    residual = psub(target, gram_poly(basis, g))
    for m, r in residual.items():
        if m not in entry_for:
            return None                       # basis cannot span target
        i, j = entry_for[m]
        delta = r if i == j else r / 2        # weight 1 (diag) or 2 (off)
        g[i][j] += delta
        if i != j:
            g[j][i] += delta
    return g
