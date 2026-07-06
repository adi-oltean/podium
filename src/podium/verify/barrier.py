"""Infinite-horizon abort-safety barrier certificates for CW drift.

Claim certified: every passive (thrust-free) CW trajectory starting in
the initial ellipsoid X0 keeps radial/cross-track-plane separation
> keep-out radius FOR ALL TIME — and since ||(x,y,z)|| >= ||(x,z)||,
that bounds full 3-D separation regardless of along-track drift. This
is the e/i-vector-separation heritage argument, machine-checked.

Architecture (the untrusted-synthesizer / trusted-checker split):

* Coordinates: time-scaled LVLH u = (X, Y, Z, VX, VY, VZ) with
  X = x, VX = vx/n, tau = n t. The CW matrix is then INTEGER
  (du/dtau = A u), so the whole problem is rational.
* Barrier basis: the CW flow invariants
      c1   = 4X + 2VY              (mean radial offset; linear)
      Ax2  = (3X + 2VY)^2 + VX^2   (in-plane amplitude^2)
      Az2  = Z^2 + VZ^2            (cross-track amplitude^2)
  B = a1*c1^2 + a2*Ax2 + a3*Az2 + a4*c1 + a5. Every basis element is a
  flow invariant with an INTEGER coefficient matrix, so dB/dtau = 0
  holds structurally and the checker verifies A'P + PA == 0 EXACTLY —
  the barrier is a conserved quantity, hence its sublevel sets are
  invariant for all time (no discretization, no horizon).
* Safety conditions as S-procedure PSD constraints on homogenized
  [1; u] quadratic forms:
      C1 = -M(B) - lam0*G0 - E(eps0)  is PSD   (B <= -eps0 on X0)
      C2 =  M(B) - lamU*GU - E(epsU)  is PSD   (B >= +epsU on the KOZ)
  with G0 the X0-ellipsoid form, GU the RN keep-out form, lam >= 0.
* Synthesis: tiny SDP over (a1..a5, lam0, lamU) via cvxpy/Clarabel —
  UNTRUSTED. Rationalized with bounded denominators afterwards.
* Verification: `verify_certificate` is the TRUSTED path — pure
  `fractions.Fraction` arithmetic, matrices rebuilt from the integer
  invariant bases, PSD decided by the all-principal-minors criterion
  (exact determinants). No floats anywhere in the checker.

If the initial set is not passively safe (e.g. a hold centered on the
V-bar axis, whose coast passes straight through the target), synthesis
is infeasible — correctly: no certificate exists for an unsafe fact.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations

Frac = Fraction
Mat = list[list[Fraction]]

_N = 6  # state dimension (scaled LVLH)
_H = 7  # homogenized [1; u]

# integer CW matrix in scaled coordinates: du/dtau = A_CW u
A_CW = [
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1],
    [3, 0, 0, 0, 2, 0],
    [0, 0, 0, -2, 0, 0],
    [0, 0, -1, 0, 0, 0],
]


def _zeros(r: int, c: int) -> Mat:
    return [[Frac(0)] * c for _ in range(r)]


def _sym_outer(v: list[Fraction]) -> Mat:
    return [[vi * vj for vj in v] for vi in v]


def _mat_add(*ms: Mat) -> Mat:
    out = _zeros(len(ms[0]), len(ms[0][0]))
    for m in ms:
        for i in range(len(out)):
            for j in range(len(out[0])):
                out[i][j] += m[i][j]
    return out


def _mat_scale(s: Fraction, m: Mat) -> Mat:
    return [[s * x for x in row] for row in m]


# --- invariant quadratic bases (6x6, integer) --------------------------
# c1 = 4X + 2VY ; u2 = 3X + 2VY ; Ax2 = u2^2 + VX^2 ; Az2 = Z^2 + VZ^2
_C1_VEC = [Frac(v) for v in (4, 0, 0, 0, 2, 0)]
_U2_VEC = [Frac(v) for v in (3, 0, 0, 0, 2, 0)]
_VX_VEC = [Frac(v) for v in (0, 0, 0, 1, 0, 0)]
_Z_VEC = [Frac(v) for v in (0, 0, 1, 0, 0, 0)]
_VZ_VEC = [Frac(v) for v in (0, 0, 0, 0, 0, 1)]

P_C1SQ = _sym_outer(_C1_VEC)
P_AX2 = _mat_add(_sym_outer(_U2_VEC), _sym_outer(_VX_VEC))
P_AZ2 = _mat_add(_sym_outer(_Z_VEC), _sym_outer(_VZ_VEC))


@dataclass(frozen=True)
class AbortSafetyCase:
    """Problem data (all rational): X0 ellipsoid + RN keep-out radius.

    X0 = { u : sum ((u_i - center_i)/radius_i)^2 <= 1 }, scaled coords
    (positions in meters, velocities in meters — i.e. v/n).
    """

    center: tuple[Fraction, ...]
    radii: tuple[Fraction, ...]
    koz_radius: Fraction


@dataclass(frozen=True)
class BarrierCertificate:
    """B = a1*c1^2 + a2*Ax2 + a3*Az2 + a4*c1 + a5; S-procedure data."""

    a: tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
    lam0: Fraction
    lam_u: Fraction
    eps0: Fraction
    eps_u: Fraction
    case: AbortSafetyCase


# --- homogenized forms --------------------------------------------------
def _homog(p: Mat, lin: list[Fraction], const: Fraction) -> Mat:
    """[1;u]' M [1;u] = u'Pu + lin.u + const."""
    m = _zeros(_H, _H)
    m[0][0] = const
    for i in range(_N):
        m[0][i + 1] = lin[i] / Frac(2)
        m[i + 1][0] = lin[i] / Frac(2)
        for j in range(_N):
            m[i + 1][j + 1] = p[i][j]
    return m


def barrier_matrices(cert: BarrierCertificate) -> tuple[Mat, Mat]:
    """(P quadratic part 6x6, M homogenized 7x7) rebuilt from the a's."""
    a1, a2, a3, a4, a5 = cert.a
    p = _mat_add(_mat_scale(a1, P_C1SQ), _mat_scale(a2, P_AX2),
                 _mat_scale(a3, P_AZ2))
    lin = [a4 * v for v in _C1_VEC]
    return p, _homog(p, lin, a5)


def _g0_matrix(case: AbortSafetyCase) -> Mat:
    """g0(u) = 1 - sum ((u_i - c_i)/r_i)^2 as a homogenized form."""
    p = _zeros(_N, _N)
    lin = [Frac(0)] * _N
    const = Frac(1)
    for i in range(_N):
        w = Frac(1) / (case.radii[i] * case.radii[i])
        p[i][i] = -w
        lin[i] = 2 * w * case.center[i]
        const -= w * case.center[i] * case.center[i]
    return _homog(p, lin, const)


def _gu_matrix(case: AbortSafetyCase) -> Mat:
    """gu(u) = R^2 - X^2 - Z^2 (RN keep-out) as a homogenized form."""
    p = _zeros(_N, _N)
    p[0][0] = Frac(-1)
    p[2][2] = Frac(-1)
    return _homog(p, [Frac(0)] * _N, case.koz_radius * case.koz_radius)


def _det(m: Mat) -> Fraction:
    """Exact determinant by fraction-free-ish Gaussian elimination."""
    n = len(m)
    a = [row[:] for row in m]
    det = Frac(1)
    for col in range(n):
        piv = None
        for r in range(col, n):
            if a[r][col] != 0:
                piv = r
                break
        if piv is None:
            return Frac(0)
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
            det = -det
        det *= a[col][col]
        inv = Frac(1) / a[col][col]
        for r in range(col + 1, n):
            f = a[r][col] * inv
            if f != 0:
                for c in range(col, n):
                    a[r][c] -= f * a[col][c]
    return det


def is_psd(m: Mat) -> bool:
    """Exact PSD test: every principal minor is >= 0 (Sylvester for
    semidefiniteness). 2^7 - 1 = 127 exact determinants at our size."""
    n = len(m)
    for i in range(n):
        for j in range(n):
            if m[i][j] != m[j][i]:
                return False
    idx = list(range(n))
    for k in range(1, n + 1):
        for sub in combinations(idx, k):
            mm = [[m[i][j] for j in sub] for i in sub]
            if _det(mm) < 0:
                return False
    return True


def verify_certificate(cert: BarrierCertificate) -> list[str]:
    """TRUSTED checker: exact rational arithmetic only.

    Returns a list of violations (empty == certificate PROVEN). Every
    quantity is rebuilt here from the certificate scalars and the
    integer invariant bases — nothing synthesized is trusted.
    """
    problems: list[str] = []
    for name, v in (("lam0", cert.lam0), ("lam_u", cert.lam_u)):
        if not isinstance(v, Fraction) or v < 0:
            problems.append(f"{name} must be a nonnegative Fraction")
    for name, v in (("eps0", cert.eps0), ("eps_u", cert.eps_u)):
        if not isinstance(v, Fraction) or v <= 0:
            problems.append(f"{name} must be a positive Fraction")
    if not all(isinstance(x, Fraction) for x in cert.a):
        problems.append("barrier coefficients must be Fractions")
    if problems:
        return problems

    p, m = barrier_matrices(cert)

    # C3: dB/dtau = u'(A'P + PA)u must vanish EXACTLY (B is conserved,
    # so its sublevel sets are invariant for all time).
    a_cw = [[Frac(v) for v in row] for row in A_CW]
    lie = _zeros(_N, _N)
    for i in range(_N):
        for j in range(_N):
            s = Frac(0)
            for k in range(_N):
                s += a_cw[k][i] * p[k][j] + p[i][k] * a_cw[k][j]
            lie[i][j] = s
    if any(lie[i][j] != 0 for i in range(_N) for j in range(_N)):
        problems.append("A'P + PA != 0: barrier is not a flow invariant")

    e0 = _zeros(_H, _H)
    e0[0][0] = cert.eps0
    eu = _zeros(_H, _H)
    eu[0][0] = cert.eps_u

    # C1: B <= -eps0 on X0
    c1m = _mat_add(_mat_scale(Frac(-1), m),
                   _mat_scale(-cert.lam0, _g0_matrix(cert.case)),
                   _mat_scale(Frac(-1), e0))
    if not is_psd(c1m):
        problems.append("C1 not PSD: B <= -eps0 on X0 not established")

    # C2: B >= +eps_u on the keep-out set
    c2m = _mat_add(m, _mat_scale(-cert.lam_u, _gu_matrix(cert.case)),
                   _mat_scale(Frac(-1), eu))
    if not is_psd(c2m):
        problems.append("C2 not PSD: B >= eps_u on KOZ not established")

    return problems


# --- untrusted synthesis (sandbox: cvxpy SDP) ---------------------------
def synthesize(case: AbortSafetyCase, margin: float = 1.0) -> "dict | None":
    """SDP synthesis of the barrier scalars (UNTRUSTED — always run
    verify_certificate on the rationalized result). Returns the float
    solution dict or None if infeasible."""
    import cvxpy as cp
    import numpy as np

    def f(x: Fraction) -> float:
        return float(x)

    a = cp.Variable(5)
    lam0 = cp.Variable(nonneg=True)
    lam_u = cp.Variable(nonneg=True)
    p_bases = [P_C1SQ, P_AX2, P_AZ2]
    p_np = [np.array([[f(x) for x in row] for row in b]) for b in p_bases]
    c1v = np.array([f(x) for x in _C1_VEC])

    def homog_expr(pq, lin, const):
        rows = []
        rows.append(cp.hstack([cp.reshape(const, (1, 1), order="C"),
                               cp.reshape(lin / 2, (1, _N), order="C")]))
        rows.append(cp.hstack([cp.reshape(lin / 2, (_N, 1), order="C"), pq]))
        return cp.vstack(rows)

    pq = a[0] * p_np[0] + a[1] * p_np[1] + a[2] * p_np[2]
    m = homog_expr(pq, a[3] * c1v, a[4])

    g0 = np.array([[f(x) for x in row] for row in _g0_matrix(case)])
    gu = np.array([[f(x) for x in row] for row in _gu_matrix(case)])
    e = np.zeros((_H, _H))
    e[0, 0] = 1.0

    t = cp.Variable(nonneg=True)
    cons = [
        -m - lam0 * g0 - margin * e >> t * np.eye(_H),
        m - lam_u * gu - margin * e >> t * np.eye(_H),
        cp.abs(a) <= 1e4,
    ]
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver=cp.CLARABEL)
    if prob.status not in ("optimal", "optimal_inaccurate") or t.value is None:
        return None
    return {"a": [float(v) for v in a.value], "lam0": float(lam0.value),
            "lam_u": float(lam_u.value), "slack": float(t.value)}


def rationalize(sol: dict, case: AbortSafetyCase, eps: float = 0.5,
                max_den: int = 10_000) -> BarrierCertificate:
    """Round the float SDP solution to bounded-denominator rationals;
    the synthesis slack absorbs the rounding. eps must be < margin."""
    fr = [Frac(x).limit_denominator(max_den) for x in sol["a"]]
    return BarrierCertificate(
        a=(fr[0], fr[1], fr[2], fr[3], fr[4]),
        lam0=max(Frac(0), Frac(sol["lam0"]).limit_denominator(max_den)),
        lam_u=max(Frac(0), Frac(sol["lam_u"]).limit_denominator(max_den)),
        eps0=Frac(eps).limit_denominator(1000),
        eps_u=Frac(eps).limit_denominator(1000),
        case=case,
    )
