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

  positive semidefinite. M >= 0 makes f0(x) - lambda f1(x) - t a GLOBAL
  Lagrangian minorant (>= 0 for all x), so for feasible x (f1(x) >= 0) and
  lambda >= 0, f0(x) >= t; hence t <= J* (weak duality / S-procedure).
  `certify_lower_bound` checks lambda >= 0 and is_psd(M) exactly -- no
  matrix inverse.

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

from podium.verify.barrier import _det, is_psd

Vec = list[F]
Mat = list[list[F]]


def _check_qcqp(p0: Mat, q0: Vec, p1: Mat, q1: Vec) -> int:
    """Validate QCQP data shapes (square n x n objectives, length-n
    vectors); raise ValueError on malformed input rather than raising an
    opaque IndexError or silently truncating later. Inputs must be exact
    Fractions -- floats/NaN break the exact-arithmetic soundness contract."""
    n = len(q0)
    if (len(q1) != n or len(p0) != n or len(p1) != n
            or any(len(row) != n for row in p0)
            or any(len(row) != n for row in p1)):
        raise ValueError("QCQP data dimensions are inconsistent")
    return n


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
        b = (q0[i] - lam * q1[i]) / F(2)      # / F(2): exact even for int data
        m[i][n] = b
        m[n][i] = b
    m[n][n] = r0 - lam * r1 - t
    return m


def certify_lower_bound(p0: Mat, q0: Vec, r0: F, p1: Mat, q1: Vec, r1: F,
                        lam: F, t: F) -> bool:
    """Exact certificate that t <= J* for the nonconvex QCQP: returns True
    iff lam >= 0 and M(lam, t) >= 0. All inputs Fractions (rationalize an
    SDP-dual solution first)."""
    _check_qcqp(p0, q0, p1, q1)
    return lam >= 0 and is_psd(lower_bound_matrix(
        p0, q0, r0, p1, q1, r1, lam, t))


def _quad(p: Mat, q: Vec, r: F, x: Vec) -> F:
    """x'P x + q'x + r, exact."""
    n = len(x)
    val = r
    for i in range(n):
        val += q[i] * x[i]
        for j in range(n):
            val += x[i] * p[i][j] * x[j]
    return val


def certify_upper_bound(p0: Mat, q0: Vec, r0: F, p1: Mat, q1: Vec, r1: F,
                        x: Vec) -> F | None:
    """Certified upper bound. If x is EXACTLY feasible (f1(x) >= 0 in exact
    rationals), return f0(x): then J* <= f0(x). Otherwise return None.

    Feasibility is checked exactly and independently -- it does NOT trust a
    solver's tolerance. A point that is only feasible-within-tolerance
    (e.g. a KKT solution with a tiny inequality violation) is rejected,
    because its objective can dip BELOW J* and would otherwise collapse the
    bracket beneath the true optimum. Any exactly-feasible x gives a valid
    upper bound; a solver/KKT step is only for choosing a good x."""
    if _check_qcqp(p0, q0, p1, q1) != len(x):
        raise ValueError("x has the wrong dimension")
    if _quad(p1, q1, r1, x) < 0:            # exactly infeasible
        return None
    return _quad(p0, q0, r0, x)


def closes(lower_t: F, upper_j: "F | None") -> bool:
    """Low-level combiner: True iff a valid lower bound and a valid
    (exactly-feasible) upper bound coincide. NOTE: this only compares two
    numbers -- it does not check they came from the same problem. Prefer
    `certified_optimum`, which binds both legs to one QCQP's data."""
    return upper_j is not None and lower_t == upper_j


def certified_optimum(p0: Mat, q0: Vec, r0: F, p1: Mat, q1: Vec, r1: F,
                      lam: F, t: F, x: Vec) -> tuple["F | None", "F | None", bool]:
    """Compose both bracket legs from ONE problem's data, so a lower-bound
    certificate cannot be paired with an upper point from a *different*
    problem (provenance binding). Returns (t_lb, j_ub, closed):
    t_lb = t if certify_lower_bound holds else None; j_ub = the exact upper
    bound if x is exactly feasible else None; closed is True iff both hold
    and t_lb == j_ub -- in which case J* = t_lb = j_ub exactly. Because
    both legs are checked against the same (P0..r1), a caller cannot form a
    false closure by mixing certificates across problems."""
    t_lb = t if certify_lower_bound(p0, q0, r0, p1, q1, r1, lam, t) else None
    j_ub = certify_upper_bound(p0, q0, r0, p1, q1, r1, x)
    return t_lb, j_ub, (t_lb is not None and j_ub is not None and t_lb == j_ub)


# --- multiple constraints (Theorem 4) -----------------------------------

Con = tuple[Mat, Vec, F]        # (P_k, q_k, r_k) for f_k(x) >= 0


def _check_multi(p0: Mat, q0: Vec, cons: list[Con]) -> int:
    """Validate an m-constraint QCQP's shapes; raise ValueError on
    malformed data instead of silently ignoring extra rows or raising an
    opaque IndexError."""
    n = len(q0)
    if len(p0) != n or any(len(row) != n for row in p0):
        raise ValueError("objective dimensions are inconsistent")
    for pk, qk, _rk in cons:
        if len(pk) != n or any(len(row) != n for row in pk) or len(qk) != n:
            raise ValueError("constraint dimensions are inconsistent")
    return n


def lower_bound_matrix_multi(p0: Mat, q0: Vec, r0: F, cons: list[Con],
                             lams: Vec, t: F) -> Mat:
    """S-procedure LMI M(lam, t) for J* = min f0 s.t. f_k >= 0 (k=1..m):
    A = P0 - sum_k lam_k P_k, off = (q0 - sum_k lam_k q_k)/2, corner =
    r0 - sum_k lam_k r_k - t."""
    n = len(q0)
    a = [[p0[i][j] for j in range(n)] for i in range(n)]
    b = [q0[i] for i in range(n)]
    d = r0 - t
    for k, (pk, qk, rk) in enumerate(cons):
        lk = lams[k]
        for i in range(n):
            b[i] -= lk * qk[i]
            for j in range(n):
                a[i][j] -= lk * pk[i][j]
        d -= lk * rk
    m: Mat = [[F(0)] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(n):
            m[i][j] = a[i][j]
        m[i][n] = b[i] / F(2)
        m[n][i] = b[i] / F(2)
    m[n][n] = d
    return m


def certify_lower_bound_multi(p0: Mat, q0: Vec, r0: F, cons: list[Con],
                              lams: Vec, t: F) -> bool:
    """Exact certificate t <= J* for the m-constraint nonconvex QCQP:
    True iff every lam_k >= 0 and M(lam, t) >= 0. Same weak-duality proof
    as the single-constraint case; here the S-lemma need not hold, so the
    best certifiable t (the Shor bound) can be strictly below J* -- an
    exactly-certified duality gap."""
    _check_multi(p0, q0, cons)
    if len(lams) != len(cons):
        raise ValueError("len(lams) must equal the number of constraints")
    if any(lk < 0 for lk in lams):
        return False
    return is_psd(lower_bound_matrix_multi(p0, q0, r0, cons, lams, t))


def certify_upper_bound_multi(p0: Mat, q0: Vec, r0: F, cons: list[Con],
                              x: Vec) -> F | None:
    """Certified upper bound J* <= f0(x): returns f0(x) if x is EXACTLY
    feasible for ALL constraints (f_k(x) >= 0), else None."""
    if _check_multi(p0, q0, cons) != len(x):
        raise ValueError("x has the wrong dimension")
    for pk, qk, rk in cons:
        if _quad(pk, qk, rk, x) < 0:
            return None
    return _quad(p0, q0, r0, x)


# --- exact rational recovery of the dual lower bound (Theorem 2) --------

def _solve(a: Mat, b: Vec) -> Vec | None:
    """Exact solution of A w = b (Fractions), or None if A is singular."""
    n = len(b)
    m = [list(a[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if m[r][col] != 0), None)
        if piv is None:
            return None
        m[col], m[piv] = m[piv], m[col]
        inv = F(1) / m[col][col]
        m[col] = [v * inv for v in m[col]]
        for r in range(n):
            if r != col and m[r][col] != 0:
                f = m[r][col]
                m[r] = [m[r][k] - f * m[col][k] for k in range(n + 1)]
    return [m[i][n] for i in range(n)]


def dual_value(p0: Mat, q0: Vec, r0: F, p1: Mat, q1: Vec, r1: F,
               lam: F) -> F | None:
    """Exact Lagrangian dual value g(lam) = min_x [f0(x) - lam f1(x)],
    defined (and equal to the largest t with M(lam, t) >= 0) when
    A = P0 - lam P1 is positive definite. Returns None otherwise (the
    singular 'hard case', A not PD). g(lam) <= J* for lam >= 0."""
    n = len(q0)
    a = [[p0[i][j] - lam * p1[i][j] for j in range(n)] for i in range(n)]
    if not (is_psd(a) and _det(a) != 0):        # need A > 0
        return None
    b = [q0[i] - lam * q1[i] for i in range(n)]
    w = _solve(a, b)
    if w is None:
        return None
    bw = sum((b[i] * w[i] for i in range(n)), F(0))     # b' A^{-1} b
    return (r0 - lam * r1) - bw / F(4)


def recover_lower_bound(p0: Mat, q0: Vec, r0: F, p1: Mat, q1: Vec, r1: F,
                        lam_float: float, max_den: int = 10**6
                        ) -> tuple[F, F] | None:
    """Round a floating-point dual multiplier (e.g. from an untrusted SDP
    solve) to a rational lam and return (lam, t) with t = g(lam) an EXACT
    certified lower bound (t <= J*, verified by certify_lower_bound), or
    None if the rounded lam is not dual-feasible. When lam_float is near a
    rational optimal multiplier lam*, low-denominator rounding recovers
    lam* exactly and t = J* (Theorem 2)."""
    lam = F(lam_float).limit_denominator(max_den)
    if lam < 0:
        return None
    t = dual_value(p0, q0, r0, p1, q1, r1, lam)
    if t is None or not certify_lower_bound(p0, q0, r0, p1, q1, r1, lam, t):
        return None
    return lam, t
