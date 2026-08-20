"""Fixed-width WDD-tier certificate verification -- the decision layer.

The trajectory-QCQP certificate protocol has two exact tiers with very different
execution profiles, and only one of them can fly. The GENERAL tier
(`podium.verify.riccati`, `podium.verify.bracket`) decides ``M(lam, t) >= 0`` by a
band ``LDL^T`` sweep: sound for any coupling, but its pivots are unbounded rationals
whose bit-width grows along the horizon, so it needs arbitrary-precision arithmetic
and admits no worst-case execution time. The **WDD tier** never forms a Schur
complement at all -- it decides acceptance by checks that are local to a three-stage
window, so every operand stays at input height forever, independent of the horizon
``N``. This module is that tier written the way flight code must be written:

    checked fixed-width integer arithmetic; overflow is a LOUD REFUSAL, never a
    wraparound.

Scope
-----
**Identity-coupled regime only** (``A_k = I``, ``R_k = I``, ``Pi_j = I``,
``W_k = w_k I``): the objective is

    f_0(x) = sum_{k<N-1} ||x_{k+1} - x_k||^2 + sum_k w_k ||x_k - g_k||^2

with one keep-out per stage, ``f_k(x) = ||x_k - c_k||^2 - rho_k^2 >= 0``, and
multipliers ``lam_k >= 0``. The Lagrangian Hessian ``H(lam) = P0 - sum_k lam_k
Pi_k^T Pi_k`` then has diagonal blocks ``(m_k + w_k - lam_k) I`` (``m_k`` = number of
neighbours: 1 at a boundary stage, 2 inside) and off-diagonal blocks ``-I``, so every
scalar row check ``H_ii >= sum_{j != i} |H_ij|`` collapses to the single per-stage
comparison ``w_k >= lam_k`` -- one cross-multiplication of two input rationals, no
division and no Schur complement anywhere. **General coupling is out of scope here**
and routes to the existing sweep (`podium.verify.riccati` / `podium.verify.bracket`)
with unbounded rationals; a general instance is typically *outside* the WDD regime,
and bringing one inside it (diagonal rescaling, goal-weight redesign) is a synthesis
question -- classical H-matrix / generalized-diagonal-dominance machinery that
deliberately has no code here. This module only *checks*: no eigenvalue computation,
no synthesis, nothing untrusted.

What is decided
---------------
Given the instance data and an emitted certificate ``(lam, x_bar)``:

  (a) ``lam_k >= 0`` for every stage;
  (b) weak diagonal dominance of ``H(lam)``, i.e. ``w_k >= lam_k`` per stage;
  (c) stationarity: ``grad(f_0 - sum_k lam_k f_k)(x_bar) = 0`` exactly, blockwise;
  (d) feasibility: ``f_k(x_bar) >= 0`` for every keep-out;
  (e) per-stage complementary slackness: ``lam_k = 0`` or ``f_k(x_bar) = 0``.

Soundness (proved in the session findings behind this module, and standing in for the
bordered sweep): (b) + Gershgorin gives ``2H >= 0``, so ``psi := f_0 - sum lam_k f_k``
is convex; by (c) ``x_bar`` is a stationary point, hence a global minimizer, so
``psi >= t := psi(x_bar)`` everywhere, which homogenizes to ``M(lam, t) >= 0`` -- the
protocol's acceptance predicate, established WITHOUT running any sweep. With (d),
``t <= J* <= f_0(x_bar)``. (e) makes every term of ``sum_k lam_k f_k(x_bar)``
vanish, so ``t = f_0(x_bar)``: the bracket **closes** and ``x_bar`` is a certified
global optimum. Failing (e) alone is therefore not a refusal -- the verdict is
ACCEPT with ``closed = False`` (a valid bound, an open bracket). Refusal asserts
nothing: the tier is sound and deliberately incomplete.

Deliberately **out of scope: the value layer.** Materializing ``t* = f_0(x_bar)``
sums ``N`` stage terms, and bounding that width needs a bounded-denominator
hypothesis on the stage data (satisfied by CW-like families, violated by
Yamanaka-Ankersen). Every check above is decided without ever forming ``t*``; a
caller that wants the number computes it itself, under its own hypothesis.

Arithmetic model
----------------
Rationals are **unnormalized ``(num, den)`` integer pairs with ``den > 0``**. No gcd,
no reduction, no `fractions.Fraction`: normalization is not fixed-width-friendly and
none of the checks need it (they are sign and zero tests, invariant under positive
rescaling). Python integers do not overflow, so this module *simulates* the
fixed-width contract instead of relying on the machine: `checked_add`, `checked_sub`
and `checked_mul` take the budget explicitly, re-check both operands and the result
against ``|x| < 2**(W-1)``, and report ``_OVERFLOW`` with value 0 rather than a
wrapped word. Faults propagate as status codes up to `verify_fixed_width`, which
turns them into a structured `Verdict` naming the operand and the stage -- never an
exception raised mid-flight, so the C shape stays ``status = step(...)``.

The width ledger (this module IS the width theorem's executable form)
---------------------------------------------------------------------
Let ``L`` be the instance height: the largest bit-length of any numerator or
denominator in the input data (`instance_height_bits`), and ``d`` the stage
dimension. Pair arithmetic costs, with ``b(.)`` a bit-length:

    add/sub  (n1,q1) +- (n2,q2) -> (n1*q2 +- n2*q1, q1*q2)
             b(num) <= max(b(n1)+b(q2), b(n2)+b(q1)) + 1,  b(den) <= b(q1)+b(q2)
    mul      (n1,q1) *  (n2,q2) -> (n1*n2, q1*q2)
             b(num) <= b(n1)+b(n2),                        b(den) <= b(q1)+b(q2)
    cmp      sign(a-b) = sign(n1*q2 - n2*q1):  costs b(n1)+b(q2) and b(n2)+b(q1),
             the difference one bit more -- no denominator product is formed.

Per check, with the ledger rows of the width theorem named alongside:

    | check                     | row | widest integer operand | horizon dep. |
    |---------------------------|-----|------------------------|--------------|
    | (a) lam_k >= 0            | O1  | L (a sign test only)   | none         |
    | (b) w_k >= lam_k          | O3' | 2L + 1                 | none         |
    | (c) stationarity residual | O4  | 9L + 7                 | none         |
    | (d) feasibility f_k       | O5  | (4d+1)L + d + 2        | none         |
    | (e) slackness             | O6  | max(O1, O5), zero tests| none         |

Row (c) accumulates five terms -- ``(m_k + w_k - lam_k) x_bar_k``, the two neighbour
witnesses, ``w_k g_k`` and ``lam_k c_k`` -- so the diagonal coefficient reaches
``(2L+3, 2L)`` and the running residual ``(9L+7, 9L)``. Row (d) squares ``d``
coordinate offsets and subtracts the squared radius: ``(4Ld + d + 1, 4Ld)`` after the
sum, ``((4d+1)L + d + 2, (4d+1)L)`` after the radius.

**No operand in the table depends on N.** That is the whole theorem: every check
reads data from at most three adjacent stages, so no denominator or magnitude can
accumulate along the horizon, and a single budget ``W`` sized from ``(L, d)`` alone
covers a trajectory of any length. `sufficient_width_bits` returns that budget.

The constants above are the *pair-model* constants. The sharp ledger of the width
theorem (``delta*L + ceil(log2 phi)``, e.g. ``2L+1`` for the dominance leg, ``4L+c(d)``
for stationarity) is stated for integers over one precomputed stage denominator; in
the pair model each addition multiplies denominators instead of sharing one, which is
what inflates ``c0`` from 4 to 9 on row (c) and adds the ``4d`` on row (d). The sharp
constants are recovered here exactly when the caller adopts the shared-denominator
rationalization convention -- present each stage's data as integers over a common
scale, i.e. ``den = 1`` -- because every check is a sign or zero test and so is
invariant under the positive rescaling that makes it so. Then row (c) falls to
``2L + 6`` and row (d) to ``2L + ceil(log2(d+1)) + 3``.

Static-subset fit (`docs/verification.md`)
------------------------------------------
Pure functions, no globals, no I/O, no RNG, no floats on the certified path. Loop
bounds come from input extents only (``N``, ``d``) -- in the C rendering they are
compile-time constants; there is no data-dependent ``while``, no recursion, and no
allocation in the step path (the five-term stationarity table is a static table).
Every division is absent by construction: the tier never divides. The one scalar
parameter, ``width_bits``, carries a range contract; the contract check is a sandbox
guard outside the certified arithmetic, which stays integer-only.
"""

from __future__ import annotations

from dataclasses import dataclass

from podium.verify.contracts import Interval, contract

# An unnormalized rational: (numerator, denominator) with denominator > 0.
Rat = tuple[int, int]

# Arithmetic status codes (the C shape is `status = op(...)`).
_OK = 0
_OVERFLOW = 1

_ZERO: Rat = (0, 1)
_ONE: Rat = (1, 1)

# Refusal reasons. Rendered as `reason@operand[stage k]` by `Verdict.__str__`.
REASON_OK = ""
REASON_BAD_INPUT = "bad-input"
REASON_MULTIPLIER = "multiplier-negative"
REASON_MARGIN = "margin-violated"
REASON_WIDTH = "width-exceeded"
REASON_STATIONARITY = "stationarity-failed"
REASON_FEASIBILITY = "feasibility-failed"


# --- the width budget --------------------------------------------------------
def width_limit(width_bits: int) -> int:
    """``2**(W-1)``: the exclusive magnitude bound of a ``W``-bit word."""
    return 1 << (width_bits - 1)


def fits(value: int, limit: int) -> bool:
    """True iff ``value`` is representable in the budget: ``|value| < limit``.

    The extreme negative two's-complement word ``-2**(W-1)`` is excluded on purpose:
    a checker that can refuse but must never trap cannot afford an operand whose
    negation overflows.
    """
    return -limit < value < limit


def sufficient_width_bits(height_bits: int, dim: int) -> int:
    """The budget the ledger proves sufficient for ANY horizon.

    ``max(9L + 7, (4d+1)L + d + 2) + 1`` bits, with ``L = height_bits`` the instance
    height and ``d = dim`` the stage dimension -- rows (c) and (d) of the module
    ledger, plus one bit for the sign. **No ``N`` appears**: this is the datasheet
    line. Verification of an instance of height ``L`` at this budget cannot refuse
    for width, at any number of stages.
    """
    return max(9 * height_bits + 7, (4 * dim + 1) * height_bits + dim + 2) + 1


# --- checked fixed-width integer arithmetic ----------------------------------
def checked_add(a: int, b: int, limit: int) -> tuple[int, int]:
    """``a + b`` in the budget. ``(_OK, sum)`` or ``(_OVERFLOW, 0)``.

    On refusal the value is 0 and must not be consumed. Both operands are
    re-checked, so an out-of-budget input cannot be smuggled through a chain.
    """
    if not fits(a, limit) or not fits(b, limit):
        return _OVERFLOW, 0
    s = a + b
    if not fits(s, limit):
        return _OVERFLOW, 0
    return _OK, s


def checked_sub(a: int, b: int, limit: int) -> tuple[int, int]:
    """``a - b`` in the budget. ``(_OK, difference)`` or ``(_OVERFLOW, 0)``."""
    if not fits(a, limit) or not fits(b, limit):
        return _OVERFLOW, 0
    s = a - b
    if not fits(s, limit):
        return _OVERFLOW, 0
    return _OK, s


def checked_mul(a: int, b: int, limit: int) -> tuple[int, int]:
    """``a * b`` in the budget. ``(_OK, product)`` or ``(_OVERFLOW, 0)``.

    This is the load-bearing helper: a wrapped product would silently turn a false
    certificate into an accepted one, so the product is formed at full precision and
    the budget is checked before it is handed back.
    """
    if not fits(a, limit) or not fits(b, limit):
        return _OVERFLOW, 0
    p = a * b
    if not fits(p, limit):
        return _OVERFLOW, 0
    return _OK, p


# --- unnormalized pair rationals ---------------------------------------------
def rat_sign(a: Rat) -> int:
    """``-1 / 0 / +1``. Denominators are positive, so this is the numerator's sign
    and costs no arithmetic at all (ledger row O1)."""
    return (a[0] > 0) - (a[0] < 0)


def rat_mul(a: Rat, b: Rat, limit: int) -> tuple[int, Rat]:
    """``a * b`` as ``(n1*n2, q1*q2)`` -- no reduction (ledger: widths add)."""
    st, n = checked_mul(a[0], b[0], limit)
    if st != _OK:
        return st, _ZERO
    st, q = checked_mul(a[1], b[1], limit)
    if st != _OK:
        return st, _ZERO
    return _OK, (n, q)


def rat_add(a: Rat, b: Rat, limit: int) -> tuple[int, Rat]:
    """``a + b`` as ``(n1*q2 + n2*q1, q1*q2)`` -- no reduction, no gcd."""
    st, left = checked_mul(a[0], b[1], limit)
    if st != _OK:
        return st, _ZERO
    st, right = checked_mul(b[0], a[1], limit)
    if st != _OK:
        return st, _ZERO
    st, num = checked_add(left, right, limit)
    if st != _OK:
        return st, _ZERO
    st, den = checked_mul(a[1], b[1], limit)
    if st != _OK:
        return st, _ZERO
    return _OK, (num, den)


def rat_sub(a: Rat, b: Rat, limit: int) -> tuple[int, Rat]:
    """``a - b`` as ``(n1*q2 - n2*q1, q1*q2)`` -- no reduction, no gcd."""
    st, left = checked_mul(a[0], b[1], limit)
    if st != _OK:
        return st, _ZERO
    st, right = checked_mul(b[0], a[1], limit)
    if st != _OK:
        return st, _ZERO
    st, num = checked_sub(left, right, limit)
    if st != _OK:
        return st, _ZERO
    st, den = checked_mul(a[1], b[1], limit)
    if st != _OK:
        return st, _ZERO
    return _OK, (num, den)


def rat_cmp(a: Rat, b: Rat, limit: int) -> tuple[int, int]:
    """``sign(a - b)`` by cross-multiplication (ledger row O3').

    Costs ``b(num a) + b(den b)`` and ``b(num b) + b(den a)`` on the two products and
    one more bit on their difference -- ``2L + 1`` at instance height ``L``. No
    denominator product is formed, which is why the dominance leg is the cheapest
    check in the tier.
    """
    st, left = checked_mul(a[0], b[1], limit)
    if st != _OK:
        return st, 0
    st, right = checked_mul(b[0], a[1], limit)
    if st != _OK:
        return st, 0
    st, diff = checked_sub(left, right, limit)
    if st != _OK:
        return st, 0
    return _OK, (diff > 0) - (diff < 0)


# --- instance and verdict -----------------------------------------------------
@dataclass(frozen=True)
class IdentityCoupledInstance:
    """Instance data plus the emitted certificate, all as ``(num, den)`` pairs.

    ``N = len(w)`` stages of dimension ``d = len(x_bar[0])``; one keep-out per stage.

    * ``w``      -- stage goal weights ``w_k`` (the ``W_k = w_k I`` diagonal mass).
    * ``lam``    -- per-stage multiplier loads ``lam_k`` (the certificate's duals).
    * ``g``      -- goal points ``g_k``, ``N x d``.
    * ``c``      -- keep-out centres ``c_k``, ``N x d``.
    * ``rho_sq`` -- SQUARED keep-out radii ``rho_k^2`` (squared, so the checker never
      forms a product it does not need -- ledger row O5 counts ``rho^2`` as data).
    * ``x_bar``  -- the emitted witness, ``N x d``.

    Every field is checked structurally before any arithmetic runs; a denominator
    that is not strictly positive is a ``bad-input`` refusal, not an assumption.
    """

    w: tuple[Rat, ...]
    lam: tuple[Rat, ...]
    g: tuple[tuple[Rat, ...], ...]
    c: tuple[tuple[Rat, ...], ...]
    rho_sq: tuple[Rat, ...]
    x_bar: tuple[tuple[Rat, ...], ...]


@dataclass(frozen=True)
class Verdict:
    """The decision-layer verdict.

    * ``accepted`` -- (a)-(d) held: ``M(lam, t) >= 0`` at ``t = psi(x_bar)``, hence
      ``t <= J* <= f_0(x_bar)``.
    * ``closed``   -- (e) held too: ``t = f_0(x_bar)``, the bracket closes and
      ``x_bar`` is a certified global optimum. Only meaningful when accepted.
    * ``reason``   -- ``REASON_*`` on refusal, ``REASON_OK`` on acceptance.
    * ``stage``    -- the stage the refusal is scoped to, ``-1`` if none.
    * ``operand``  -- the operand or field named by the refusal, ``""`` if none.
    """

    accepted: bool
    closed: bool
    reason: str
    stage: int
    operand: str

    def __str__(self) -> str:
        if self.accepted:
            return "ACCEPT (bracket closed)" if self.closed else "ACCEPT (bracket open)"
        where = f"@{self.operand}" if self.operand else ""
        at = f"[stage {self.stage}]" if self.stage >= 0 else ""
        return f"REFUSE {self.reason}{where}{at}"


def _refuse(reason: str, stage: int, operand: str) -> Verdict:
    return Verdict(accepted=False, closed=False, reason=reason, stage=stage,
                   operand=operand)


# --- input screening ----------------------------------------------------------
def _is_rat(value: object) -> bool:
    """A well-formed unnormalized rational: an ``(int, int)`` pair, denominator
    strictly positive. Booleans are rejected -- ``bool`` is an ``int`` subclass in
    Python but has no place in certificate data. Floats are rejected by the same
    test, keeping the certified path exactly integer."""
    return (isinstance(value, tuple) and len(value) == 2
            and isinstance(value[0], int) and not isinstance(value[0], bool)
            and isinstance(value[1], int) and not isinstance(value[1], bool)
            and value[1] > 0)


def _shape_fault(inst: IdentityCoupledInstance) -> str:
    """Structural screen. Returns ``""`` or the malformed field's name."""
    n = len(inst.w)
    if n < 1:
        return "stage-count"
    if len(inst.lam) != n or len(inst.rho_sq) != n:
        return "stage-vector-length"
    if len(inst.g) != n or len(inst.c) != n or len(inst.x_bar) != n:
        return "stage-block-count"
    dim = len(inst.x_bar[0])
    if dim < 1:
        return "stage-dimension"
    for k in range(n):
        if (len(inst.g[k]) != dim or len(inst.c[k]) != dim
                or len(inst.x_bar[k]) != dim):
            return "block-dimension"
    return ""


def _rational_fault(inst: IdentityCoupledInstance) -> str:
    """Representation screen. Returns ``""`` or the offending datum's name."""
    for name, vec in (("w", inst.w), ("lam", inst.lam), ("rho_sq", inst.rho_sq)):
        for k in range(len(vec)):
            if not _is_rat(vec[k]):
                return f"{name}[{k}]"
    for name, blk in (("g", inst.g), ("c", inst.c), ("x_bar", inst.x_bar)):
        for k in range(len(blk)):
            for i in range(len(blk[k])):
                if not _is_rat(blk[k][i]):
                    return f"{name}[{k}][{i}]"
    return ""


def _input_width_fault(inst: IdentityCoupledInstance, limit: int) -> str:
    """Budget screen on the DATA (the ledger's ``L``): every numerator and every
    denominator must already fit, before a single operation runs. Returns ``""`` or
    the offending datum's name."""
    for name, vec in (("w", inst.w), ("lam", inst.lam), ("rho_sq", inst.rho_sq)):
        for k in range(len(vec)):
            if not fits(vec[k][0], limit) or not fits(vec[k][1], limit):
                return f"{name}[{k}]"
    for name, blk in (("g", inst.g), ("c", inst.c), ("x_bar", inst.x_bar)):
        for k in range(len(blk)):
            for i in range(len(blk[k])):
                if not fits(blk[k][i][0], limit) or not fits(blk[k][i][1], limit):
                    return f"{name}[{k}][{i}]"
    return ""


def instance_height_bits(inst: IdentityCoupledInstance) -> int:
    """``L``: the largest bit-length of any numerator or denominator in the data.

    The ledger's single input parameter. Feed it to `sufficient_width_bits` together
    with the stage dimension to size the budget before running -- the answer does not
    depend on the number of stages.
    """
    height = 0
    for vec in (inst.w, inst.lam, inst.rho_sq):
        for k in range(len(vec)):
            height = max(height, abs(vec[k][0]).bit_length(), vec[k][1].bit_length())
    for blk in (inst.g, inst.c, inst.x_bar):
        for k in range(len(blk)):
            for i in range(len(blk[k])):
                height = max(height, abs(blk[k][i][0]).bit_length(),
                             blk[k][i][1].bit_length())
    return height


def _screen(inst: IdentityCoupledInstance, width_bits: int) -> Verdict | None:
    """All refusals that precede arithmetic: shape, representation, data width."""
    fault = _shape_fault(inst)
    if fault:
        return _refuse(REASON_BAD_INPUT, -1, fault)
    fault = _rational_fault(inst)
    if fault:
        return _refuse(REASON_BAD_INPUT, -1, fault)
    fault = _input_width_fault(inst, width_limit(width_bits))
    if fault:
        return _refuse(REASON_WIDTH, -1, fault)
    return None


# --- the five checks ----------------------------------------------------------
def _check_multipliers(inst: IdentityCoupledInstance) -> Verdict | None:
    """(a) ``lam_k >= 0``. Denominators are positive after the screen, so this is the
    numerator's sign: ledger row O1, width ``L``, no arithmetic."""
    for k in range(len(inst.lam)):
        if inst.lam[k][0] < 0:
            return _refuse(REASON_MULTIPLIER, k, "lam")
    return None


def _check_dominance(inst: IdentityCoupledInstance, limit: int) -> Verdict | None:
    """(b) weak diagonal dominance of ``H(lam)``, identity-coupled form.

    Row ``i`` of stage ``k`` is ``[-1, m_k + w_k - lam_k, -1]`` (one ``-1`` per
    neighbour), so ``H_ii >= sum_{j != i} |H_ij|`` is exactly ``w_k >= lam_k``. WDD is
    *weak*: equality is in-regime and accepts. Dominance also makes the diagonal
    nonnegative for free (``m_k + (w_k - lam_k) >= 0``), which is the other half of
    the Gershgorin hypothesis, so no separate check is needed. Ledger row O3',
    ``2L + 1`` bits.
    """
    for k in range(len(inst.w)):
        st, sgn = rat_cmp(inst.w[k], inst.lam[k], limit)
        if st != _OK:
            return _refuse(REASON_WIDTH, k, "dominance.cross-product")
        if sgn < 0:
            return _refuse(REASON_MARGIN, k, "w-minus-lam")
    return None


def _diag_coefficient(w_k: Rat, lam_k: Rat, m_k: int,
                      limit: int) -> tuple[int, str, Rat]:
    """``m_k + w_k - lam_k``: the (scalar) diagonal of ``H(lam)`` at stage ``k``.
    Width ``(2L + 3, 2L)``."""
    st, a = rat_sub(w_k, lam_k, limit)
    if st != _OK:
        return st, "stationarity.diagonal", _ZERO
    st, a = rat_add(a, (m_k, 1), limit)
    if st != _OK:
        return st, "stationarity.diagonal", _ZERO
    return _OK, "", a


def _residual(inst: IdentityCoupledInstance, a: Rat, k: int, i: int,
              left: Rat, right: Rat, limit: int) -> tuple[int, str, Rat]:
    """Coordinate ``i`` of the stationarity residual at stage ``k``, i.e. half of
    ``grad(f_0 - sum_j lam_j f_j)(x_bar)``:

        (m_k + w_k - lam_k) x_k[i] - x_{k-1}[i] - x_{k+1}[i]
                                   - w_k g_k[i] + lam_k c_k[i]

    Evaluated as a fixed five-term table (a static table in the C rendering);
    boundary stages pass an exact zero for the missing neighbour, so the shape is the
    same at every stage and no branch changes the operation count. Ledger row O4,
    ``(9L + 7, 9L)``.
    """
    terms = (
        (a, inst.x_bar[k][i], 1),
        (_ONE, left, -1),
        (_ONE, right, -1),
        (inst.w[k], inst.g[k][i], -1),
        (inst.lam[k], inst.c[k][i], 1),
    )
    acc = _ZERO
    for coef, val, sgn in terms:
        st, term = rat_mul(coef, val, limit)
        if st != _OK:
            return st, "stationarity.product", _ZERO
        if sgn > 0:
            st, acc = rat_add(acc, term, limit)
        else:
            st, acc = rat_sub(acc, term, limit)
        if st != _OK:
            return st, "stationarity.accumulate", _ZERO
    return _OK, "", acc


def _check_stationarity(inst: IdentityCoupledInstance, limit: int) -> Verdict | None:
    """(c) the residual must be EXACTLY zero at every stage and coordinate -- a zero
    test on an integer numerator, never a tolerance."""
    n = len(inst.w)
    dim = len(inst.x_bar[0])
    for k in range(n):
        m_k = (1 if k > 0 else 0) + (1 if k < n - 1 else 0)
        st, operand, a = _diag_coefficient(inst.w[k], inst.lam[k], m_k, limit)
        if st != _OK:
            return _refuse(REASON_WIDTH, k, operand)
        for i in range(dim):
            left = inst.x_bar[k - 1][i] if k > 0 else _ZERO
            right = inst.x_bar[k + 1][i] if k < n - 1 else _ZERO
            st, operand, resid = _residual(inst, a, k, i, left, right, limit)
            if st != _OK:
                return _refuse(REASON_WIDTH, k, operand)
            if resid[0] != 0:
                return _refuse(REASON_STATIONARITY, k, "gradient")
    return None


def _feasibility_value(x_k: tuple[Rat, ...], c_k: tuple[Rat, ...], rho_sq_k: Rat,
                       limit: int) -> tuple[int, str, Rat]:
    """``f_k(x_bar) = sum_i (x_k[i] - c_k[i])^2 - rho_sq_k``. Ledger row O5,
    ``((4d+1)L + d + 2, (4d+1)L)``."""
    acc = _ZERO
    for i in range(len(x_k)):
        st, offset = rat_sub(x_k[i], c_k[i], limit)
        if st != _OK:
            return st, "feasibility.centre-offset", _ZERO
        st, square = rat_mul(offset, offset, limit)
        if st != _OK:
            return st, "feasibility.square", _ZERO
        st, acc = rat_add(acc, square, limit)
        if st != _OK:
            return st, "feasibility.accumulate", _ZERO
    st, value = rat_sub(acc, rho_sq_k, limit)
    if st != _OK:
        return st, "feasibility.radius", _ZERO
    return _OK, "", value


def _check_feasibility_and_slackness(
        inst: IdentityCoupledInstance, limit: int) -> tuple[Verdict | None, bool]:
    """(d) ``f_k(x_bar) >= 0`` and (e) ``lam_k = 0 or f_k(x_bar) = 0``.

    (e) is folded into the same stage pass because it is a pair of zero tests on
    operands (d) and (a) have already produced -- no multiplication, and crucially no
    sum over stages, which is exactly why the decision layer carries no ``log N``.
    Returns ``(refusal or None, closed)``.
    """
    closed = True
    for k in range(len(inst.w)):
        st, operand, value = _feasibility_value(
            inst.x_bar[k], inst.c[k], inst.rho_sq[k], limit)
        if st != _OK:
            return _refuse(REASON_WIDTH, k, operand), False
        if rat_sign(value) < 0:
            return _refuse(REASON_FEASIBILITY, k, "keep-out"), False
        if inst.lam[k][0] != 0 and value[0] != 0:
            closed = False          # accepted, but the bracket stays open
    return None, closed


# --- the entry point ----------------------------------------------------------
@contract(width_bits=Interval(8, 4096))
def verify_fixed_width(inst: IdentityCoupledInstance,
                       width_bits: int = 64) -> Verdict:
    """Verify an identity-coupled WDD-tier certificate at a fixed width budget.

    ``width_bits`` is the word width ``W``: every operand on the certified path,
    inputs included, must satisfy ``|x| < 2**(W-1)`` or the verdict is
    ``width-exceeded`` naming the operand and the stage. The default 64 is int64
    semantics. Size it in advance with
    ``sufficient_width_bits(instance_height_bits(inst), d)`` -- a bound in the data
    height and the stage dimension alone, with no dependence on the number of stages.

    Returns a `Verdict`. ACCEPT means ``t <= J* <= f_0(x_bar)`` for
    ``t = psi(x_bar)``; ACCEPT with ``closed`` means additionally
    ``f_0(x_bar) = J*`` -- ``x_bar`` is a certified global optimum. REFUSE asserts
    nothing about the instance: the tier is sound and incomplete, and a refusal for
    ``width-exceeded`` in particular says only that this budget was too small.

    Never raises on certificate data: malformed input is the ``bad-input`` verdict.
    The one declared contract is on ``width_bits`` -- an out-of-range budget is a
    caller programming error, not a certificate defect, and raises `ContractError`
    in the sandbox (the range is rendered as an ACSL precondition under translation).
    """
    early = _screen(inst, width_bits)
    if early is not None:
        return early
    limit = width_limit(width_bits)
    verdict = _check_multipliers(inst)
    if verdict is not None:
        return verdict
    verdict = _check_dominance(inst, limit)
    if verdict is not None:
        return verdict
    verdict = _check_stationarity(inst, limit)
    if verdict is not None:
        return verdict
    verdict, closed = _check_feasibility_and_slackness(inst, limit)
    if verdict is not None:
        return verdict
    return Verdict(accepted=True, closed=closed, reason=REASON_OK, stage=-1,
                   operand="")
