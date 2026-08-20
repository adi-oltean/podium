"""Fixed-width WDD-tier certificate verification -- podium.verify.wdd_fixed.

The trusted claim has two halves and both are pinned here.

*Decision correctness*: on the identity-coupled closing family the checker accepts
with the bracket closed, and each way a certificate can be wrong (negative
multiplier, exhausted WDD margin, nonzero stationarity residual, infeasible witness,
unexploited slack) produces the right refusal at the right stage -- cross-checked
against an independent exact-`Fraction` oracle over a fixed-seed randomized sweep.
Every assertion is EXACT: the tier is a sign-and-zero-test machine, so there is no
tolerance anywhere, not even a derived one.

*Fixed-width discipline*: the budget is the load-bearing semantic. Overflow must be a
loud refusal, never a wraparound -- the wraparound test below builds an instance whose
stationarity residual is exactly ``2**12``, so a 12-bit machine that wrapped would see
a zero residual and ACCEPT a false certificate; the checker refuses instead, and the
same instance at 64 bits is correctly refused on its merits. The horizon-independence
test runs one budget, sized from the data height alone, over horizons 1..97.
"""

import dataclasses
import random
from fractions import Fraction

import pytest

from podium.verify import wdd_fixed as wf
from podium.verify.contracts import ContractError

# Integer points with integer squared norms: (3,4) -> 25, (5,12) -> 169,
# (8,15) -> 289.  Cycling with period 3 keeps the instance height flat in N,
# which is what makes the horizon-independence assertion meaningful.
_PTS = ((3, 4), (5, 12), (8, 15))


def closing_family(n, lam=None, w=None):
    """The identity-coupled closing family, exactly and in-regime by construction.

    Witnesses sit ON each keep-out sphere (``f_k(x_bar) = 0``), multipliers are
    positive, and the goals are recovered from stationarity so the residual is
    exactly zero:

        g_k[i] = ((m_k + w_k - lam_k) x_k[i] - x_{k-1}[i] - x_{k+1}[i]
                  + lam_k c_k[i]) / w_k

    with centres at the origin. ``lam`` and ``w`` are integer sequences (defaults:
    all 1 and all 2, so the WDD margin ``w_k - lam_k`` is exactly 1 at every stage).
    """
    lam = [1] * n if lam is None else list(lam)
    w = [2] * n if w is None else list(w)
    pts = [_PTS[k % 3] for k in range(n)]
    x_bar = tuple(tuple((v, 1) for v in pts[k]) for k in range(n))
    c = tuple(((0, 1), (0, 1)) for _ in range(n))
    rho_sq = tuple((pts[k][0] ** 2 + pts[k][1] ** 2, 1) for k in range(n))
    g = []
    for k in range(n):
        m_k = (1 if k > 0 else 0) + (1 if k < n - 1 else 0)
        row = []
        for i in range(2):
            num = (m_k + w[k] - lam[k]) * pts[k][i]
            if k > 0:
                num -= pts[k - 1][i]
            if k < n - 1:
                num -= pts[k + 1][i]
            row.append((num, w[k]))
        g.append(tuple(row))
    return wf.IdentityCoupledInstance(
        w=tuple((v, 1) for v in w), lam=tuple((v, 1) for v in lam),
        g=tuple(g), c=c, rho_sq=rho_sq, x_bar=x_bar)


def _set_scalar(inst, field, k, value):
    vec = list(getattr(inst, field))
    vec[k] = value
    return dataclasses.replace(inst, **{field: tuple(vec)})


def _set_entry(inst, field, k, i, value):
    blk = [list(row) for row in getattr(inst, field)]
    blk[k][i] = value
    return dataclasses.replace(inst, **{field: tuple(tuple(r) for r in blk)})


# ---- the closing family -------------------------------------------------------

def test_closing_family_accepts_with_closure():
    verdict = wf.verify_fixed_width(closing_family(3))
    assert verdict.accepted is True
    assert verdict.closed is True
    assert verdict.reason == wf.REASON_OK
    assert verdict.stage == -1
    assert verdict.operand == ""
    assert str(verdict) == "ACCEPT (bracket closed)"


def test_closing_family_margin_is_exactly_one():
    # w_k - lam_k = 1 at every stage: strictly inside the tier, by exactly one
    # unit of goal-weight mass.  Shaving that unit lands on the boundary; shaving
    # more leaves the tier.
    inst = closing_family(3)
    limit = wf.width_limit(64)
    for k in range(3):
        status, sign = wf.rat_cmp(inst.w[k], inst.lam[k], limit)
        assert status == 0 and sign == 1
        status, diff = wf.rat_sub(inst.w[k], inst.lam[k], limit)
        assert status == 0 and diff == (1, 1)


def test_horizon_independent_width_budget():
    """The datasheet line: ONE budget, sized from the data height and the stage
    dimension alone, verifies every horizon.  The height itself is flat in N once
    the period-3 data cycle is complete (N >= 3)."""
    height = wf.instance_height_bits(closing_family(3))
    assert height == 9                       # rho^2 = 289 is the tallest datum
    budget = wf.sufficient_width_bits(height, 2)
    for n in (1, 2, 3, 4, 8, 40, 97):
        inst = closing_family(n)
        if n >= 3:
            assert wf.instance_height_bits(inst) == height
        verdict = wf.verify_fixed_width(inst, width_bits=budget)
        assert verdict.accepted is True and verdict.closed is True
        # int64 is more than enough for this family at every horizon
        assert wf.verify_fixed_width(inst, width_bits=64).accepted is True


def test_single_stage_scalar_dimension():
    # N = 1, d = 1: no neighbours, so the diagonal is w_0 - lam_0 alone.
    # a = 2, residual = 2*5 - 3*(10/3) = 0, f = 25 - 25 = 0, lam > 0 -> closed.
    inst = wf.IdentityCoupledInstance(
        w=((3, 1),), lam=((1, 1),), g=(((10, 3),),), c=(((0, 1),),),
        rho_sq=((25, 1),), x_bar=(((5, 1),),))
    verdict = wf.verify_fixed_width(inst)
    assert verdict.accepted is True and verdict.closed is True


# ---- refusals: one per way a certificate can be wrong --------------------------

def test_margin_violation_names_the_stage():
    # lam_1 = 3 against w_1 = 2: the WDD margin at stage 1 is -1.  The goals are
    # rebuilt for the new multiplier, so stationarity still holds exactly -- the
    # ONLY thing wrong is the margin, and the refusal must say so.
    inst = closing_family(3, lam=[1, 3, 1])
    verdict = wf.verify_fixed_width(inst)
    assert verdict.accepted is False
    assert verdict.reason == wf.REASON_MARGIN
    assert verdict.stage == 1
    assert verdict.operand == "w-minus-lam"
    assert str(verdict) == "REFUSE margin-violated@w-minus-lam[stage 1]"
    # stages 0 and 2 are untouched: violating stage 0 instead moves the report
    assert wf.verify_fixed_width(closing_family(3, lam=[3, 1, 1])).stage == 0
    assert wf.verify_fixed_width(closing_family(3, lam=[1, 1, 3])).stage == 2


def test_zero_margin_is_in_regime():
    """WDD is WEAK dominance: ``H_ii >= sum |H_ij|`` with equality is dominant, the
    Gershgorin discs touch the origin, and ``H >= 0`` still follows.  So
    ``w_k == lam_k`` -- zero margin at every stage -- must ACCEPT, and one unit
    less must not."""
    inst = closing_family(3, lam=[1, 1, 1], w=[1, 1, 1])
    verdict = wf.verify_fixed_width(inst)
    assert verdict.accepted is True and verdict.closed is True
    limit = wf.width_limit(64)
    for k in range(3):
        assert wf.rat_cmp(inst.w[k], inst.lam[k], limit) == (0, 0)   # exactly zero
    below = closing_family(3, lam=[2, 2, 2], w=[1, 1, 1])
    assert wf.verify_fixed_width(below).reason == wf.REASON_MARGIN


def test_negative_multiplier_refused_before_anything_else():
    inst = closing_family(3, lam=[1, -1, 1])
    verdict = wf.verify_fixed_width(inst)
    assert verdict.accepted is False
    assert verdict.reason == wf.REASON_MULTIPLIER
    assert verdict.stage == 1
    # a negative numerator with a positive denominator is still negative
    inst2 = _set_scalar(closing_family(3), "lam", 2, (-1, 7))
    assert wf.verify_fixed_width(inst2).reason == wf.REASON_MULTIPLIER
    assert wf.verify_fixed_width(inst2).stage == 2


def test_stationarity_residual_must_be_exactly_zero():
    inst = _set_entry(closing_family(3), "g", 1, 0, (5, 2))   # was (4, 2)
    verdict = wf.verify_fixed_width(inst)
    assert verdict.accepted is False
    assert verdict.reason == wf.REASON_STATIONARITY
    assert verdict.stage == 1
    assert verdict.operand == "gradient"
    # the smallest possible perturbation still refuses: no tolerance exists
    tiny = _set_entry(closing_family(3), "g", 2, 1, (18 * 10 ** 6 + 1, 2 * 10 ** 6))
    assert wf.verify_fixed_width(tiny).reason == wf.REASON_STATIONARITY
    # ... and the unperturbed rewrite of the same number, unreduced, still accepts
    same = _set_entry(closing_family(3), "g", 2, 1, (18 * 10 ** 6, 2 * 10 ** 6))
    assert wf.verify_fixed_width(same).accepted is True


def test_infeasible_witness():
    # rho^2 = 170 against ||x_bar_1||^2 = 169: the witness is INSIDE the keep-out.
    inst = _set_scalar(closing_family(3), "rho_sq", 1, (170, 1))
    verdict = wf.verify_fixed_width(inst)
    assert verdict.accepted is False
    assert verdict.reason == wf.REASON_FEASIBILITY
    assert verdict.stage == 1
    assert verdict.operand == "keep-out"


def test_slack_keeps_the_bracket_open_without_refusing():
    """Complementary slackness is the CLOSURE leg, not an acceptance leg: with
    ``lam_1 > 0`` and ``f_1(x_bar) = 1 > 0`` the bound ``t <= J* <= f_0(x_bar)``
    still holds, so the verdict is ACCEPT with the bracket open."""
    inst = _set_scalar(closing_family(3), "rho_sq", 1, (168, 1))
    verdict = wf.verify_fixed_width(inst)
    assert verdict.accepted is True
    assert verdict.closed is False
    assert str(verdict) == "ACCEPT (bracket open)"
    # lam_1 = 0 makes the same slack harmless: the product lam_1 f_1 vanishes.
    # The goals must follow the multiplier change (it moves H's diagonal), so the
    # family is rebuilt rather than patched.
    closed = _set_scalar(closing_family(3, lam=[1, 0, 1]), "rho_sq", 1, (168, 1))
    verdict2 = wf.verify_fixed_width(closed)
    assert verdict2.accepted is True and verdict2.closed is True


# ---- the fixed-width contract --------------------------------------------------

def test_no_silent_wraparound_a_wrapping_machine_would_accept():
    """The instance is built so that a 12-bit machine WITH wraparound accepts a
    false certificate, and the checker must not.

    ``N = 1``, ``d = 1``, ``lam = 0``, ``w = 64``, ``x_bar = 64``, ``g = c = 0``:
    the stationarity residual is ``64 * 64 = 4096 = 2**12``, which is exactly
    ``0`` in 12-bit two's complement, and the feasibility value ``64**2 - 0`` wraps
    to ``0`` as well -- so a wrapping verifier would report ACCEPT with the bracket
    closed.  The truth is that the residual is nonzero and the certificate is
    garbage.
    """
    inst = wf.IdentityCoupledInstance(
        w=((64, 1),), lam=((0, 1),), g=(((0, 1),),), c=(((0, 1),),),
        rho_sq=((0, 1),), x_bar=(((64, 1),),))
    # the two quantities a wrapping machine would see as zero
    assert (64 * 64) % (1 << 12) == 0

    refused = wf.verify_fixed_width(inst, width_bits=12)
    assert refused.accepted is False
    assert refused.reason == wf.REASON_WIDTH
    assert refused.stage == 0
    assert refused.operand == "stationarity.product"

    # at a budget that can hold the residual, the same instance is refused on its
    # merits -- proving the 12-bit refusal was the guard firing, not the verdict
    honest = wf.verify_fixed_width(inst, width_bits=64)
    assert honest.accepted is False
    assert honest.reason == wf.REASON_STATIONARITY
    assert honest.stage == 0


def test_width_budget_refuses_then_accepts_the_same_instance():
    """Refusal for width is a statement about the BUDGET, not the instance.

    ``d = 2``, one stage, witness ``(45, 45)`` on a keep-out of squared radius
    2000.  Each squared coordinate (2025) fits a 12-bit word; their sum (4050) does
    not, so the accumulate step refuses.  At 16 bits the identical instance
    verifies and closes.
    """
    inst = wf.IdentityCoupledInstance(
        w=((2, 1),), lam=((0, 1),), g=(((90, 2), (90, 2)),),
        c=(((0, 1), (0, 1)),), rho_sq=((2000, 1),),
        x_bar=(((45, 1), (45, 1)),))
    refused = wf.verify_fixed_width(inst, width_bits=12)
    assert refused.accepted is False
    assert refused.reason == wf.REASON_WIDTH
    assert refused.stage == 0
    assert refused.operand == "feasibility.accumulate"
    assert str(refused) == "REFUSE width-exceeded@feasibility.accumulate[stage 0]"

    accepted = wf.verify_fixed_width(inst, width_bits=16)
    assert accepted.accepted is True and accepted.closed is True


def test_width_refusal_from_the_input_screen():
    # The data itself must fit before a single operation runs: rho^2 = 169 does not
    # fit an 8-bit budget (limit 128), so the screen refuses with the field named
    # and no stage attached.
    verdict = wf.verify_fixed_width(closing_family(3), width_bits=8)
    assert verdict.accepted is False
    assert verdict.reason == wf.REASON_WIDTH
    assert verdict.operand == "rho_sq[1]"
    assert verdict.stage == -1
    assert wf.verify_fixed_width(closing_family(3), width_bits=64).accepted is True


def test_checked_mul_never_wraps():
    """Exhaustive over a small budget: the helper either returns the EXACT product
    or refuses with value 0.  It never returns the wrapped word, and the sweep
    confirms that wrapping would in fact have produced a different number."""
    bits = 6
    limit = wf.width_limit(bits)          # 32
    modulus = 1 << bits
    wraps_seen = 0
    for a in range(-40, 41):
        for b in range(-40, 41):
            status, value = wf.checked_mul(a, b, limit)
            exact = a * b
            representable = (wf.fits(a, limit) and wf.fits(b, limit)
                             and wf.fits(exact, limit))
            if representable:
                assert status == 0 and value == exact
            else:
                assert status != 0 and value == 0
                wrapped = (exact + limit) % modulus - limit
                wraps_seen += wrapped != exact
    assert wraps_seen > 0                 # the danger the refusal averts is real


def test_checked_add_and_sub_never_wrap():
    limit = wf.width_limit(6)
    for a in range(-40, 41):
        for b in range(-40, 41):
            for op, exact in ((wf.checked_add, a + b), (wf.checked_sub, a - b)):
                status, value = op(a, b, limit)
                if wf.fits(a, limit) and wf.fits(b, limit) and wf.fits(exact, limit):
                    assert status == 0 and value == exact
                else:
                    assert status != 0 and value == 0


def test_fits_excludes_the_extreme_negative_word():
    limit = wf.width_limit(8)
    assert wf.fits(127, limit) and wf.fits(-127, limit)
    assert not wf.fits(128, limit)
    assert not wf.fits(-128, limit)       # so that negation is always total


# ---- the pair-rational arithmetic ---------------------------------------------

def test_pair_rationals_are_never_normalized():
    """gcd is not fixed-width-friendly and the checks do not need it, so the
    operations must return the raw cross-multiplied pair."""
    limit = wf.width_limit(64)
    assert wf.rat_add((1, 2), (1, 2), limit) == (0, (4, 4))
    assert wf.rat_sub((1, 2), (1, 2), limit) == (0, (0, 4))
    assert wf.rat_mul((2, 4), (3, 6), limit) == (0, (6, 24))
    # comparison is still exact on unreduced operands
    assert wf.rat_cmp((1, 2), (2, 4), limit) == (0, 0)
    assert wf.rat_cmp((1, 2), (1, 3), limit) == (0, 1)
    assert wf.rat_cmp((1, 3), (1, 2), limit) == (0, -1)
    assert wf.rat_sign((-3, 7)) == -1 and wf.rat_sign((0, 7)) == 0


def test_rational_helpers_propagate_overflow_as_status():
    limit = wf.width_limit(8)             # 128
    assert wf.rat_mul((20, 1), (20, 1), limit) == (1, (0, 1))
    assert wf.rat_add((100, 3), (100, 3), limit) == (1, (0, 1))
    assert wf.rat_cmp((100, 3), (1, 100), limit) == (1, 0)   # 100 * 100 overflows
    assert wf.rat_cmp((100, 3), (1, 1), limit) == (0, 1)     # 100 - 3 does not


def test_sufficient_width_bits_is_the_ledger():
    # max(9L + 7, (4d+1)L + d + 2) + 1, the two tallest ledger rows plus a sign bit
    assert wf.sufficient_width_bits(1, 1) == max(16, 8) + 1
    assert wf.sufficient_width_bits(10, 4) == max(97, 176) + 1
    # monotone in both arguments, and free of any horizon term by construction
    assert wf.sufficient_width_bits(11, 4) > wf.sufficient_width_bits(10, 4)
    assert wf.sufficient_width_bits(10, 5) > wf.sufficient_width_bits(10, 4)


def test_instance_height_bits():
    inst = closing_family(3)
    assert wf.instance_height_bits(inst) == (289).bit_length()
    # a taller datum raises the height, and only the bit-length matters
    taller = _set_entry(inst, "x_bar", 0, 0, (3, 1024))
    assert wf.instance_height_bits(taller) == 11


# ---- malformed input -----------------------------------------------------------

def test_denominator_must_be_strictly_positive():
    for bad in ((1, 0), (1, -2)):
        inst = _set_scalar(closing_family(3), "w", 1, bad)
        verdict = wf.verify_fixed_width(inst)
        assert verdict.accepted is False
        assert verdict.reason == wf.REASON_BAD_INPUT
        assert verdict.operand == "w[1]"
    inst = _set_entry(closing_family(3), "x_bar", 2, 1, (4, 0))
    assert wf.verify_fixed_width(inst).operand == "x_bar[2][1]"


def test_floats_and_bools_are_not_certificate_data():
    inst = _set_scalar(closing_family(3), "rho_sq", 0, (25.0, 1))
    assert wf.verify_fixed_width(inst).reason == wf.REASON_BAD_INPUT
    inst = _set_entry(closing_family(3), "g", 0, 0, (True, 1))
    assert wf.verify_fixed_width(inst).reason == wf.REASON_BAD_INPUT
    inst = _set_scalar(closing_family(3), "lam", 0, (1,))
    assert wf.verify_fixed_width(inst).reason == wf.REASON_BAD_INPUT


def test_shape_faults():
    base = closing_family(3)
    empty = dataclasses.replace(base, w=(), lam=(), rho_sq=(), g=(), c=(), x_bar=())
    assert wf.verify_fixed_width(empty).operand == "stage-count"
    short = dataclasses.replace(base, lam=base.lam[:2])
    assert wf.verify_fixed_width(short).operand == "stage-vector-length"
    fewer = dataclasses.replace(base, g=base.g[:2])
    assert wf.verify_fixed_width(fewer).operand == "stage-block-count"
    flat = dataclasses.replace(base, x_bar=((),) * 3)
    assert wf.verify_fixed_width(flat).operand == "stage-dimension"
    ragged = dataclasses.replace(base, c=(base.c[0], ((0, 1),), base.c[2]))
    assert wf.verify_fixed_width(ragged).operand == "block-dimension"


def test_width_bits_carries_a_range_contract():
    inst = closing_family(3)
    with pytest.raises(ContractError, match="width_bits"):
        wf.verify_fixed_width(inst, width_bits=4)
    with pytest.raises(ContractError, match="width_bits"):
        wf.verify_fixed_width(inst, width_bits=8192)
    assert getattr(wf.verify_fixed_width, "__podium_contract__")["width_bits"].lo == 8


# ---- soundness against an independent exact-Fraction oracle --------------------

def _fr(pair):
    return Fraction(pair[0], pair[1])


def _oracle_stationarity(inst):
    """Stage of the first nonzero exact residual, or -1."""
    n = len(inst.w)
    dim = len(inst.x_bar[0])
    for k in range(n):
        m_k = (1 if k > 0 else 0) + (1 if k < n - 1 else 0)
        a = m_k + _fr(inst.w[k]) - _fr(inst.lam[k])
        for i in range(dim):
            resid = a * _fr(inst.x_bar[k][i])
            if k > 0:
                resid -= _fr(inst.x_bar[k - 1][i])
            if k < n - 1:
                resid -= _fr(inst.x_bar[k + 1][i])
            resid -= _fr(inst.w[k]) * _fr(inst.g[k][i])
            resid += _fr(inst.lam[k]) * _fr(inst.c[k][i])
            if resid != 0:
                return k
    return -1


def _oracle_feasibility(inst):
    """(stage of the first violated keep-out or -1, complementary-slackness flag)."""
    closed = True
    for k in range(len(inst.w)):
        value = -_fr(inst.rho_sq[k])
        for i in range(len(inst.x_bar[k])):
            off = _fr(inst.x_bar[k][i]) - _fr(inst.c[k][i])
            value += off * off
        if value < 0:
            return k, False
        if _fr(inst.lam[k]) != 0 and value != 0:
            closed = False
    return -1, closed


def _oracle(inst):
    """The five checks re-implemented in `fractions.Fraction`, in the protocol's
    order, deliberately sharing nothing with the module under test."""
    for k in range(len(inst.lam)):
        if _fr(inst.lam[k]) < 0:
            return (wf.REASON_MULTIPLIER, k, False)
    for k in range(len(inst.w)):
        if _fr(inst.w[k]) < _fr(inst.lam[k]):
            return (wf.REASON_MARGIN, k, False)
    stage = _oracle_stationarity(inst)
    if stage >= 0:
        return (wf.REASON_STATIONARITY, stage, False)
    stage, closed = _oracle_feasibility(inst)
    if stage >= 0:
        return (wf.REASON_FEASIBILITY, stage, False)
    return (wf.REASON_OK, -1, closed)


def _perturb(rng, inst):
    """One random defect (or none) on top of a closing instance."""
    kind = rng.randrange(7)
    k = rng.randrange(len(inst.w))
    i = rng.randrange(2)
    if kind == 1:
        out = _set_scalar(inst, "lam", k, (-rng.randint(1, 4), rng.randint(1, 4)))
    elif kind == 2:
        out = _set_scalar(inst, "w", k, (rng.randint(0, 1), 1))
    elif kind == 3:
        out = _set_entry(inst, "g", k, i, (rng.randint(-9, 9), rng.randint(1, 4)))
    elif kind == 4:
        out = _set_scalar(inst, "rho_sq", k, (rng.randint(0, 400), 1))
    elif kind == 5:
        out = _set_entry(inst, "x_bar", k, i, (rng.randint(-9, 9), rng.randint(1, 3)))
    elif kind == 6:
        out = _set_scalar(inst, "lam", k, (0, 1))
    else:
        out = inst
    return out


def test_matches_an_exact_fraction_oracle_over_random_instances():
    rng = random.Random(20260819)
    seen = {}
    for _ in range(600):
        n = rng.randint(1, 5)
        lam = [rng.randint(0, 3) for _ in range(n)]
        w = [lam[k] + rng.randint(0, 3) for k in range(n)]
        w = [v if v > 0 else 1 for v in w]
        inst = _perturb(rng, closing_family(n, lam=lam, w=w))
        reason, stage, closed = _oracle(inst)
        # 512 bits is far above the ledger bound for this data, so the fixed-width
        # guard never fires and the two verdicts are directly comparable.
        verdict = wf.verify_fixed_width(inst, width_bits=512)
        assert verdict.reason == reason, (reason, str(verdict))
        if reason == wf.REASON_OK:
            assert verdict.accepted is True
            assert verdict.closed is closed
        else:
            assert verdict.accepted is False
            assert verdict.stage == stage
        seen[reason] = seen.get(reason, 0) + 1
    # every branch of the protocol is genuinely exercised by the sweep
    for reason in (wf.REASON_OK, wf.REASON_MULTIPLIER, wf.REASON_MARGIN,
                   wf.REASON_STATIONARITY, wf.REASON_FEASIBILITY):
        assert seen.get(reason, 0) > 0, (reason, seen)
    assert seen[wf.REASON_OK] > 50


def test_random_instances_verify_at_the_ledger_budget():
    """The width theorem, exercised: sizing the budget from the instance height
    alone is enough -- no random instance ever refuses for width at
    `sufficient_width_bits`, at any horizon in the sweep."""
    rng = random.Random(20260820)
    for _ in range(200):
        n = rng.randint(1, 24)
        lam = [rng.randint(0, 3) for _ in range(n)]
        w = [max(1, lam[k] + rng.randint(0, 3)) for k in range(n)]
        inst = _perturb(rng, closing_family(n, lam=lam, w=w))
        budget = wf.sufficient_width_bits(wf.instance_height_bits(inst), 2)
        verdict = wf.verify_fixed_width(inst, width_bits=budget)
        assert verdict.reason != wf.REASON_WIDTH, str(verdict)


# ---------------------------------------------------------------------------
# Overflow propagation out of every intermediate operation.
#
# Overflow is a loud refusal, never a wraparound -- that is the tier's
# load-bearing semantic, and it has to hold at EVERY intermediate step, not
# just at the ones an in-regime instance happens to reach.  Each test below
# picks operands that fit the budget individually and pushes exactly one
# interior operation over it, so the refusal is attributable to that step.
# ---------------------------------------------------------------------------

# 6 bits: |x| < 32.  Numerators/denominators below stay inside that; their
# products do not.
_W = 6
_LIM = 1 << (_W - 1)


def test_rat_mul_refuses_on_the_numerator_product():
    st, _ = wf.rat_mul((7, 1), (7, 1), _LIM)
    assert st == wf._OVERFLOW


def test_rat_mul_refuses_on_the_denominator_product():
    st, _ = wf.rat_mul((1, 7), (1, 7), _LIM)
    assert st == wf._OVERFLOW


def test_rat_add_refuses_on_each_interior_step():
    # left cross-product n1*q2
    assert wf.rat_add((7, 1), (1, 7), _LIM)[0] == wf._OVERFLOW
    # right cross-product n2*q1
    assert wf.rat_add((1, 7), (7, 1), _LIM)[0] == wf._OVERFLOW
    # the sum of two in-budget cross-products
    assert wf.rat_add((31, 1), (31, 1), _LIM)[0] == wf._OVERFLOW
    # the denominator product, with numerators that stay small
    assert wf.rat_add((1, 7), (1, 7), _LIM)[0] == wf._OVERFLOW


def test_rat_sub_refuses_on_each_interior_step():
    assert wf.rat_sub((7, 1), (1, 7), _LIM)[0] == wf._OVERFLOW
    assert wf.rat_sub((1, 7), (7, 1), _LIM)[0] == wf._OVERFLOW
    assert wf.rat_sub((-16, 1), (16, 1), _LIM)[0] == wf._OVERFLOW
    assert wf.rat_sub((1, 7), (1, 7), _LIM)[0] == wf._OVERFLOW


def test_rat_cmp_refuses_on_each_interior_step():
    assert wf.rat_cmp((7, 1), (1, 7), _LIM)[0] == wf._OVERFLOW
    assert wf.rat_cmp((1, 7), (7, 1), _LIM)[0] == wf._OVERFLOW
    assert wf.rat_cmp((16, 1), (-16, 1), _LIM)[0] == wf._OVERFLOW


def test_input_width_fault_names_the_offending_datum():
    """The budget screen runs on the DATA, before any operation."""
    base = closing_family(2)
    wide = 10**5   # clears every real datum, rejects the injected 10**6
    assert wf._input_width_fault(base, wide) == ""
    for field, name in (("w", "w[0]"), ("lam", "lam[0]"), ("rho_sq", "rho_sq[0]")):
        vec = list(getattr(base, field))
        vec[0] = (10**6, 1)
        inst = dataclasses.replace(base, **{field: tuple(vec)})
        assert wf._input_width_fault(inst, wide) == name
    for field, name in (("g", "g[0][0]"), ("c", "c[0][0]"), ("x_bar", "x_bar[0][0]")):
        blk = [list(row) for row in getattr(base, field)]
        blk[0][0] = (10**6, 1)
        inst = dataclasses.replace(
            base, **{field: tuple(tuple(r) for r in blk)})
        assert wf._input_width_fault(inst, wide) == name


def _one_stage(w, lam):
    """A single-stage instance with trivial geometry, so the only thing that can
    overflow is the w-vs-lam arithmetic the test is aiming at."""
    return wf.IdentityCoupledInstance(
        w=(w,), lam=(lam,), g=(((0, 1),),), c=(((0, 1),),),
        rho_sq=((0, 1),), x_bar=(((0, 1),),))


def test_dominance_refuses_for_width_on_the_cross_product():
    """Every datum fits the budget and the cross-product still does not: the tier
    refuses for WIDTH naming the operand, rather than asserting a margin verdict
    it cannot compute."""
    inst = _one_stage((100, 1), (1, 100))     # 100*100 = 10_000, budget is 128
    verdict = wf.verify_fixed_width(inst, width_bits=8)
    assert verdict.reason == wf.REASON_WIDTH
    assert verdict.operand == "dominance.cross-product"


def test_stationarity_refuses_for_width_on_the_diagonal():
    """Dominance clears (cross-products and their difference all fit) and the
    diagonal's denominator product does not -- the refusal has to come from the
    stationarity leg, naming its own operand."""
    inst = _one_stage((1, 100), (1, 100))     # cross 100, diff 0; den 100*100
    verdict = wf.verify_fixed_width(inst, width_bits=8)
    assert verdict.reason == wf.REASON_WIDTH
    assert verdict.operand == "stationarity.diagonal"


def test_diag_coefficient_refuses_on_each_interior_step():
    # w - lam overflows on the cross-product
    st, operand, _ = wf._diag_coefficient((1, 7), (1, 7), 0, _LIM)
    assert st == wf._OVERFLOW and operand == "stationarity.diagonal"
    # w - lam fits, adding the integer m_k does not
    st, operand, _ = wf._diag_coefficient((20, 1), (2, 1), 20, _LIM)
    assert st == wf._OVERFLOW and operand == "stationarity.diagonal"


def test_residual_refuses_on_product_and_on_accumulate():
    inst = closing_family(2)
    # the very first term's product overflows
    st, operand, _ = wf._residual(inst, (10**6, 1), 0, 0, wf._ZERO, wf._ZERO, _LIM)
    assert st == wf._OVERFLOW and operand == "stationarity.product"
    # products fit, the running sum does not
    st, operand, _ = wf._residual(inst, (1, 1), 0, 0, (7, 1), (1, 7), _LIM)
    assert st == wf._OVERFLOW and operand in (
        "stationarity.product", "stationarity.accumulate")


def test_feasibility_value_refuses_on_each_interior_step():
    big = (10**6, 1)
    # centre offset
    st, operand, _ = wf._feasibility_value((big,), ((1, 1),), (1, 1), _LIM)
    assert st == wf._OVERFLOW and operand == "feasibility.centre-offset"
    # the square of an in-budget offset
    st, operand, _ = wf._feasibility_value(((7, 1),), ((0, 1),), (1, 1), _LIM)
    assert st == wf._OVERFLOW and operand == "feasibility.square"
    # the radius subtraction, with the accumulated sum in budget
    st, operand, _ = wf._feasibility_value(((2, 1),), ((0, 1),), (1, 31), _LIM)
    assert st == wf._OVERFLOW and operand == "feasibility.radius"
