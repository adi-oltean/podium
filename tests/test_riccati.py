"""Exact block-tridiagonal (Riccati) PSD certificate -- podium.verify.riccati.

The trusted claim: `block_tridiag_psd` gives the SAME verdict as the dense
`barrier.is_psd` on the assembled matrix, but exploits the band (O(N d^3)); and
`riccati_storage` returns the arrival-cost / information-form pivot blocks for the
benign positive-definite regime (blocks are eliminated in INCREASING index order, so
the pivots accumulate the cost of the eliminated head x_1..x_{k-1} -- they are not the
backward cost-to-go, which eliminates in decreasing order and gives different numbers).
These tests pin soundness (equivalence to is_psd
over random instances), every edge branch, and the no-float discipline.
"""

import itertools
import random
from fractions import Fraction as F

import pytest

from podium.verify import riccati
from podium.verify.barrier import is_psd


def _sym(block):
    d = len(block)
    return [[block[i][j] + block[j][i] for j in range(d)] for i in range(d)]


def _rand_block(rng, d):
    return [[F(rng.randint(-3, 3)) for _ in range(d)] for _ in range(d)]


def _rand_chain(rng, N, d, shift):
    """Random symmetric block-tridiagonal chain. `shift` added to each diagonal
    block's diagonal pushes it toward (large shift) or away from PSD."""
    diag = []
    for _ in range(N):
        a = _rand_block(rng, d)
        s = _sym(a)
        for i in range(d):
            s[i][i] += shift
        diag.append(s)
    off = [_rand_block(rng, d) for _ in range(N - 1)]
    return diag, off


# ---- soundness: equivalence to the dense exact PSD test ---------------------

def test_matches_is_psd_over_random_instances():
    rng = random.Random(20260711)
    seen_true = seen_false = 0
    for _ in range(200):
        N = rng.randint(1, 6)
        d = rng.randint(1, 4)
        shift = rng.choice([F(0), F(3), F(8), F(-2), F(20)])
        diag, off = _rand_chain(rng, N, d, shift)
        dense = riccati.assemble(diag, off)
        expect = is_psd(dense)
        assert riccati.block_tridiag_psd(diag, off) == expect
        seen_true += expect
        seen_false += not expect
    assert seen_true > 0 and seen_false > 0        # both verdicts exercised


def test_assemble_is_symmetric_banded():
    diag = [[[F(2), F(0)], [F(0), F(2)]], [[F(3), F(0)], [F(0), F(3)]]]
    off = [[[F(1), F(0)], [F(0), F(1)]]]
    m = riccati.assemble(diag, off)
    assert m[0][2] == F(1) and m[2][0] == F(1)     # (0,1) block symmetric
    assert m[0][3] == F(0)
    for i in range(4):
        for j in range(4):
            assert m[i][j] == m[j][i]


# ---- positive-definite benign chain -----------------------------------------

def test_benign_pd_chain_certifies_and_storage():
    # strongly diagonally dominant -> PD -> PSD, storage pivots all PSD
    diag = [[[F(5), F(1)], [F(1), F(5)]] for _ in range(4)]
    off = [[[F(-1), F(0)], [F(0), F(-1)]] for _ in range(3)]
    assert riccati.block_tridiag_psd(diag, off) is True
    pivots = riccati.riccati_storage(diag, off)
    assert len(pivots) == 4
    assert all(is_psd(s) for s in pivots)
    # first pivot equals the first diagonal block
    assert pivots[0] == diag[0]


def test_scalar_chain_matches_manual():
    # tridiag(-1, c, -1), c = 5/2, N=3 : pivots 5/2, 21/10, ...
    diag = [[[F(5, 2)]] for _ in range(3)]
    off = [[[F(-1)]] for _ in range(2)]
    pivots = riccati.riccati_storage(diag, off)
    assert pivots[0] == [[F(5, 2)]]
    assert pivots[1] == [[F(5, 2) - F(1) / F(5, 2)]]      # 21/10
    assert riccati.block_tridiag_psd(diag, off) is True


# ---- non-PSD detection ------------------------------------------------------

def test_indefinite_chain_refused():
    diag = [[[F(-1)]], [[F(1)]]]           # negative pivot
    off = [[[F(0)]]]
    assert riccati.block_tridiag_psd(diag, off) is False
    with pytest.raises(ValueError, match="not PSD"):
        riccati.riccati_storage(diag, off)


def test_zero_pivot_psd_and_nonpsd():
    # PSD with a zero pivot: [[0,0],[0,1]] block, no coupling -> still PSD
    diag = [[[F(0), F(0)], [F(0), F(1)]]]
    off = []
    assert riccati.block_tridiag_psd(diag, off) is True
    # zero pivot with a nonzero in its row -> NOT psd
    diag2 = [[[F(0), F(2)], [F(2), F(1)]]]
    assert riccati.block_tridiag_psd(diag2, off) is False
    assert riccati.block_tridiag_psd(diag2, off) == is_psd(riccati.assemble(diag2, off))


def test_storage_singular_pivot_raises():
    # first pivot PSD but singular -> next Schur step cannot invert
    diag = [[[F(0), F(0)], [F(0), F(0)]], [[F(1), F(0)], [F(0), F(1)]]]
    off = [[[F(1), F(0)], [F(0), F(1)]]]
    with pytest.raises(ValueError, match="singular"):
        riccati.riccati_storage(diag, off)


# ---- input validation & no-float discipline ---------------------------------

def test_rejects_floats():
    diag = [[[1.0]]]        # float, not Fraction
    with pytest.raises(TypeError, match="Fraction"):
        riccati.block_tridiag_psd(diag, [])
    with pytest.raises(TypeError, match="Fraction"):
        riccati.block_tridiag_psd([[[F(1)]], [[F(1)]]], [[[2.0]]])   # float in off


def test_shape_validation():
    with pytest.raises(ValueError, match="at least one"):
        riccati.block_tridiag_psd([], [])
    with pytest.raises(ValueError, match="block dimension"):
        riccati.block_tridiag_psd([[]], [])
    with pytest.raises(ValueError, match="off-diagonal block count"):
        riccati.block_tridiag_psd([[[F(1)]], [[F(1)]]], [])          # N=2, off=0
    with pytest.raises(ValueError, match="square and uniform"):
        riccati.block_tridiag_psd([[[F(1), F(0)]]], [])              # non-square
    with pytest.raises(ValueError, match="match the block size"):
        riccati.block_tridiag_psd([[[F(1)]], [[F(1)]]], [[[F(1), F(0)]]])


def test_rejects_asymmetric_diagonal_block():
    # Asymmetric diagonal block: the band sweep's zero-pivot logic and border
    # carry assume symmetry, and barrier.is_psd refuses such input -- the band
    # checkers must refuse rather than silently diverge from the dense verdict.
    asym = [[[F(1), F(5)], [F(0), F(1)]]]
    with pytest.raises(ValueError, match="symmetric"):
        riccati.block_tridiag_psd(asym, [])
    with pytest.raises(ValueError, match="symmetric"):
        riccati.border_band_psd(asym, [], [F(0), F(0)], F(0))
    # the symmetric counterpart is accepted (verdict matches the dense checker)
    sym = [[[F(1), F(5)], [F(5), F(1)]]]
    assert riccati.block_tridiag_psd(sym, []) == is_psd(riccati.assemble(sym, []))
    assert riccati.border_band_psd(sym, [], [F(0), F(0)], F(0)) is False


def test_single_block():
    assert riccati.block_tridiag_psd([[[F(4)]]], []) is True
    assert riccati.riccati_storage([[[F(4)]]], []) == [[[F(4)]]]


# ---- bordered sweep: the full S-procedure LMI M(lam, t) ---------------------

def _assemble_bordered(diag, off, lin, corner):
    m = riccati.assemble(diag, off)
    n = len(m)
    full = [m[i][:] + [lin[i] / F(2)] for i in range(n)]
    full.append([lin[j] / F(2) for j in range(n)] + [corner])
    return full


def test_border_matches_dense_M_over_random_instances():
    rng = random.Random(20260811)
    seen_true = seen_false = 0
    for _ in range(200):
        N = rng.randint(1, 5)
        d = rng.randint(1, 3)
        shift = rng.choice([F(0), F(2), F(6), F(-1), F(12)])
        diag, off = _rand_chain(rng, N, d, shift)
        n = N * d
        lin = [F(rng.randint(-3, 3)) for _ in range(n)]
        corner = F(rng.randint(-2, 12))
        expect = is_psd(_assemble_bordered(diag, off, lin, corner))
        assert riccati.border_band_psd(diag, off, lin, corner) == expect
        seen_true += expect
        seen_false += not expect
    assert seen_true > 0 and seen_false > 0


def test_border_singular_H_range_condition():
    # H singular, lin outside range(H) -> refused regardless of corner
    assert riccati.border_band_psd([[[F(0)]]], [], [F(1)], F(5)) is False
    # lin in range (zero) -> verdict is the corner sign
    assert riccati.border_band_psd([[[F(0)]]], [], [F(0)], F(0)) is True
    assert riccati.border_band_psd([[[F(0)]]], [], [F(0)], F(-1)) is False
    # block case: kernel direction with a consistent (zero) border entry
    diag = [[[F(0), F(0)], [F(0), F(2)]]]
    assert riccati.border_band_psd(diag, [], [F(0), F(2)], F(1)) is True
    assert riccati.border_band_psd(diag, [], [F(1), F(2)], F(1)) is False
    # band-propagated refusal: zero pivot at k>0 whose border entry is nonzero
    # only AFTER the k=0 pivot's carried update (original entry is 0). H=[[1,1],
    # [1,1]] is PSD rank 1; lin/2=[2,0] is outside range(H)=span{(1,1)}: refuse
    # for any corner.
    diag2 = [[[F(1)]], [[F(1)]]]
    off2 = [[[F(1)]]]
    assert riccati.border_band_psd(diag2, off2, [F(4), F(0)], F(5)) is False
    # same chain, in-range border u=-2Hy (y=(2,0)): originally nonzero border
    # entry at the zero pivot must cancel exactly by propagation; accepted,
    # terminal pivot = 1.
    assert riccati.border_band_psd(diag2, off2, [F(-4), F(-4)], F(5)) is True
    assert riccati.border_terminal_pivot(diag2, off2, [F(-4), F(-4)], F(5)) == F(1)
    # two zero pivots: the FIRST consistent, the SECOND (band-propagated) not.
    # A "singularity already detected" shortcut that stops checking the border
    # after one consistent zero pivot must still refuse: H = diag(0)+[[1,1],[1,1]]
    # is PSD rank 1, lin/2 = (0,0,1) outside range(H) = span{(0,1,1)}.
    diag3 = [[[F(0)]], [[F(1)]], [[F(1)]]]
    off3 = [[[F(0)]], [[F(1)]]]
    assert riccati.border_band_psd(diag3, off3, [F(0), F(0), F(2)], F(5)) is False
    # same chain, in-range border lin/2 = H w, w=(0,2,0): BOTH zero pivots must be
    # found consistent (the second only after the k=1 carry); terminal pivot 1.
    assert riccati.border_band_psd(diag3, off3, [F(0), F(4), F(4)], F(5)) is True
    assert riccati.border_terminal_pivot(diag3, off3, [F(0), F(4), F(4)], F(5)) == F(1)
    # two CONSECUTIVE zero pivots (H = 0 block, 2x2): a shortcut that only checks
    # the border once a kernel run is "already established" (i.e. skips the check
    # at the second of two adjacent zero pivots) must still refuse an inconsistent
    # second entry.
    diag4 = [[[F(0)]], [[F(0)]]]
    off4 = [[[F(0)]]]
    assert riccati.border_band_psd(diag4, off4, [F(0), F(2)], F(5)) is False
    assert riccati.border_band_psd(diag4, off4, [F(0), F(0)], F(5)) is True
    assert riccati.border_terminal_pivot(diag4, off4, [F(0), F(0)], F(5)) == F(5)
    assert riccati.border_band_psd(diag4, off4, [F(0), F(0)], F(-1)) is False
    # THREE zero pivots followed by a band-propagated inconsistent zero: a shortcut
    # that stops checking the border after two consistent zero pivots ("kernel run
    # established") must still refuse at the third.
    diag5 = [[[F(0)]], [[F(0)]], [[F(1)]], [[F(1)]]]
    off5 = [[[F(0)]], [[F(0)]], [[F(1)]]]
    assert riccati.border_band_psd(diag5, off5, [F(0), F(0), F(0), F(2)], F(5)) is False
    assert riccati.border_band_psd(diag5, off5, [F(0), F(0), F(4), F(4)], F(5)) is True
    assert riccati.border_terminal_pivot(diag5, off5, [F(0), F(0), F(4), F(4)], F(5)) == F(1)


# ---- exhaustive closure of the zero-pivot check against history-based skips --
#
# The cases above each killed ONE unsound shortcut of the recurring family
# "skip the zero-pivot consistency check (band-row part, border part, or both)
# once the zero-pivot history satisfies some predicate".  The two tests below
# close a BOUNDED slice of that family, not the whole thing -- see the scope
# note after them (round-5 extended the bound; a threshold predicate can
# always be set past whatever bound is tested, so no finite suite closes the
# family outright).
#
# Key reachability fact (itself pinned in test_border_singular_H_range_condition):
# a consistent zero pivot forces its super-diagonal entry to 0, so the next pivot
# receives no update from it.  Hence on d=1 chains the pivot history at any step
# of a non-refused sweep prefix is exactly a string over three symbols --
#   P  positive pivot,
#   Zr raw zero (diagonal entry 0, preceding off 0),
#   Zd derived zero (eliminated to 0; only reachable right after a P),
# and every such string is realizable with entries in {0, 1}.  Enumerating ALL
# of them up to length 6 (in both border-column variants), each followed by a
# zero pivot violating each sub-check separately, means: any implementation that
# skips either sub-check at ANY reachable history of length <= 6 accepts an
# instance the dense oracle refuses.  The small-grid test below verdict-matches
# the dense oracle on every d=1 chain over an INTEGER value grid only -- it does
# not cover d>1 block width or fractional entries; see the round-5 tests further
# down for those.


def _zero_pivot_histories(max_len):
    """All reachable pivot-history strings over {P, Zr, Zd} (Zd only after P)."""
    out, frontier = [()], [()]
    for _ in range(max_len):
        nxt = []
        for s in frontier:
            for sym in ("P", "Zr", "Zd"):
                if sym == "Zd" and (not s or s[-1] != "P"):
                    continue
                nxt.append(s + (sym,))
        out.extend(nxt)
        frontier = nxt
    return out


def _realize_history(s, hot):
    """Realize history string ``s`` as a d=1 chain.  ``hot=False``: border column
    identically zero along the prefix; ``hot=True``: border entry 1 at every P
    pivot (and the exactly-cancelling carry value at every Zd) -- nonzero border
    history, still consistent at every zero pivot."""
    diag, off, lin = [], [], []
    for i, t in enumerate(s):
        if t == "Zd":
            off.append([[F(1)]])            # preceding P has pivot exactly 1
            diag.append([[F(1)]])           # 1 - 1^2/1 = 0 after elimination
            lin.append(F(2) if hot else F(0))   # col = lin/2 - col_prev = 0
        else:
            if i > 0:
                off.append([[F(0)]])
            diag.append([[F(1) if t == "P" else F(0)]])
            lin.append(F(2) if (hot and t == "P") else F(0))
    return diag, off, lin


def _append_violating_zero(prefix, s, prov, mode, tail, carry):
    """Extend a realized history prefix with a zero pivot violating the given
    sub-check(s): ``prov`` raw/derived provenance, ``mode`` border/band/both,
    optional clean tail block.  Returns the extended (diag, off, lin)."""
    pd, po, pl = prefix
    d2, o2, l2 = list(pd), list(po), list(pl)
    if prov == "raw":
        if s:
            o2.append([[F(0)]])
        d2.append([[F(0)]])
    else:
        o2.append([[F(1)]])
        d2.append([[F(1)]])
    want = F(1) if mode in ("border", "both") else F(0)
    l2.append(2 * (want + carry))
    if tail:
        o2.append([[F(1) if mode in ("band", "both") else F(0)]])
        d2.append([[F(1)]])
        l2.append(F(0))
    return d2, o2, l2


def test_border_zero_pivot_check_closure_exhaustive():
    for s in _zero_pivot_histories(6):
        for hot in (False, True):
            prefix = _realize_history(s, hot)
            npos = sum(1 for t in s if t == "P") if hot else 0
            # clean endings: each hot P pivot adds col^2/piv = 1 to the arrival
            # cost, so corner 5+npos keeps the terminal pivot at exactly +5.
            if s:
                for corner, expect in ((F(5 + npos), True), (F(-1), False)):
                    got = riccati.border_band_psd(*prefix, corner)
                    assert got == is_psd(_assemble_bordered(*prefix, corner))
                    assert got == expect, (s, hot, corner)
            # violating zero pivot at position len(s): raw or derived, violating
            # the border entry, the band row, or both, with/without a tail block.
            # Corner 5+npos again, so an implementation that skips the violated
            # sub-check runs to completion and ACCEPTS -- every instance here
            # genuinely distinguishes a skipping mutant from the correct sweep.
            for prov in ("raw", "der"):
                if prov == "der" and not (s and s[-1] == "P"):
                    continue
                carry = F(1) if (hot and prov == "der") else F(0)
                for mode in ("border", "band", "both"):
                    for tail in (True, False):
                        if mode in ("band", "both") and not tail:
                            continue    # band violation needs a next block
                        d2, o2, l2 = _append_violating_zero(
                            prefix, s, prov, mode, tail, carry)
                        corner = F(5 + npos)
                        assert riccati.border_band_psd(d2, o2, l2, corner) is False, (
                            s, hot, prov, mode, tail)
                        assert is_psd(_assemble_bordered(d2, o2, l2, corner)) is False


def test_border_matches_dense_M_exhaustive_small_grid():
    for n in range(1, 4):
        for dv in itertools.product((F(-1), F(0), F(1)), repeat=n):
            diag = [[[x]] for x in dv]
            for ov in itertools.product((F(-1), F(0), F(1)), repeat=n - 1):
                off = [[[x]] for x in ov]
                for lv in itertools.product((F(0), F(2)), repeat=n):
                    lin = list(lv)
                    for corner in (F(-1), F(0), F(5)):
                        assert riccati.border_band_psd(diag, off, lin, corner) == \
                            is_psd(_assemble_bordered(diag, off, lin, corner))


# Round-5 closure extensions (the two tests below).  The exhaustive closure above
# is scoped to d=1, histories of length <= 6, and integer entry values; six unsound
# mutants OUTSIDE that scope survived the then-full suite (each skips or weakens a
# zero-pivot sub-check conditioned on something the closure does not vary):
#   1. skip the border sub-check after >= 7 consistent zero pivots,
#   2. skip the whole zero-pivot check at positions k >= 7,
#   3. skip the band-row sub-check whenever d > 1,
#   4. skip the border sub-check at the last pivot when n >= 8,
#   5. band-row window truncated by one column (misses column k+w) when d > 1,
#   6. band-row check blind to nonzero proper-fraction entries.
# Every instance below is refused by inspection: it has a zero diagonal entry with
# a nonzero entry elsewhere in that row of the bordered M, so x = e_k -+ t e_j
# gives x^T M x = -2 t |m_kj| + t^2 m_jj < 0 for small t > 0.


def test_border_zero_pivot_depth_closure():
    # DEPTH: violating zero pivot at every position up to 16, after straight
    # all-P and all-Zr prefixes, in both border-column variants, with and without
    # a clean tail (so the violation also sits at the LAST pivot of the chain).
    # Kills any skip predicate thresholded on k, n, or zero-pivot count <= 16.
    for sym in ("P", "Zr"):
        for m in range(16):
            s = (sym,) * m
            for hot in (False, True):
                prefix = _realize_history(s, hot)
                npos = sum(1 for t in s if t == "P") if hot else 0
                for prov in ("raw", "der"):
                    if prov == "der" and not (s and s[-1] == "P"):
                        continue
                    carry = F(1) if (hot and prov == "der") else F(0)
                    for mode in ("border", "band"):
                        for tail in (True, False):
                            if mode == "band" and not tail:
                                continue
                            d2, o2, l2 = _append_violating_zero(
                                prefix, s, prov, mode, tail, carry)
                            assert riccati.border_band_psd(
                                d2, o2, l2, F(5 + npos)) is False, (
                                sym, m, hot, prov, mode, tail)


def test_border_zero_pivot_width_and_value_closure():
    # WIDTH x VALUE: for d = 1, 2, 3, a zero pivot at k=0 whose band row is
    # violated at exactly ONE window column j (1 <= j <= w = 2d-1, spanning both
    # the in-block columns and every off-block column), with integer AND
    # proper-fraction violating values.  Everything else benign, border column
    # zero, corner positive: any mutant that skips or narrows the band-row check
    # for that (d, j), or cannot see the fractional value, runs on and accepts.
    for d in (1, 2, 3):
        w = 2 * d - 1
        for j in range(1, w + 1):
            for v in (F(1), F(1, 2), F(-1, 2)):
                diag0 = [[F(0)] * d for _ in range(d)]
                for i in range(1, d):
                    diag0[i][i] = F(1)
                off0 = [[F(0)] * d for _ in range(d)]
                if j < d:
                    diag0[0][j] = v
                    diag0[j][0] = v
                else:
                    off0[0][j - d] = v
                diag = [diag0, [[F(1) if i == jj else F(0) for jj in range(d)]
                                for i in range(d)]]
                off = [off0]
                lin = [F(0)] * (2 * d)
                got = riccati.border_band_psd(diag, off, lin, F(5))
                assert got is False, (d, j, v)
                assert is_psd(_assemble_bordered(diag, off, lin, F(5))) is False
    # border sub-check at a zero pivot strictly INSIDE a d>1 block (k=1),
    # complementing the existing k=0 block case above.
    assert riccati.border_band_psd(
        [[[F(1), F(0)], [F(0), F(0)]]], [], [F(0), F(2)], F(5)) is False
    assert riccati.border_band_psd(
        [[[F(1), F(0)], [F(0), F(0)]]], [], [F(0), F(0)], F(5)) is True


def test_border_corner_is_arrival_cost():
    # H = [[2]], lin = [2]: M = [[2, 1], [1, c]] is PSD iff c >= 1/2
    assert riccati.border_band_psd([[[F(2)]]], [], [F(2)], F(1, 2)) is True
    assert riccati.border_band_psd([[[F(2)]]], [], [F(2)], F(499, 1000)) is False


def _mm(A, B):
    """Exact matrix product (small blocks)."""
    return [[sum((A[i][t] * B[t][j] for t in range(len(B))), F(0))
             for j in range(len(B[0]))] for i in range(len(A))]


def _mt(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def _rand_frac(rng):
    return F(rng.randint(-3, 3), rng.choice((1, 2, 3)))


def _bidiag_psd_chain(rng, N, d, want_singular):
    """``H = L L^T`` with ``L`` block lower-bidiagonal (``L_kk = Dk``, ``L_{k,k-1} =
    Ek``). Then ``H`` is PSD by construction and exactly block-tridiagonal:
    ``H_kk = Dk Dk^T + Ek Ek^T`` and ``H_{k,k+1} = Dk E_{k+1}^T``, with the (k, k+2)
    blocks identically zero. Zeroing a row of some ``Dk`` makes ``L`` -- hence ``H`` --
    genuinely singular, which is the boundary case the bordered sweep must handle via
    the range condition rather than a strict pivot."""
    Dk = [[[_rand_frac(rng) for _ in range(d)] for _ in range(d)] for _ in range(N)]
    Ek = [[[_rand_frac(rng) for _ in range(d)] for _ in range(d)] for _ in range(N)]
    Ek[0] = [[F(0)] * d for _ in range(d)]
    if want_singular:
        Dk[rng.randrange(N)][rng.randrange(d)] = [F(0)] * d
    diag = [_add(_mm(Dk[k], _mt(Dk[k])), _mm(Ek[k], _mt(Ek[k]))) for k in range(N)]
    off = [_mm(Dk[k], _mt(Ek[k + 1])) for k in range(N - 1)]
    return diag, off


def _add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _matvec(A, x):
    return [sum((A[i][j] * x[j] for j in range(len(x))), F(0)) for i in range(len(A))]


def _solve_exact(A, b):
    """Exact Gauss-Jordan solve of ``A x = b`` for a square, possibly SINGULAR ``A``.
    Free variables are set to zero. Returns ``(x, rank)``, or ``(None, rank)`` when the
    system is inconsistent. Written here, deliberately independent of anything in
    `podium.verify`, so the dual value below is an outside computation."""
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]
    piv_col = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [v / pv for v in M[r]]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][j] - f * M[r][j] for j in range(n + 1)]
        piv_col.append(c)
        r += 1
    if any(all(M[i][j] == 0 for j in range(n)) and M[i][n] != 0 for i in range(r, n)):
        return None, r
    x = [F(0)] * n
    for i, c in enumerate(piv_col):
        x[c] = M[i][n]
    return x, r


def test_border_terminal_pivot_is_dual_value():
    """CLAIM-050: the terminal corner pivot of the bordered sweep is ``g(lam) - t``.

    With ``corner = rho - t`` (``rho = r0 - sum_j lam_j r_j``, the paper's constant
    term) and ``lin = u`` the linear border coupling, the Lagrangian dual value is

        g(lam) = rho + min_x (x^T H x + u^T x) = rho + x*^T H x* + u^T x*,   2 H x* = -u

    -- a genuine minimum, hence well defined, exactly when ``H >= 0`` and ``u`` lies in
    ``range(H)``. Both hold here by construction: ``H = L L^T`` is PSD, and ``u = -2 H y``
    puts ``u`` in ``range(H)`` even when ``H`` is singular. ``g`` is recomputed from an
    independent exact solve (`_solve_exact`), and cross-checked against the closed form
    ``g = rho - y^T H y`` that ``u = -2 H y`` implies -- neither uses the sweep."""
    rng = random.Random(20260815)
    n_singular = n_total = 0
    for i in range(320):
        N = rng.randint(1, 4)
        d = rng.randint(1, 3)
        diag, off = _bidiag_psd_chain(rng, N, d, want_singular=(i % 2 == 0))
        n = N * d
        dense = riccati.assemble(diag, off)
        assert is_psd(dense)                      # L L^T: PSD by construction
        y = [_rand_frac(rng) for _ in range(n)]
        u = [-2 * v for v in _matvec(dense, y)]   # u in range(H), singular H included
        rho = _rand_frac(rng)
        t = _rand_frac(rng)
        corner = rho - t                          # main.tex: r0 - sum lam_j r_j - t

        xstar, rank = _solve_exact([[2 * v for v in row] for row in dense],
                                   [-v for v in u])
        assert xstar is not None                  # consistent: u = -2 H y
        singular = rank < n
        n_singular += singular
        n_total += 1
        g = (rho + sum((xstar[a] * dense[a][b] * xstar[b]
                        for a in range(n) for b in range(n)), F(0))
             + sum((u[a] * xstar[a] for a in range(n)), F(0)))
        # closed form implied by u = -2 H y: min value is -y^T H y
        assert g == rho - sum((y[a] * dense[a][b] * y[b]
                               for a in range(n) for b in range(n)), F(0))

        pivot = riccati.border_terminal_pivot(diag, off, u, corner)
        assert pivot is not None                  # H >= 0 and u in range(H)
        assert pivot == g - t                     # the claim, exactly
        # the accessor and the boolean verdict are the same sweep
        assert riccati.border_band_psd(diag, off, u, corner) == (pivot >= 0)
        assert riccati.border_band_psd(diag, off, u, corner) == is_psd(
            _assemble_bordered(diag, off, u, corner))
    assert n_total >= 300
    assert n_singular >= 100                      # the boundary case is really exercised


def test_border_terminal_pivot_refuses_like_the_verdict():
    # lin outside range(H): no dual value, so no pivot -- and the verdict is False
    assert riccati.border_terminal_pivot([[[F(0)]]], [], [F(1)], F(5)) is None
    assert riccati.border_band_psd([[[F(0)]]], [], [F(1)], F(5)) is False
    # H not PSD: refused before the corner
    assert riccati.border_terminal_pivot([[[F(-1)]]], [], [F(0)], F(5)) is None
    # H = [[2]], lin = [2]: arrival cost = corner - 2^2/(4*2) = corner - 1/2
    assert riccati.border_terminal_pivot([[[F(2)]]], [], [F(2)], F(1, 2)) == F(0)
    assert riccati.border_terminal_pivot([[[F(2)]]], [], [F(2)], F(3)) == F(5, 2)


def test_border_validation_and_no_floats():
    with pytest.raises(ValueError, match="border vector"):
        riccati.border_band_psd([[[F(1)]]], [], [F(1), F(1)], F(0))
    with pytest.raises(TypeError, match="Fraction"):
        riccati.border_band_psd([[[F(1)]]], [], [1.0], F(0))
    with pytest.raises(TypeError, match="Fraction"):
        riccati.border_band_psd([[[F(1)]]], [], [F(1)], 0.5)
