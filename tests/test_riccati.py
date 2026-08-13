"""Exact block-tridiagonal (Riccati) PSD certificate -- podium.verify.riccati.

The trusted claim: `block_tridiag_psd` gives the SAME verdict as the dense
`barrier.is_psd` on the assembled matrix, but exploits the band (O(N d^3)); and
`riccati_storage` returns the cost-to-go / storage-function pivot blocks for the
benign positive-definite regime. These tests pin soundness (equivalence to is_psd
over random instances), every edge branch, and the no-float discipline.
"""

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


def test_border_corner_is_arrival_cost():
    # H = [[2]], lin = [2]: M = [[2, 1], [1, c]] is PSD iff c >= 1/2
    assert riccati.border_band_psd([[[F(2)]]], [], [F(2)], F(1, 2)) is True
    assert riccati.border_band_psd([[[F(2)]]], [], [F(2)], F(499, 1000)) is False


def test_border_validation_and_no_floats():
    with pytest.raises(ValueError, match="border vector"):
        riccati.border_band_psd([[[F(1)]]], [], [F(1), F(1)], F(0))
    with pytest.raises(TypeError, match="Fraction"):
        riccati.border_band_psd([[[F(1)]]], [], [1.0], F(0))
    with pytest.raises(TypeError, match="Fraction"):
        riccati.border_band_psd([[[F(1)]]], [], [F(1)], 0.5)
