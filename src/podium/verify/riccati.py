"""Exact block-tridiagonal (Riccati) PSD certificates for trajectory QCQPs.

For a block-banded trajectory QCQP the S-procedure Lagrangian Hessian
``H = P0 - sum_k lam_k P_k`` is **block-tridiagonal** (each stage couples only its
neighbours; each per-stage keep-out touches only its own block). Its exact PSD test --
the trusted acceptance step of the optimality-gap certificate (`podium.verify.bracket`) --
then needs only an ``LDL^T`` that eliminates **within the band**: ``O(n * w^2) = O(N d^3)``
exact rational operations instead of the dense ``O(n^3)`` of `barrier.is_psd`.

Restricted to the band this ``LDL^T`` *is* the discrete **Riccati** / block-Schur recursion
of an LQ problem (classical, Wright 1993; Rao-Wright-Rawlings 1998): eliminating blocks in
increasing order, the pivot blocks are the accumulated arrival-cost / information-form Hessians

    S_0 = H_00,   S_k = H_kk - H_{k,k-1} S_{k-1}^{-1} H_{k-1,k}.

`block_tridiag_psd` is the sound verdict (provably equal to `barrier.is_psd` on the
assembled matrix -- see tests) and exploits the band; `riccati_storage` returns the
storage-function pivot blocks for the benign positive-definite regime.

Same "check the answer, not the run", no-float-in-the-trusted-path discipline as
barrier / bracket / kkt: every entry must be an exact `Fraction`.
"""

from __future__ import annotations

from fractions import Fraction as F

from podium.verify.barrier import is_psd

Vec = list[F]
Mat = list[list[F]]


def _require_fractions(blocks: list[Mat], what: str) -> None:
    """Reject any non-Fraction entry: floats/NaN break exact-arithmetic soundness."""
    for b in blocks:
        for row in b:
            for v in row:
                if not isinstance(v, F):
                    raise TypeError(
                        f"{what} requires exact Fraction entries; "
                        f"got {type(v).__name__}")


def _validate(diag: list[Mat], off: list[Mat]) -> int:
    """Shape checks: N square d x d diagonal blocks, N-1 d x d off blocks (the
    (k, k+1) super-diagonal). Returns d. Raises ValueError on malformed input."""
    if not diag:
        raise ValueError("need at least one diagonal block")
    d = len(diag[0])
    if d == 0:
        raise ValueError("block dimension must be positive")
    if len(off) != len(diag) - 1:
        raise ValueError("off-diagonal block count must be len(diag) - 1")
    for b in diag:
        if len(b) != d or any(len(row) != d for row in b):
            raise ValueError("diagonal blocks must be square and uniform")
    for b in off:
        if len(b) != d or any(len(row) != d for row in b):
            raise ValueError("off-diagonal blocks must match the block size")
    _require_fractions(diag, "block-tridiagonal PSD certificate")
    _require_fractions(off, "block-tridiagonal PSD certificate")
    return d


def assemble(diag: list[Mat], off: list[Mat]) -> Mat:
    """Dense symmetric block-tridiagonal matrix from diagonal blocks ``diag[k]`` and
    super-diagonal blocks ``off[k]`` (the (k, k+1) block; (k+1, k) is its transpose).
    Only band entries are filled, so the result has scalar bandwidth ``2d - 1``."""
    d = _validate(diag, off)
    N = len(diag)
    n = N * d
    m: Mat = [[F(0)] * n for _ in range(n)]
    for k in range(N):
        for i in range(d):
            for j in range(d):
                m[k * d + i][k * d + j] = diag[k][i][j]
    for k in range(N - 1):
        for i in range(d):
            for j in range(d):
                m[k * d + i][(k + 1) * d + j] = off[k][i][j]
                m[(k + 1) * d + j][k * d + i] = off[k][i][j]
    return m


def _band_ldlt_psd(a: Mat, w: int) -> bool:
    """Exact PSD test by ``LDL^T`` eliminating only within scalar bandwidth ``w``.
    For a symmetric matrix whose nonzeros satisfy ``|i-j| <= w`` the factor keeps the
    same bandwidth, so this returns the same verdict as the dense `barrier.is_psd` in
    ``O(n w^2)``. Mirrors `is_psd`'s exact zero-pivot handling."""
    n = len(a)
    for k in range(n):
        if a[k][k] < 0:
            return False
        hi = min(n, k + w + 1)
        if a[k][k] == 0:
            if any(a[k][j] != 0 for j in range(k + 1, hi)):
                return False
            continue
        piv = a[k][k]
        for i in range(k + 1, hi):
            if a[i][k] != 0:
                f = a[i][k] / piv
                for j in range(k, hi):
                    a[i][j] -= f * a[k][j]
    return True


def block_tridiag_psd(diag: list[Mat], off: list[Mat]) -> bool:
    """Exact PSD certificate for the block-tridiagonal ``H`` given by its blocks, via
    the band-restricted (Riccati) ``LDL^T`` in ``O(N d^3)``. Sound: returns the same
    verdict as `barrier.is_psd(assemble(diag, off))`. Refuses non-Fraction entries."""
    d = _validate(diag, off)
    m = assemble(diag, off)
    return _band_ldlt_psd([row[:] for row in m], 2 * d - 1)


def border_band_psd(diag: list[Mat], off: list[Mat], lin: Vec, corner: F) -> bool:
    """Exact PSD certificate for the bordered matrix

        M = [[ H,        lin/2 ],
             [ lin^T/2,  corner ]]

    with ``H`` the block-tridiagonal matrix given by its blocks -- the full
    S-procedure LMI ``M(lam, t)`` of `podium.verify.bracket` when ``H`` is the
    Lagrangian Hessian, ``lin = q0 - sum_k lam_k q_k`` and ``corner =
    r0 - sum_k lam_k r_k - t``. One band ``LDL^T`` sweep carries the border
    column: each pivot updates its band window plus its border entry, a zero
    pivot demands a zero row *including the border entry* (this is exactly the
    ``lin in range(H)`` condition of the generalized Schur complement, detected
    as a consistent zero-pivot row), and the terminal corner pivot is the
    arrival cost ``corner - lin^T H^+ lin / 4``. Same verdict as the dense
    `barrier.is_psd` on the assembled ``M``, in ``O(N d^3)`` exact rational
    operations. Refuses non-Fraction entries."""
    d = _validate(diag, off)
    n = len(diag) * d
    if len(lin) != n:
        raise ValueError("border vector length must match the matrix dimension")
    if any(not isinstance(v, F) for v in lin) or not isinstance(corner, F):
        raise TypeError(
            "bordered PSD certificate requires exact Fraction entries")
    w = 2 * d - 1
    a = assemble(diag, off)
    col = [v / F(2) for v in lin]
    c = corner
    for k in range(n):
        if a[k][k] < 0:
            return False
        hi = min(n, k + w + 1)
        if a[k][k] == 0:
            if any(a[k][j] != 0 for j in range(k + 1, hi)) or col[k] != 0:
                return False
            continue
        piv = a[k][k]
        for i in range(k + 1, hi):
            if a[i][k] != 0:
                f = a[i][k] / piv
                for j in range(k, hi):
                    a[i][j] -= f * a[k][j]
                col[i] -= f * col[k]
        c -= col[k] * col[k] / piv
    return c >= 0


def _inv(a: Mat) -> Mat:
    """Exact inverse via Gauss-Jordan (small d). Raises ValueError if singular."""
    n = len(a)
    aug = [[a[i][j] for j in range(n)] + [F(i == j) for j in range(n)]
           for i in range(n)]
    for c in range(n):
        piv = next((r for r in range(c, n) if aug[r][c] != 0), None)
        if piv is None:
            raise ValueError("singular pivot block (not positive definite)")
        aug[c], aug[piv] = aug[piv], aug[c]
        p = aug[c][c]
        aug[c] = [v / p for v in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [aug[r][j] - f * aug[c][j] for j in range(2 * n)]
    return [row[n:] for row in aug]


def riccati_storage(diag: list[Mat], off: list[Mat]) -> list[Mat]:
    """The Riccati / block-Schur pivot blocks ``S_k`` (the accumulated arrival-cost /
    information-form Hessians) for the benign **positive-definite** regime:

        S_0 = diag[0],  S_k = diag[k] - off[k-1]^T S_{k-1}^{-1} off[k-1].

    Returns ``[S_0, ..., S_{N-1}]``, each an exact ``d x d`` block. Raises ValueError if
    a pivot is singular (use `block_tridiag_psd` for the general semidefinite verdict)."""
    d = _validate(diag, off)
    pivots: list[Mat] = []
    prev: Mat | None = None
    for k in range(len(diag)):
        if prev is None:
            s = [row[:] for row in diag[0]]
        else:
            b = off[k - 1]
            bt = [[b[j][i] for j in range(d)] for i in range(d)]
            inv_prev = _inv(prev)
            # S_k = diag[k] - b^T prev^{-1} b
            tmp = [[sum((inv_prev[i][t] * b[t][j] for t in range(d)), F(0))
                    for j in range(d)] for i in range(d)]
            corr = [[sum((bt[i][t] * tmp[t][j] for t in range(d)), F(0))
                     for j in range(d)] for i in range(d)]
            s = [[diag[k][i][j] - corr[i][j] for j in range(d)] for i in range(d)]
        if not is_psd(s):
            raise ValueError(f"pivot block {k} is not PSD; chain is not benign")
        pivots.append(s)
        prev = s
    return pivots
