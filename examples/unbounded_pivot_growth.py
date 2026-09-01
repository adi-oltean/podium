"""Exact verification cost is bounded neither by instance height nor by order.

The exact checkers in ``podium.verify`` decide positive semidefiniteness of a
banded (and bordered) rational matrix by an exact ``LDL^T`` sweep. A natural
hope is that the cost of that sweep -- the bit-size of the largest pivot it must
carry -- is controlled by something cheap to read off the instance: the size of
the numbers in the matrix, or the order of the matrix. This example exhibits two
families showing that neither controls it, and machine-checks a proof of each.

Both families are built from the same one-dimensional template: stage states
coupled by a path Laplacian, one unit-radius keep-out per stage, and stage goals
recovered from stationarity, so every instance carries a designed optimum with
its multipliers and the bracket closes exactly.

FAMILY 1 -- height stays constant, the pivot grows without bound.

    Entries of the bordered matrix lie in {-1, 0, 1, 2, 3}; the optimum is the
    origin; unit keep-outs are contacted at c_k = -1; multipliers alternate 1
    and 0. The normalized height stays L_M = 2 at EVERY horizon, while the
    largest interior pivot goes from 12 bits at order 9 to 483 bits at order
    257. Proof that this never stops, every link checked below:

      (a) The border is identically zero and the corner is zero: g_k = -lam_k
          gives q0 = 2 lam, so lin_i = q0_i - lam_i(-2 c_i) = 0, and t = f0(0)
          = r0 makes the corner 0. The interior pivots of the bordered matrix
          are therefore exactly the pivots of H.
      (b) H = tridiag(-1, d, -1) with d_1 = 1, d_j = 3 for even j <= N-1,
          d_j = 2 for odd 3 <= j <= N-1. Its pivots are p_j = D_j / D_{j-1}
          with the leading principal minors obeying the continuant recurrence
              D_j = d_j D_{j-1} - D_{j-2},   D_0 = D_1 = 1,
          and gcd(D_j, D_{j-1}) = gcd(D_1, D_0) = 1 all the way down, so
          bit(p_j) = bitlen(D_j) + bitlen(D_{j-1}) EXACTLY.
      (c) d_j >= 2 for 2 <= j <= N-1 gives, by induction from D_2 = 2 > D_1 = 1,
              D_j >= 2 D_{j-1} - D_{j-2} > D_{j-1} > 0,
          so no pivot vanishes and the sweep runs pivoting-free.
      (d) Composing two steps,
              D_{2m+1} = 5 D_{2m-1} - 2 D_{2m-2} >= 3 D_{2m-1},
          using D_{2m-2} <= D_{2m-1} from (c); hence D_{2m+1} >= 3^m, and since
          3^m >= 2^m, bitlen(D_{2m+1}) >= m+1.
      (e) For even N the largest interior pivot carries at least
          bitlen(D_{N-1}) >= (N-2)/2 + 1 = N/2 bits: linear in the order.

          (The two-step transfer matrix [[5,-2],[3,-1]] has determinant 1 and
          trace 4, so its dominant eigenvalue is 2 + sqrt(3) ~ 3.732 and the
          true rate is ~0.95 bits per index; the measured 483 bits at N = 256
          sits on that line. The proof needs only the factor 3.)

    So no function of the height bounds the pivot: the height is 2 at every
    horizon while the pivot exceeds N/2 bits.

FAMILY 2 -- order stays fixed, the pivot grows without bound.

    At the fixed order 9, with the optimum cycling through (0, 1, -1, 2), unit
    keep-outs contacted at c = x-1, and multipliers 1/q with q = 2^m + 1:
    H = P0 - (1/q) I, so the FIRST pivot is exactly

        p_1 = 2 - 1/q = (2q-1)/q,    gcd(2q-1, q) = 1,

    whose bit-size at q = 2^m + 1 is bitlen(2^{m+1}+1) + bitlen(2^m+1) = 2m+3
    EXACTLY. Admissibility holds for every m >= 1, since P0 is a path Laplacian
    plus the identity, so H >= (1 - 1/q) I > 0, and the contact is exact, so the
    bracket closes. Given any bound B, taking m = ceil(B/2) makes the largest
    interior pivot exceed B at order 9. So no function of the order bounds the
    pivot either.

Only the trusted checkers in ``podium.verify`` are imported; the instances,
heights, pivots and minors are all rebuilt here from scratch, so the checkers
are exercised rather than trusted.

Run from the repository root:

    python3 examples/unbounded_pivot_growth.py
"""

from fractions import Fraction as F
from math import gcd

from podium.verify import bracket, riccati

# ---------------------------------------------------------------- measures --

def bitlen(a: int) -> int:
    return abs(a).bit_length()


def bit(x: F) -> int:
    """Per-rational bit-size: numerator bits + denominator bits."""
    return bitlen(x.numerator) + bitlen(x.denominator)


def ht(values):
    """Normalized height of a tuple: clear one common denominator D; the height
    is max(bitlen(D), max_i bitlen(a_i))."""
    values = [F(v) for v in values]
    D = 1
    for v in values:
        D = D * v.denominator // gcd(D, v.denominator)
    h = bitlen(D)
    for v in values:
        a = v * D
        assert a.denominator == 1
        h = max(h, bitlen(a.numerator))
    return h


def dense_ldl_pivots(M):
    """Exact dense symmetric elimination; returns the n+1 pivots of M."""
    n = len(M)
    A = [[F(M[i][j]) for j in range(n)] for i in range(n)]
    piv = []
    for k in range(n):
        p = A[k][k]
        piv.append(p)
        if p == 0:
            continue
        for i in range(k + 1, n):
            if A[i][k] != 0:
                f = A[i][k] / p
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return piv


# ------------------------------------------------- shared instance builder --

def lap_apply(x):
    n = len(x)
    out = []
    for k in range(n):
        v = F(0)
        if k > 0:
            v += x[k] - x[k - 1]
        if k < n - 1:
            v += x[k] - x[k + 1]
        out.append(v)
    return out


def assemble(x, lam, c, rho):
    """d = 1, A = R = W = 1, one unit-radius keep-out per stage, goals
    recovered from stationarity so the designed point is stationary."""
    N = len(x)
    lx = lap_apply(x)
    g = [lx[k] + x[k] - lam[k] * (x[k] - c[k]) for k in range(N)]
    P0 = [[F(0)] * N for _ in range(N)]
    for k in range(N):
        deg = (1 if k > 0 else 0) + (1 if k < N - 1 else 0)
        P0[k][k] = F(deg) + F(1)
        if k < N - 1:
            P0[k][k + 1] = F(-1)
            P0[k + 1][k] = F(-1)
    q0 = [-2 * g[k] for k in range(N)]
    r0 = sum((g[k] * g[k] for k in range(N)), F(0))
    cons = []
    for k in range(N):
        Pk = [[F(0)] * N for _ in range(N)]
        Pk[k][k] = F(1)
        qk = [F(0)] * N
        qk[k] = -2 * c[k]
        cons.append((Pk, qk, c[k] * c[k] - rho[k] * rho[k]))
    return dict(N=N, x=x, lam=lam, c=c, rho=rho, g=g, P0=P0, q0=q0, r0=r0,
                cons=cons)


def quad(P, q, r, x):
    n = len(x)
    tot = F(0)
    for i in range(n):
        for j in range(n):
            tot += x[i] * P[i][j] * x[j]
        tot += q[i] * x[i]
    return tot + r


def h_blocks(inst):
    N = inst["N"]
    H = [[inst["P0"][i][j] for j in range(N)] for i in range(N)]
    for k in range(N):
        H[k][k] -= inst["lam"][k]
    return ([[[H[k][k]]] for k in range(N)],
            [[[H[k][k + 1]]] for k in range(N - 1)])


def admissible(inst):
    """Every clause of the construction, checked directly: the keep-out is
    active at the designed point, the multipliers are nonnegative, H is PSD,
    and stationarity holds stage by stage."""
    N, x, lam, c, rho = inst["N"], inst["x"], inst["lam"], inst["c"], inst["rho"]
    for k in range(N):
        assert quad(*inst["cons"][k], x) == 0, f"stage {k}: keep-out inactive"
        assert rho[k] > 0 and x[k] != c[k] and lam[k] >= 0
    diag, off = h_blocks(inst)
    assert riccati.block_tridiag_psd(diag, off), "H is not PSD"
    lx = lap_apply(x)
    for k in range(N):
        assert 2 * (lx[k] + x[k] - inst["g"][k]) == lam[k] * 2 * (x[k] - c[k]), \
            f"stage {k}: stationarity violated"
    return True


def certify(inst):
    """Assemble the certificate, run the trusted checkers, build M, measure."""
    N, x, lam, cons = inst["N"], inst["x"], inst["lam"], inst["cons"]
    t = quad(inst["P0"], inst["q0"], inst["r0"], x) - sum(
        (lam[k] * quad(*cons[k], x) for k in range(N)), F(0))
    lb = bracket.certify_lower_bound_multi(
        inst["P0"], inst["q0"], inst["r0"], cons, lam, t)
    ub = bracket.certify_upper_bound_multi(
        inst["P0"], inst["q0"], inst["r0"], cons, x)
    closed = bool(lb) and ub is not None and bracket.closes(t, ub)
    diag, off = h_blocks(inst)
    band_h = riccati.block_tridiag_psd(diag, off)
    lin = [inst["q0"][i] - sum((lam[k] * cons[k][1][i] for k in range(N)), F(0))
           for i in range(N)]
    corner = inst["r0"] - sum((lam[k] * cons[k][2] for k in range(N)), F(0)) - t
    band_m = riccati.border_band_psd(diag, off, lin, corner)
    M = [[F(0)] * (N + 1) for _ in range(N + 1)]
    for i in range(N):
        for j in range(N):
            M[i][j] = inst["P0"][i][j] - (lam[i] if i == j else 0)
        M[i][N] = M[N][i] = lin[i] / 2
    M[N][N] = corner
    m_entries = ([M[i][j] for i in range(N) for j in range(N)]
                 + [lin[i] / 2 for i in range(N)] + [corner])
    return dict(t=t, closed=closed, band_h=band_h, band_m=band_m, M=M,
                lin=lin, corner=corner, L_M=ht(m_entries))


# --------------------------------------- family 1: constant height, growing --

def const_height_family(N):
    x = [F(0)] * N
    lam = [F(1) if k % 2 == 0 else F(0) for k in range(N)]
    c = [F(-1)] * N
    rho = [F(1)] * N
    return assemble(x, lam, c, rho)


def diag_pattern(N):
    """The diagonal pattern of H used by step (b) of the proof."""
    d = []
    for k in range(N):                      # stage k <-> 1-indexed j = k+1
        deg = (1 if k > 0 else 0) + (1 if k < N - 1 else 0)
        d.append(F(deg + 1) - (F(1) if k % 2 == 0 else F(0)))
    return d


def minors(dvec):
    """Continuant recurrence D_j = d_j D_{j-1} - D_{j-2}, step (b)."""
    D = [F(1), dvec[0]]
    for j in range(2, len(dvec) + 1):
        D.append(dvec[j - 1] * D[j - 1] - D[j - 2])
    return D                                 # D[0..N]


def family_1():
    print("FAMILY 1: constant height, pivot growing without bound")
    print(f"  {'N':>4} {'n+1':>4} {'L_M':>4} {'max piv bits':>13} "
          f"{'proved >= N/2':>14} {'closed':>7}")
    measured = {}
    for N in (8, 16, 32, 64, 128, 256):
        inst = const_height_family(N)
        admissible(inst)
        r = certify(inst)
        # -- (a): border identically zero, corner zero
        assert all(v == 0 for v in r["lin"]), "border not identically zero"
        assert r["corner"] == 0, "corner nonzero"
        piv = dense_ldl_pivots(r["M"])
        interior, corner_piv = piv[:-1], piv[-1]
        assert corner_piv == 0
        mx = max(bit(p) for p in interior)
        # -- the measured values
        assert r["L_M"] == 2, f"height drifted at N={N}: {r['L_M']}"
        assert r["closed"] and r["band_h"] and r["band_m"]
        entries = {e for row in r["M"] for e in row}
        assert entries <= {F(v) for v in (-1, 0, 1, 2, 3)}, entries
        # -- (b): pivots are exactly D_j / D_{j-1}, consecutive minors coprime
        d = diag_pattern(N)
        assert d[0] == 1 and d[-1] == 2 and all(
            d[j - 1] == (3 if j % 2 == 0 else 2) for j in range(2, N)), \
            "diagonal pattern deviates from the proof's"
        D = minors(d)
        assert all(v.denominator == 1 for v in D)
        Di = [int(v) for v in D]
        for j in range(1, N + 1):
            assert interior[j - 1] == D[j] / D[j - 1], f"pivot {j} != D_j/D_j-1"
            assert gcd(Di[j], Di[j - 1]) == 1, "consecutive minors not coprime"
        # -- (c): strict monotone positivity through j = N-1
        assert all(Di[j] > Di[j - 1] > 0 for j in range(2, N)), "monotonicity"
        # -- (d): two-step identity, factor >= 3, and D_{2m+1} >= 3^m
        for m in range(1, (N - 2) // 2 + 1):        # 2m+1 <= N-1
            assert Di[2 * m + 1] == 5 * Di[2 * m - 1] - 2 * Di[2 * m - 2]
            assert Di[2 * m + 1] >= 3 * Di[2 * m - 1]
            assert Di[2 * m + 1] >= 3 ** m
        # -- (e): the certified linear lower bound on the max pivot
        assert bit(interior[N - 2]) >= bitlen(Di[N - 1]) >= N // 2
        assert mx >= N // 2, "certified bound violated"
        measured[N] = mx
        print(f"  {N:>4} {N + 1:>4} {r['L_M']:>4} {mx:>13} {N // 2:>14} "
              f"{str(r['closed']):>7}")
    assert measured[8] == 12 and measured[256] == 483, \
        "quoted figures (12 bits at order 9, 483 at order 257) drifted"
    # 'at every horizon' includes odd horizons: height 2, admissible, closing,
    # and the same certified bound with N-2 in place of N-1.
    for N in (9, 33, 257):
        inst = const_height_family(N)
        admissible(inst)
        r = certify(inst)
        assert r["L_M"] == 2 and r["closed"] and r["band_h"] and r["band_m"]
        piv = dense_ldl_pivots(r["M"])
        mx = max(bit(p) for p in piv[:-1])
        assert mx >= (N - 1) // 2
    print("  odd horizons 9/33/257: height 2, admissible and closing, bound holds")
    # the proved geometric floor D_{2m+1} >= 3^m, pushed far past the measured
    # window (a pure integer recurrence, no matrix needed):
    d = diag_pattern(4096)
    Di = [int(v) for v in minors(d)]
    for m in (500, 1000, 2000):
        assert Di[2 * m + 1] >= 3 ** m
    assert bitlen(Di[4095]) >= 2047
    print("  beyond-window floor: D_4095 >= 3^2047 -> pivot >= 2047 bits at "
          "order 4097 (proved)")
    return True


# ------------------------------------------ family 2: fixed order, growing --

FIXED_N = 8
X_CYCLE = (0, 1, -1, 2)


def fixed_order_family(m):
    q = (1 << m) + 1
    x = [F(X_CYCLE[k % 4]) for k in range(FIXED_N)]
    lam = [F(1, q)] * FIXED_N
    c = [x[k] - 1 for k in range(FIXED_N)]
    rho = [F(1)] * FIXED_N
    return assemble(x, lam, c, rho), q


def family_2():
    print("\nFAMILY 2: fixed order 9, pivot growing without bound")
    print(f"  {'m':>5} {'L_M':>5} {'max piv bits':>13} {'p_1 bits':>9} "
          f"{'2m+3':>5} {'closed':>7}")
    rows = {}
    for m in (4, 16, 64, 256, 1024):
        inst, q = fixed_order_family(m)
        admissible(inst)
        r = certify(inst)
        piv = dense_ldl_pivots(r["M"])
        interior, corner_piv = piv[:-1], piv[-1]
        assert corner_piv == 0
        mx = max(bit(p) for p in interior)
        # -- the first pivot exactly, and the certified bound
        assert interior[0] == F(2 * q - 1, q), "first pivot != (2q-1)/q"
        assert gcd(2 * q - 1, q) == 1
        assert bit(interior[0]) == 2 * m + 3, "p_1 bit-size != 2m+3"
        assert mx >= 2 * m + 3, "certified fixed-order bound violated"
        # the height is unbounded too: the common denominator of M's entries is
        # a multiple of q, so L_M >= bitlen(q) = m+1.
        assert r["L_M"] >= m + 1
        assert r["closed"] and r["band_h"] and r["band_m"]
        rows[m] = (r["L_M"], mx)
        print(f"  {m:>5} {r['L_M']:>5} {mx:>13} {bit(interior[0]):>9} "
              f"{2 * m + 3:>5} {str(r['closed']):>7}")
    assert rows[4] == (10, 85) and rows[256] == (262, 4116), \
        "quoted figures (height 10->262, pivot 85->4116) drifted"
    print("  for any bound B: m = ceil(B/2) gives pivot > B at order 9 "
          "(proved; m=1024 shown)")
    return True


def run():
    ok = family_1() and family_2()
    print("\nPASS: both growth results are proved, not sampled:" if ok
          else "FAIL")
    print("  - constant height: max interior pivot >= N/2 bits "
          "(D_{2m+1} >= 3^m via the continuant recurrence), height 2 at every "
          "horizon;")
    print("  - fixed order 9: the first pivot (2q-1)/q carries exactly 2m+3 "
          "bits at q = 2^m+1, so the pivot is arbitrarily large at fixed "
          "order.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(run())
