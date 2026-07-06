# Exact optimality-gap certificates for nonconvex QCQPs

This document states the four results implemented in
[`src/podium/verify/bracket.py`](../src/podium/verify/bracket.py) and maps each
to the function that realizes it and the test that checks it. The numbering here
is the one used in the accompanying paper. Every bound below is an exact rational
(`fractions.Fraction`); the trusted checker uses no floating point.

## Problem

For a quadratically-constrained quadratic program (QCQP)

```
J* = inf_x  f0(x)   subject to   f_k(x) >= 0,  k = 1..m,
     f_j(x) = x' P_j x + q_j' x + r_j   (P_j symmetric),
```

a keep-out constraint `||x - c|| >= R` makes the feasible set nonconvex, so a
solver can return a feasible point but not, on its own, a proof that no better
point exists. We bracket `J*` between two exact-rational certificates.

For `lambda >= 0` (componentwise) and `t` define the `(n+1) x (n+1)` symmetric
matrix

```
M(lambda, t) = [ P0 - sum_k lambda_k P_k      (q0 - sum_k lambda_k q_k)/2 ]
               [ (q0 - sum_k lambda_k q_k)'/2  r0 - sum_k lambda_k r_k - t ]
```

so that, with `z = [x; 1]`, `z' M(lambda, t) z = f0(x) - sum_k lambda_k f_k(x) - t`.

## Theorem 1 (lower bound; soundness)

If `lambda >= 0` and `M(lambda, t) >= 0` (positive semidefinite), then
`t <= J*`. Symmetrically, any exactly feasible `x-bar` gives `J* <= f0(x-bar)`.

*Idea.* `M >= 0` makes `f0(x) - sum lambda_k f_k(x) - t >= 0` for all `x`; at a
feasible point each `f_k >= 0` and `lambda_k >= 0`, so `f0(x) >= t`. This is
Lagrangian weak duality and needs no convexity. The upper leg checks feasibility
of `x-bar` exactly, since a point feasible only within a solver tolerance could
have objective below `J*` and place the upper bound below the optimum.

- **Code:** `certify_lower_bound`, `certify_upper_bound`, `lower_bound_matrix`
  (single constraint); `certify_lower_bound_multi`, `certify_upper_bound_multi`,
  `lower_bound_matrix_multi` (several constraints).
- **Tests:** `test_lower_bound_certificate_is_exact_and_sound`,
  `test_upper_bound_rejects_tolerance_feasible_point`.

### Corollary (exact optimum)

If a certified lower bound `t` equals a certified upper bound `u = f0(x-bar)`,
then `J* = t = u`: an exact certificate of the global optimum of a nonconvex
program. `certified_optimum` binds both bounds to a single problem's data so
that bounds from different instances are never compared.

- **Code:** `certified_optimum`, `closes`.
- **Tests:** `test_bracket_closes_to_exact_global_optimum`,
  `test_bracket_is_a_real_exact_gap_for_a_suboptimal_cut`,
  `test_certified_optimum_binds_provenance_across_problems`.

## Theorem 2 (nonsingular recovery)

For one constraint, let `A(lambda) = P0 - lambda P1` and
`g(lambda) = inf_x (f0 - lambda f1)` on `U = { lambda : A(lambda) > 0 }`. If the
optimal multiplier `lambda*` lies in `U` (so `A(lambda*)` is nonsingular), then
`g` is smooth and concave near `lambda*` with `g'(lambda*) = 0`, so rounding an
approximate dual `lambda ~ lambda*` to denominator at most `D` gives an exactly
certified lower bound with `J* - g(lambda) = O(1/D^2)` (quadratic). If `lambda*`
and `J*` are rational, the bracket closes exactly.

*Idea.* `g(lambda) = c - (1/4) b' A(lambda)^{-1} b` by completing the square, a
rational function of `lambda`, computed exactly with a rational linear solve. It
is the largest `t` with `M(lambda, t) >= 0` (Schur complement).

- **Code:** `dual_value` (exact `g(lambda)` via `_solve`), `recover_lower_bound`
  (rounds a float dual, then re-certifies).
- **Tests:** `test_dual_value_gives_certified_lower_bounds`,
  `test_recovery_converges_to_optimum_as_lambda_approaches_lam_star`,
  `test_recover_exact_optimum_from_float_dual`.

## Theorem 3 (singular hard case)

If `A(lambda*)` is positive semidefinite but singular (the trust-region hard
case): (a) if `lambda*` and `J*` are rational, `(lambda*, J*)` is still an exact
certificate, since `M` is a rank-deficient PSD matrix that the exact test
accepts, so soundness is unaffected; (b) approaching `lambda*` from within `U`
converges to `J*` only linearly; (c) exact closure requires `lambda*, J*`
rational, which is not generic (`J*` may be irrational).

Soundness is retained; only the recovery rate falls from quadratic to linear.
`dual_value` returns `None` at a singular `lambda*` (it requires `A > 0`), yet
`certify_lower_bound(lambda*, J*)` still verifies.

- **Test:** `test_hard_case_singular_optimum`.

## Theorem 4 (several constraints; certified duality gap)

Theorem 1 and its corollary hold for any `m`. For `m >= 2` the S-lemma need not
apply, so the largest certifiable lower bound (the Shor relaxation value) may be
strictly below `J*`. The bracket then certifies a duality gap `u - t > 0` whose
endpoints are both exact rationals. Higher levels of the moment hierarchy give
nondecreasing certified lower bounds converging to `J*` under a compactness
condition, each verified by an exact PSD test.

- **Code:** `certify_lower_bound_multi`, `certify_upper_bound_multi`.
- **Test:** `test_multi_constraint_bracket`.

## Notes

- The exact PSD test is a symmetric (`LDL^T`) Gaussian elimination in exact
  rationals, in `O(n^3)`: a symmetric matrix is PSD iff every pivot is `>= 0` and
  each zero pivot leaves the rest of its row (hence, by symmetry, its column) in
  the Schur complement zero. The factorization `M = L D L^T` with `D >= 0` is
  itself the nonnegativity witness. This avoids eigenvalues (in general
  irrational) and the leading-minor Sylvester test (which certifies strict
  definiteness only), and replaces an earlier all-principal-minors test that was
  exponential in `n`.
- Full proofs are given in the accompanying paper; this document states the
  results and their code/test realization.
