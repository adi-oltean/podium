# 40 — Verified-KKT checker (exact-rational optimality certificate)

GitHub issue: https://github.com/adi-oltean/podium/issues/40

## Fix (landed) — `podium.verify.kkt`

The R4 "certificate-carrying optimization" pattern, checker half. An
online convex solver (Clarabel, ECOS, an embedded QOCOGEN kernel) is
UNTRUSTED; `verify_qp(P,q,G,h,A,b, x,mu,nu)` re-verifies its claimed
solution in `fractions.Fraction` arithmetic and returns a `KKTReport`
of EXACT residuals:

- stationarity  max|P x + q + G' mu + A' nu|
- eq_residual   max|A x - b|
- ineq_violation max(0, (G x - h)_i)   (primal infeasibility)
- dual_violation max(0, -mu_i)         (dual infeasibility)
- duality_gap   sum_i mu_i (h - G x)_i (bounds p(x) - p* when feasible)

For a CONVEX QP (P symmetric PSD) these conditions are necessary and
sufficient for global optimality, so `report.certified(tol)` is a
rigorous near-optimality certificate — and the arithmetic being exact
means the residuals ARE the KKT violations, with no float uncertainty
in the check itself. `rationalize_vec/_mat` convert solver floats to
exact rationals before checking (trusted-checker discipline, as in
verify.barrier: no floats in the verification path).

## Receipts (tests/test_kkt.py)

- Equality QP min 1/2||x||^2 s.t. x1+x2=2 -> x=(1,1), nu=-1: every
  residual EXACTLY Fraction(0), certified at tol=0.
- Inequality QP min 1/2 x^2 s.t. x>=1 -> x=1, mu=1 active: comp
  slackness and stationarity EXACTLY zero.
- Untrusted Clarabel solution of a Layer-0 min-energy rendezvous QP
  (CW-STM reachability, 3 impulses) certifies within 1e-5, and the
  certified objective matches the solver's.
- Perturbed primal: exact eq_residual = 1/10, rejected.
- Flipped dual sign: dual_violation = 1, rejected.
- Non-symmetric P: flagged structurally.

## Why exactness matters

A float KKT check can only say "residual < eps in float"; rounding
could hide a real violation or invent one. Rationalizing first and
checking in Fraction makes the residual a mathematical fact about the
exact point — the same discipline as the barrier certificates (#20).

## Deferred

Feeding an embedded CVXPYgen/QOCOGEN solver's KKT dump directly to the
checker; SOCP cone complementarity
(second-order-cone membership in exact arithmetic) for the conic
Layer-0 problems; an interval-arithmetic variant for problems without
a clean rational optimum.

## Push/merge instructions

Single commit on main: `40 — Verified-KKT checker (#40)`; push; close.
