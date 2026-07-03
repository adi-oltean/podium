# 27 — Emitter v1: bounded loops, YA Kepler + ROE kernels

GitHub issue: https://github.com/adi-oltean/podium/issues/27

## Fix (landed)

cemit v1 additions: `for i in range(N)` with N a literal or module-level
integer constant (resolved via the function's globals — compile-time
loop bounds, per the subset), module-level float-constant inlining
(_TWO_PI), tuple assignment (sequenced only when provably safe: no
target read in the RHS), augmented assignment, loop-variable subscripts,
np.eye return allocation with identity init.

New emitted kernels (17 total now): ya.kepler_eccentric (bounded
20-step Newton — the flagship loop), true/eccentric anomaly
conversions, propagate_true_anomaly, roe.stm_keplerian +
map_roe_to_lvlh + map_lvlh_to_roe + control_matrix — the ROE kernels
are the first CONTRACTED functions through the ACSL rendering path
end-to-end.

## Tier-1 calibration (measured, not assumed)

- arithmetic+sqrt kernels: strict BIT-exact, unchanged (9 kernels).
- trig-bearing kernels across two libms (conda CPython vs system
  glibc): divergence enters at ~1 ulp of each trig RESULT and
  propagates through arithmetic whose intermediates dwarf cancelling
  outputs — measured: an stm entry at 10 output-ulps, a
  map_roe_to_lvlh velocity diverging by 1 ulp of its ~14-magnitude
  intermediate while the output is 0.031. Per-value bound is therefore
  1e-12 x output-vector scale (translation bugs are O(value) and still
  fail loudly), incidence <= 1% (measured 0.03-0.15%). The Newton
  receipt: convergence contracts the divergence — kepler_eccentric
  passes at far tighter incidence than its per-iteration trig calls
  would suggest.

## Rejection set updated

while loops and data-dependent range bounds now carry the rejection
receipts (bounded for-range became legal).

## Deferred

Matmul lowering for the EKF kernels, function-typed parameters
(integrators), CompCert-subset audit.

## Push/merge instructions

Single commit on main: `27 — Emitter v1 (#27)`; push; close.
