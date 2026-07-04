# 30 — Emitter v2: matmul lowering, EKF to C

GitHub issue: https://github.com/adi-oltean/podium/issues/30

## Fix (landed)

Emitter v2 (`podium.emit.cemit`):
- `@shapes` decorator (podium.verify.contracts, passive metadata):
  pure matrix-algebra kernels never reveal shapes via subscripts.
- 2-D array parameters (`const double p[6][6]`), tuple returns
  lowered to multiple out-parameters (out0, out1), multiple local
  array allocations, AnnAssign normalization, loop-var +/- constant
  subscripts (q[i, i+3]).
- Array-expression lowering: static shape inference (`_ashape`), then
  matmul chains / transposes / elementwise add-sub / scalar scaling
  emitted as fixed-bound loops with temporaries; transposes handled by
  index swap at use sites. `ekf.predict` emits AS WRITTEN
  (phi @ p @ phi.T + q, symmetrized).

New flight kernel `ekf.update_sequential`: sequential scalar Joseph
updates (H = [I3 0], R = r_var*I3) — one component at a time, scalar
innovation variance, DIVISION instead of np.linalg.solve; receipt
proves equivalence to batch Joseph for diagonal R (200 random PSD
cases, 1e-9 relative). A covariance-repair clamp (s = max(s, r_var),
written as if/assign because max() is outside the C subset) makes the
division provably safe.

## Golden-vector classes (test_cemit)

- matmul kernels (`predict`): NumPy @ uses BLAS accumulation order;
  naive row-major loops CANNOT be bit-exact by construction. Bounded
  at 1e-13 x output-vector scale (6-term dot reassociation).
- `update_sequential`, `process_noise_wna`: explicit scalar loops in
  both languages -> strict BIT-exact class. 20 kernels total.

## EVA gate findings (the gate earned its keep)

1. Flat `(double *)` cast fills in the driver are not recognized as
   initializing a 2-D object -> 4 spurious uninitialized-value alarms;
   driver now writes per-element nested loops.
2. REAL envelope finding: p in [-1e6,1e6] with r_var >= 1e-6 admits
   K ~ p/r_var ~ 1e12, and the interval semantics overflow to inf by
   the third sequential measurement. The stated envelope is now the
   physical one — covariance <= 1e4 ((100 m std)^2), r_var >= 1e-2
   ((10 cm std)^2) — documented at the range declaration and in the
   r_var contract. Re-proven: 0 alarms, 41 functions, 643 statements,
   138/138 preconditions valid.

## Deferred

Function-typed parameters (integrators), CompCert-subset audit,
moving DEFAULT_RANGES into @contracts.

## Push/merge instructions

Single commit on main: `30 — Emitter v2: matmul lowering, EKF to C
(#30)`; push; close.
