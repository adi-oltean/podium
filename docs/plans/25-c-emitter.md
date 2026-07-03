# 25 — C emitter v0 + tier-1 golden vectors

GitHub issue: https://github.com/adi-oltean/podium/issues/25

## Fix (landed) — `podium.emit.cemit`

AST-based C99 emitter; the supported-subset checker REJECTS anything
outside (EmitError), which makes the emitter the operational definition
of StaticPy: pure functions, float scalars + fixed arrays,
constant-index subscripts (1-D/2-D), Python-association-preserving
arithmetic (fully parenthesized), whitelisted math.*, np.empty/zeros
only as return-array allocation, if/else + conditional expressions,
cross-kernel calls lowered through explicit temporaries (array-ness
propagates through call sites in a fixpoint pass). Contracts render as
ACSL requires clauses (with \forall for array params) plus the
analyzer [spec] block. Array returns become out-parameters.

## Receipts (all green; 2000 seeded vectors per kernel incl. branches)

- BIT-exact Python<->C for every arithmetic+sqrt kernel (quaternion
  family, mean_motion, cw_deriv) under gcc -std=c99 -O2
  -ffp-contract=off — sqrt is IEEE correctly-rounded, so same inputs
  give same bits across toolchains;
- the measured exception, asserted rather than hidden: sin/cos in the
  interpreter's (conda) libm vs system glibc differ on 21/72,000 stm
  values, <=4 ulp after propagation, <0.1% incidence — precisely the
  gap the roadmap's CORE-MATH item exists to close;
- subset rejection (loops, lists) raises EmitError;
- ACSL + [spec] rendering pinned textually;
- emission is byte-deterministic.

Debug war stories recorded: NumPy 2's repr broke the first vector
transport (np.float64(x) fed to strtod — switched to hex floats,
exact both ways); parameters used only via cross-kernel calls needed
array-ness propagation.

## Deferred (v0.6 items, authored in the roadmap)

Bounded-for-loop emission for the remaining core, CompCert-subset
audit, EVA gate over the emitted C, tier-2 on-target vectors,
CORE-MATH transcendentals.

## Push/merge instructions

Single commit on main: `25 — C emitter v0 + tier-1 golden vectors
(#25)`; push; close.
