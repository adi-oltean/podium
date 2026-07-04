# 37 — CORE-MATH correctly-rounded transcendentals

GitHub issue: https://github.com/adi-oltean/podium/issues/37
First delivery of v0.7.

## What it retires

The emitted kernels are bit-exact vs CPython for all arithmetic and
sqrt (IEEE-exact), but sin/cos are not correctly rounded in either the
conda libm or glibc, and the two disagree — the ONE measured tier-1
tolerance (21/72000 stm values at 1 ULP, `_TRANSCENDENTAL` class in
test_cemit). Correct rounding removes the ambiguity: every correctly-
rounded implementation returns the unique nearest double, so they all
agree bit-for-bit.

## Fix (landed)

- `third_party/core-math/`: vendored CORE-MATH sin.c, cos.c (MIT,
  unmodified) exporting cr_sin/cr_cos, plus coremath.h, LICENSE,
  README.
- `emit_module(correctly_rounded=True)`: emits cr_sin/cr_cos in place
  of libm sin/cos and includes "coremath.h". Everything else (double
  arithmetic, sqrt, layout) is unchanged, so it stays inside the
  CompCert subset and the EVA envelope.
- `podium.emit.croracle`: an mpmath correctly-rounded oracle
  (evaluate at 120-bit precision on the exact input double, round once
  to double). mpmath added to the dev extra; imported lazily so the
  core stays dependency-free.

## Receipts (tests/test_coremath.py, both slow)

1. Compiled cr_sin/cr_cos == the mpmath oracle BIT-EXACT over 8000
   arguments spanning the kernels' ranges (unit band, anomaly band,
   n*t up to n*20000, wide band) — proves both sides are correctly
   rounded and agree.
2. stm emitted in CR mode (linked against the vendored sources) ==
   the CR-oracle stm reference BIT-EXACT, 0/3000 — no tolerance, no
   incidence bound. The 21/72000 exception is gone.

The default (libm) golden vectors and their honest tolerance classes
are unchanged; CR mode is the opt-in that removes the last one.

## Deferred

Wiring CR mode into the standard 20-kernel golden suite and the
release/EVA/CompCert lanes (they run default mode today); atan2/other
CORE-MATH functions if future kernels need them; tier-2 qemu-aarch64
using CR mode for a bit-exact cross-target claim.

## Push/merge instructions

Single commit on main: `37 — CORE-MATH correctly-rounded
transcendentals (#37)`; push; close. Opens v0.7.
