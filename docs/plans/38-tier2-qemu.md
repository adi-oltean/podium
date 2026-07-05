# 38 — Tier-2 cross-architecture golden vectors (qemu-aarch64)

GitHub issue: https://github.com/adi-oltean/podium/issues/38

## The claim

Tier-1 shows the emitted kernels reproduce CPython bit-for-bit on the
host (x86). Tier-2 shows they do so on a DIFFERENT ISA: aarch64 (the
architecture class of most flight processors). IEEE-754 binary64
arithmetic with `-ffp-contract=off` is architecture-independent — no
x87 excess precision, no FMA contraction — so the result must be
identical bit-pattern on x86 SSE2 and ARM NEON. This retires "but does
it behave the same on the target?" with a receipt.

## Fix (landed)

- `tools/tier2_build_run.sh`: inside a debian container, apt-installs
  `gcc-aarch64-linux-gnu` + `qemu-user-static`, cross-compiles the
  emitted kernels + the golden driver STATICALLY for aarch64, and
  replays each kernel's recorded vectors under qemu-aarch64-static.
  Optionally builds the correctly-rounded stm unit against the vendored
  CORE-MATH sources.
- `tests/test_tier2_qemu.py`: REUSES the exact test_cemit harness (same
  `_DRIVER`, `_dispatch_block`, `_vectors`, `_py_call`, tolerance
  sets), writes per-kernel inputs + the CPython reference, runs the
  container, and compares:
  - arithmetic + sqrt kernels: BIT-EXACT vs CPython on aarch64;
  - libm-trig kernels: the same cross-libm class as x86 (bounded by
    output scale, <=1% incidence — aarch64 glibc vs CPython conda);
  - stm in CORE-MATH mode: BIT-EXACT vs the mpmath correctly-rounded
    oracle on aarch64 — the trig tolerance is gone cross-arch, not
    merely bounded.
- `.github/workflows/tier2.yml`: emit paths + weekly + dispatch.

## Result

Probe and full run confirm bit-identical `0.1+0.2`, `sqrt(2)`, and
every arithmetic/sqrt kernel across ISAs; stm-CR bit-exact on ARM.

## Deferred

Wiring CR mode into the DEFAULT tier-2 comparison for all trig kernels
(stm is the demonstrator today); big-endian and 32-bit targets; a real
RTOS/flight-board run.

## Push/merge instructions

Single commit on main: `38 — Tier-2 cross-architecture golden vectors
(#38)`; push; close.
