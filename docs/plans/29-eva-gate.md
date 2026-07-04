# 29 — Sound value gate: Frama-C/EVA over the emitted C

GitHub issue: https://github.com/adi-oltean/podium/issues/29

## Unblocking story

`sudo apt install frama-c` cannot work: frama-c is absent from the
Ubuntu 22.04 archive entirely (dropped in an OCaml transition; the
universe index was verified complete at 58,923 packages with zero
frama* entries). The gate runs the official framac/frama-c:dev image
(Frama-C 34 Selenium) instead; on this WSL host the Docker Desktop
credential helper breaks anonymous pulls, worked around with a clean
DOCKER_CONFIG.

## Fix (landed)

- `podium.emit.evagen`: generates the EVA driver — one check function
  per kernel constructing every input with Frama_C_double_interval
  from the @contract interval where declared, else from DEFAULT_RANGES
  (each such entry is a CONTRACT GAP, and the generator writes the gap
  list into the driver header so the assumption is visible in the
  artifact).
- `tools/eva_gate.py`: emits kernels + driver, runs EVA in the
  container (-eva -main eva_main -eva-precision 2), saves the full
  report, parses the alarm count; exit 0 only at zero alarms.
- `.github/workflows/eva.yml`: re-proves on core/emit changes +
  weekly; uploads eva_report.txt as the audit artifact.
- Emitter fix found by the gate: ACSL decimal literals are exact
  REALS, so `1e-05 <= n` is unprovable for the nearest double (it
  sits below the real) — 4 preconditions stuck at 'unknown' until
  the ACSL renderer switched to hex float literals.

## Result (the headline)

0 alarms. 35 functions, 364 statements, 100% coverage. 115/115 ACSL
preconditions valid. 100% of reached logical properties proven. That
is a SOUND float-interval proof — no division by zero, no overflow to
infinity, no invalid access — for every input in the contracted/stated
ranges, over all 17 emitted kernels including the 20-iteration Newton
solve (EVA bounds 1 - e*cos(E) >= 0.1 from the e-contract and proves
the division safe).

## Deferred

Moving DEFAULT_RANGES entries into kernel @contracts (tracked gap
list); WP/RTE for functional postconditions; memory gate extension as
kernels gain pointers (none today).

## Push/merge instructions

Single commit on main: `29 — EVA sound value gate (#29)`; push;
close; dispatch eva.yml once to confirm green in CI.
