# 10 — ARCH rendezvous reachability regression in CI

GitHub issue: https://github.com/adi-oltean/podium/issues/10

## Problem

Closed-loop safety is currently argued by simulation only. The ARCH-COMP
spacecraft rendezvous benchmark (Chan & Mitra ARCH17; planar CW, switched
approach/attempt/abort controller, 30-deg LOS cone, 0.055 m/s velocity
octagon, passive-abort target avoidance) is the community-standard hybrid
verification case — solved annually and manually by tool authors, never
(per the July 2026 sweep) wired as a third-party CI regression. Podium
ships it as: executable model + machine-readable export + JuliaReach
verification gated in CI.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/guidance/arch.py` | MISS | model, deterministic hybrid sim, spec margins, JSON export |
| `tools/reach/arch_rendezvous.jl` + `Project.toml` | MISS | consumes the export, proves SRNA01+SRA01 |
| `.github/workflows/reach.yml` | MISS | path-filtered + manual + weekly |
| `tests/test_arch.py` | MISS | python-side receipts |
| roadmap | PARTIAL | check off the benchmark item |

## Fix

- Model data with provenance (matrices from the ARCH repeatability
  packages; units: meters, minutes). Receipt tying it to Podium's own
  kernel: the abort-mode matrix IS planar CW at GEO mean motion
  (2n = 0.00876276 /min) — asserted against `cw.cw_deriv` structure.
- `export_model(abort_time)` -> JSON hybrid automaton (modes: A, b,
  invariants; transitions: guards; initial set; properties). The Julia
  side consumes this file — Podium exports, the tool verifies.
- `simulate()` deterministic hybrid RK4 with urgent guard semantics;
  `spec_margins()` for LOS/velocity/target-avoidance per mode.
- Julia script mirrors the ARCH-COMP 2025 JuliaReach settings (BOX
  delta=0.04, LazyClustering, box-template intersection), checks the
  three properties, exits nonzero on violation. Run locally once before
  shipping (receipt), then in CI on guidance/control/sim/core changes.

## Acceptance Criteria

- [ ] Python receipts green (structure, stability, export schema,
      simulation spec margins from initial-set corners)
- [ ] Local JuliaReach run verifies SRNA01 and SRA01 from the export
- [ ] CI workflow in place (path-filtered + manual + weekly)
- [ ] Roadmap updated

## Push/merge instructions

Single commit on main: `10 — ARCH reachability regression (#10)`; push;
close; verify the reach workflow runs green via workflow_dispatch.

## Verification steps

`gh workflow run reach.yml` and confirm the verification step passes;
locally: `julia --project=tools/reach tools/reach/arch_rendezvous.jl
tmp/arch_model.json`.
