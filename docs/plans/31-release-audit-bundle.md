# 31 — Per-tag release: audit bundle + fresh EVA proof

GitHub issue: https://github.com/adi-oltean/podium/issues/31

## Fix (landed)

- `podium.emit.kernels.FLIGHT_KERNELS`: the canonical emitted-kernel
  list (was duplicated across the golden tests and the EVA gate; now
  one source of truth for tests, gate, and release builder).
- `tools/build_audit_bundle.py`: builds <outdir>/{bundle.json,
  kernels.c, eva_driver.c, meta.json}. bundle.json is the
  byte-deterministic reference-mission audit (seed 7); meta.json
  carries tag/commit/toolchain stamps SEPARATELY so bundle.json stays
  byte-comparable across rebuilds of identical code. HARD GATE: exits
  nonzero unless the mission captures, every IDSS and STL margin is
  positive, and the barrier certificate verifies exactly.
- `.github/workflows/release.yml`: on v* tags (and workflow_dispatch
  as an artifact-only dry run): build the bundle (gate), run a FRESH
  EVA proof on the runner via tools/eva_gate.py (gate), upload the
  artifact, and `gh release create` with all five evidence files
  attached. permissions: contents: write.
- pyproject version 0.0.1 -> 0.5.0; tag v0.5.0 exercises the pipeline
  for real (v0.5 milestone is complete; cFS example moved to v0.6).

## Local receipts before shipping

Builder ran end-to-end (captured t=3184 s, dv 12.38 m/s, barrier
verified) and the EVA gate re-proved 0 alarms into the same audit dir;
all five files present.

## Deferred

Release notes automation from the roadmap; attaching golden-vector
transcripts; signing the bundle (checksums/attestations).

## Push/merge instructions

Commit `31 — Per-tag release: audit bundle + fresh EVA proof (#31)`;
push; watch ci; then `git tag v0.5.0` + push the tag; watch the
release workflow create the evidence-backed release; close.
