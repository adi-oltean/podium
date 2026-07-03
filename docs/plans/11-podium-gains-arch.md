# 11 — Podium-synthesized gains through the ARCH reachability gate

GitHub issue: https://github.com/adi-oltean/podium/issues/11

## Problem

The reachability gate (#10) proves the benchmark's *published* reference
gains. The point of R1 is code<->model traceability for OUR controller:
gains synthesized by the library, safety re-proven by the gate. Until the
gate proves a Podium-synthesized controller, the story is a demo, not a
workflow.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/control/lqr.py` | PARTIAL | add continuous CARE solver + clqr gain |
| `src/podium/guidance/arch.py` | PARTIAL | podium-gains closed-loop model + export variant |
| `tests/test_arch.py`, `tests/test_guidance_control.py` | PARTIAL | synthesis + margin receipts |
| `tools/reach` / `.github/workflows/reach.yml` | PARTIAL | prove the podium variant too |

## Fix

1. `care(a, b, q, r)` in `podium.control.lqr` (sandbox side): continuous
   algebraic Riccati via the Hamiltonian stable-subspace method (numpy
   eig on the 2n Hamiltonian; P = X2 X1^-1 symmetrized); `clqr` returns
   the gain. Receipts: Riccati residual ~ machine precision; closed loop
   Hurwitz.
2. `arch.podium_variant()`: plant = planar CW at the benchmark's GEO mean
   motion (from `core.cw` structure, not transcribed); per-mode
   acceleration gains from `clqr` with documented Q/R choices (approach
   mode weighted to squeeze cross-track error before the corridor
   handoff; attempt mode heavily velocity-damped for the octagon);
   closed-loop A matrices assembled as A_cw - B K. Same FSM, guards,
   initial set, and properties as the reference model.
3. Export `arch_srna01_podium.json` / `arch_sra01_podium.json`; the gate
   proves four models. Local PROVEN before shipping (iterate Q/R against
   the 12 s proof loop if needed).
4. Document the implied acceleration-gain magnitudes of the published
   controller (recovered as A_published - A_cw) for comparison.

## Acceptance Criteria

- [x] CARE receipts green (residual < 1e-10 scaled, P symmetric positive
      definite); gains differ from published at rtol 1e-3 (genuine
      synthesis); closed loops Hurwitz
- [x] Simulation margins positive from all initial-set corners, both
      scenarios
- [x] JuliaReach PROVEN locally for both podium-variant scenarios
      (~13 s total reach time, first attempt — the Q/R choices below
      landed in the verified regime without iteration)
- [x] CI gate extended to four models; run green (verified post-push)
- [x] Roadmap updated

Q/R record: Q_approach = diag(0.0033, 0.0044, 8.3, 8.3),
Q_attempt = diag(0.33, 0.33, 369, 369), R = I (per-minute units) —
chosen via the double-integrator asymptotes so closed-loop
bandwidth/damping land in the regime of the published controller
(position gains ~0.06 / 0.58 per min^2, velocity damping ~2.9 / 19 per
min), with the y-weight in approach mode raised to squeeze cross-track
error before the corridor handoff.

## Push/merge instructions

Single commit on main: `11 — Podium-synthesized gains through the ARCH
gate (#11)`; push (triggers reach workflow); confirm green; close.

## Verification steps

`julia --project=tools/reach tools/reach/arch_rendezvous.jl` over all
four exports; CI reach run green.
