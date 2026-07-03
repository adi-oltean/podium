# 08 — Layer-0 convex guidance: DPP transcription, constraint library v0, ROE planner

GitHub issue: https://github.com/adi-oltean/podium/issues/8

## Problem

The library's third pillar (convex trajectory optimization) is design-only.
All prerequisites now exist: exact CW/YA/ROE discretizations, a closed-loop
engine to fly plans in, and specs to judge them.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/guidance/convex.py` | MISS | planners (sandbox layer: full Python + cvxpy) |
| `tests/test_convex.py` | MISS | receipts incl. engine integration |
| `pyproject.toml` | PARTIAL | cvxpy/clarabel into dev extras (opt extra exists) |
| `docs/roadmap.md` | PARTIAL | check off transcription/constraints items |

## Fix

`podium.guidance.convex` (sandbox side — prototyping layer; the embedded
path arrives with the v0.5 codegen):

1. `RendezvousPlanner`: impulsive planning on a fixed node grid, states
   propagated by exact per-interval STMs (CW constant; YA per-interval via
   anomaly propagation). Burns at every node incl. arrival. DPP-compiled
   once per (grid, constraint topology): STMs, boundary states, and KOZ
   normals are `cp.Parameter`s, so re-solves never rebuild. Objectives:
   L2 (sum of burn norms, SOCP) and L1 (LP-representable). Optional
   per-burn L2 cap.
2. Constraint library v0: approach cone (exact SOC) on selected nodes;
   rotating-hyperplane KOZ (linear; normals from a reference — straight-
   line by default, refined by one fixed re-solve pass with normals from
   the previous solution; hyperplane implies true distance >= R since
   ||r|| >= n^T r); plume half-space on burns at nodes whose reference
   range is inside the plume zone.
3. `RoePlanner`: ROE-space impulsive reconfiguration (states in ROE;
   dynamics roe+ = Phi (roe + Gamma dv) with Keplerian or J2 STM and the
   control-input matrix); terminal ROE equality; L1/L2 objective;
   post-solve passive-safety verification via `guidance.safety.rn_margin`.
4. `plan_to_controller`: adapter turning a plan into an engine controller.

Receipts:
- K=1 planner reproduces `cw.two_impulse` exactly (unique feasible point);
- multi-node L2 cost <= two-impulse cost (more freedom can't cost more);
- transcription exactness: replaying planned burns through the same STMs
  hits the target at solver tolerance; through the YA STM for e=0.15;
- planned trajectory flown closed-loop through the sim engine against the
  nonlinear truth, judged by specs (arrival, corridor);
- cone/KOZ/plume bite tests: unconstrained plan violates, constrained plan
  satisfies with margin, objective increases;
- DPP re-solve with new parameters returns the correct new solution;
- ROE plan reaches the target reconfiguration; rn_margin verification demo.

Deferred (filed in the issue): LCvx thrust-annulus + discrete-time
validity conditions; Breger-How passive-safety scenarios; convex e/i
passive-safety constraints.

## Acceptance Criteria

- [ ] Planners implemented; suite green (pytest, ruff, mypy)
- [ ] All receipts above green, incl. the engine integration flight
- [ ] Roadmap updated; follow-ups filed

## Push/merge instructions

Single commit on main: `08 — Layer-0 convex guidance (#8)`; push; close.

## Verification steps

Full suite; fly the planned approach in the engine and inspect spec
margins; compare planner dv against the glideslope baseline.
