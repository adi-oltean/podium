# 21 — PTR/SCvx* SCP core with exact-flow CTCS cuts

GitHub issue: https://github.com/adi-oltean/podium/issues/21

## Problem

Layer-0 approximates the keep-out sphere with reference hyperplanes.
Layer-1 must solve the true nonconvex problem, and hold constraints in
continuous time, not only at nodes.

## Fix (landed) — `podium.guidance.scp.PtrDockingPlanner`

- PTR loop: exact-STM dynamics, node KOZ re-linearized about the
  reference, virtual buffers with penalty, quadratic trust region.
- SCvx*-style updates: penalty ramp when infeasibility stalls with
  buffers active; trust-region EXPANSION when a feasible iterate's
  linearization validates (without it, nearly-flat fuel valleys force a
  1-(mu/w_tr) crawl — observed, diagnosed, fixed).
- Convergence: feasibility (buffers empty + clean dense check) AND
  stationarity, where fuel stationarity is accepted alongside small
  steps because min-fuel valleys let the state drift at constant cost.
- CTCS, exact-flow form: coast positions are linear in the decision
  variables, so continuous-time violations become exact linear cuts.
  Two implementation lessons paid for and encoded:
  1. cut TIMES persist, cut DIRECTIONS re-linearize each iteration —
     stale fixed-direction cuts from early references can permanently
     exclude the optimum (observed: slack pinned at 380+ while the true
     violation was zero);
  2. at most one new cut per arc per iteration, keyed by rounded time.

## Receipts (all green)

- convex reduction: no active KOZ -> Layer-0 optimum within 1e-4,
  <= 6 iterations;
- far-side passage: converged, buffers ~1e-9, true sphere satisfied at
  nodes AND on a 1000-sample independent dense check, cost within 0.1%
  of the two-pass hyperplane heuristic (trust-region bias, measured);
- coarse-grid dip: nodes clear the sphere but the coast arc dips inside
  between nodes — cuts catch it (n_cuts > 0) and the final dense check
  is clean;
- penalty ramp engaged from w_pen0 = 1e-6 and still converged;
- bitwise deterministic across runs;
- engine flight: min range >= R - 2 m against the nonlinear truth.

## Deferred

6-DOF attitude-coupled PTR (with MuJoCo contact, next work item),
integral-augmentation CTCS for non-coast dynamics, state-triggered
constraints, STL-in-SCP smooth robustness encodings.

## Push/merge instructions

Single commit on main: `21 — PTR/SCvx* SCP core (#21)`; push; close.
