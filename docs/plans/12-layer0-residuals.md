# 12 — Layer-0 residuals: eccentric finite-burn, primer certificate, coast arcs, MIB bridge

GitHub issue: https://github.com/adi-oltean/podium/issues/12

## Problem

Five residuals from #9: finite-burn planning was CW-only; the LCvx audit
was a posteriori with no dual-side normality signal; coast-heavy
problems had no supported path; passive-safety drift coverage between
sparse samples was unverified; impulsive plans ignored thruster
minimum-impulse-bit quantization.

## Fix (all landed)

1. `control.lqr.ya_discrete`: YA ZOH (Ad exact STM, Bd by Simpson).
   Receipts: equals cw_discrete at e=0; two-interval composition
   identity Ad2@Bd1+Bd2 == Bd_full. FiniteBurnPlanner takes e/theta0
   with per-interval Parameter maps; controllability via the lifted
   time-varying reachability matrix; eccentric plan replays exactly.
2. Primer normality certificate from the dynamics-constraint duals:
   primer_k = Bd_k' lambda_k; Gamma-stationarity gives ||primer|| == dt
   at tight interior-slack nodes (measured exactly — pins the dual
   convention), so primer/dt is the scale-free certificate: >0.01 on
   the normal problem, <1e-4 (measured ~1e-7) on the degenerate one.
3. `find_min_time`: bounded bisection on feasibility (12 iterations);
   solve() now returns an infeasible-status plan instead of asserting.
   Receipt: 0.9x t_min infeasible; LCvx at 1.15x t_min with a 30%
   throttle floor passes the losslessness audit clean.
4. Dense passive-safety verification: Plan.ps_margins reports the TRUE
   minimum drift distance (200-sample grid) minus the radius per
   failure node — inter-sample dips become visible numbers (test allows
   a small documented dip, > -10 m, vs the 200 m radius).
5. `quantize_plan`: per-axis sigma-delta thruster-click quantization;
   cumulative residual <= q/2 per axis; the quantized plan flown through
   the engine against the nonlinear truth still arrives (q = 0.002 m/s
   ~ 1 N s on a 500 kg vehicle; measured arrival error grows ~4 m over
   the unquantized baseline).

## Acceptance Criteria

- [x] All receipts green (suite 127 tests, ruff, mypy)
- [x] Roadmap updated

## Push/merge instructions

Single commit on main: `12 — Layer-0 residuals (#12)`; push; close.
