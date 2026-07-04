# 33 — 6-DOF attitude-coupled PTR (body-fixed thruster)

GitHub issue: https://github.com/adi-oltean/podium/issues/33

## Fix (landed) — `podium.guidance.sixdof`

13-state joint planning (r, v LVLH/CW; q body->LVLH; w body rates),
4 controls (throttle of one +x-body thruster, 3-axis torque). PTR:
central-difference discrete Jacobians around the RK4-propagated,
quaternion-renormalized nonlinear reference; virtual-control slack
with penalty ramp; HARD per-block trust-region boxes with
accept/reject and contract-on-stall; Clarabel with SCS fallback.
Converges on the braking-slew scenario (30 m -> 2 m V-bar, terminal
attitude retro) in 22 iterations: slack ~1e-16, defect 2.7e-5, fuel
45.25 N s, physical two-burn profile.

## Three measured lessons (each cost one debugging round)

1. **The bilinearity vanishes at lazy references.** At a zero-thrust,
   constant-attitude reference, d(T R(q)e1) has no dq component (T=0)
   and dT enters only along the frozen axis — the linear model cannot
   produce off-axis force, and PTR wedges at irreducible slack no
   matter the penalty. The initial reference must exercise the
   coupling: SLERP attitude + a small thrust seed. (Diagnosed by
   printing the B thrust columns: they only sweep when q_ref slews.)
2. **The scenario itself was backwards.** First terminal attitude
   pointed the thruster PROGRADE — a braking approach cannot end
   that way, and the min-slack probe isolated the deficit to exactly
   the terminal v_y row (0.097 m/s). Physics: camera on the -x body
   face looks at the target while the +x thruster brakes. Also the
   TOF must be consistent with decel-only closing.
3. **Trust-region direction matters by regime.** Soft quadratic
   weights alone let the subproblem take steps that blow up the
   nonlinear attitude propagation (NaN reference). And scp.py's
   decay-on-progress rule is BACKWARDS near the solution here — it
   produced a period-2 limit cycle at defect 0.465 (measured,
   alternating fuel 44.33/44.60). Hard boxes + contract-when-defect-
   stalls converges cleanly.

## Receipts (tests/test_sixdof.py, all green)

Independent nonlinear replay of the planned controls hits the
terminal box; the planner DISCOVERS the ~90 deg slew and the final
burn's thrust axis is retro (axis_y < -0.7) — nobody encoded a slew,
only the terminal attitude box and the body-fixed thruster; all
bounds respected; thrust acceleration is exactly (T/m) R(q)e1 in the
dynamics; quaternion norms preserved to 1e-9.

## Deferred

Torque allocation to discrete thruster pairs; plume/pointing PATH
constraints (terminal-only today); use in the reference mission's
endgame; contact-attitude coupling into the MuJoCo handoff.

## Push/merge instructions

Single commit on main: `33 — 6-DOF attitude-coupled PTR (#33)`;
push; close.
