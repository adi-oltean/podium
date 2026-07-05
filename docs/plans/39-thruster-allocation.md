# 39 — Thruster/torque allocation for the 6-DOF PTR

GitHub issue: https://github.com/adi-oltean/podium/issues/39

## Fix (landed) — `podium.control.allocation`

The 6-DOF guidance (#33) commands a continuous body wrench (thrust
along +x_body, 3-axis torque). Hardware realizes it through N fixed,
PUSH-ONLY thrusters: thruster i produces force along a unit direction
d_i at body position r_i with magnitude u_i >= 0, so its wrench column
is [d_i ; r_i x d_i] and the cluster wrench is B u.

- `ThrusterConfig(positions, directions)` -> `effectiveness()` builds
  B (6 x N).
- `allocate(cfg, wrench, u_max)`: minimum-propellant LP
  (min sum u  s.t.  B u = w, 0 <= u <= u_max; scipy linprog "highs").
  On infeasibility, an NNLS fallback returns the closest achievable
  wrench and `feasible=False` with the residual — no silent lying.
- `standard_cluster(half)`: 24 thrusters, 3 inward-axis per cube
  corner; verified rank-6 B (full 6-DOF authority) — opposing corners
  give both force signs, corner offsets give both torque signs.

## Why non-negativity is the whole point

A plain pseudoinverse minimizes ||u|| but freely returns NEGATIVE
"thrust" no thruster can produce — `test_pseudoinverse_would_go_negative`
exhibits a realizable wrench whose min-norm solution is infeasible,
while the LP reproduces the same wrench with u >= 0.

## Receipts (tests/test_allocation.py)

- Cluster B is rank 6; directions unit; torque rows equal r x d.
- 200 random feasible wrenches reproduce to 1e-9 with u >= 0 and
  propellant no worse than the witness combination.
- A pure couple (torque, zero net force) is realized — the classic RCS
  requirement — with the net force cancelling to <1e-9.
- Beyond-authority demand flags infeasible and reports the closest
  achievable wrench + residual.
- Every node of a converged sixdof PTR plan's (thrust, tau) command
  allocates feasibly (slow): the guidance is hardware-realizable.

## Deferred

Path pointing/plume constraints in the planner; on-time quantization
(minimum-impulse-bit) of the allocated commands; wiring allocation into
the reference-mission endgame; wheel/CMG momentum management.

## Push/merge instructions

Single commit on main: `39 — Thruster/torque allocation (#39)`; push;
close.
