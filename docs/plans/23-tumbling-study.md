# 23 — Tumbling-target terminal guidance (scoped study)

GitHub issue: https://github.com/adi-oltean/podium/issues/23

## Scope

First-principles study + working planner for terminal approach to a
target tumbling at a KNOWN rate (planar tumble about the cross-track
axis, port at radius rho_p). Uncertain/torque-coupled tumbles are the
follow-on (they are where SCP re-enters).

## Structural result

Known tumble => the port state p(t), pdot(t) is deterministic
kinematics, so the planning problem stays CONVEX on the exact CW STM:
- terminal capture = ordinary boundary condition (position AND velocity
  of the rotating port at t_f — grapple/berth condition);
- rotating approach corridor = per-node second-order cones with KNOWN
  time-varying axes (the port direction at each node time).
`podium.guidance.tumbling.plan_tumbling_dock` implements exactly this.

## Study findings (each pinned by a test)

1. **Fuel vs rate is phase-confounded.** With a free tumble phase, the
   fuel-vs-rate curve is non-monotone: the port's orientation at
   arrival changes with the rate, and arrival geometry dominates the
   co-rotation cost. Clean envelopes must FIX the arrival phase per
   rate (phase0 = phi_arrival - w * t_f). Pinned by
   test_free_phase_is_nonmonotone_study_finding.
2. **The co-rotation intuition fails at low rates.** Even phase-fixed,
   dv(0.002 rad/s) < dv(0): matching a slowly-moving port is CHEAPER
   than nulling all relative motion, because CW natural drift supplies
   arrival velocity for free. "dv grows by rho*w" is not a theorem for
   the full trajectory problem. Pinned as the dip assertion in
   test_envelope_findings_and_closure.
3. **Measured envelope** (150 m approach, 480 s, 25 nodes, rho_p = 10 m,
   corridor 20 deg over the final half, dv cap 0.35 m/s per burn,
   arrival phase -pi/2):

   | w_spin [rad/s] | dv [m/s]  |
   |---------------:|-----------|
   | 0.000          | 0.629     |
   | 0.002          | 0.601 (dip) |
   | 0.005          | 1.211     |
   | 0.010          | 1.855     |
   | 0.020          | 2.140     |
   | 0.050          | infeasible |
   | 0.200          | infeasible |

   ~1 deg/s (0.02 rad/s) is flyable at ~2.1 m/s; ~3 deg/s closes the
   envelope under this burn cap — consistent with the flight-community
   rule of thumb that few-deg/s tumbles need dedicated capture systems.

## Receipts

Port kinematics exact (radius/speed/perpendicularity); planner matches
the rotating port to 1e-5 m / 1e-6 m/s through its own dynamics;
corridor cones verified at nodes; engine flight against the nonlinear
truth arrives ON the independently-recomputed rotating port within
1 m / 1 cm/s; envelope table regression-pinned.

## Follow-on scope (not this study)

Uncertain tumble (estimator in the loop, robust corridors), 3-D
tumbles (nutation; port cone axis no longer planar), torque-free
Euler propagation of the target from podium.dynamics.attitude, plume
constraints against the rotating body, SCP once any of these break
convexity.

## Push/merge instructions

Single commit on main: `23 — Tumbling-target study (#23)`; push; close.
