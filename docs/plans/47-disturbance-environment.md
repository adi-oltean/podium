# 47 — Disturbance environment aggregator + attitude-hold rejection

GitHub issue: https://github.com/adi-oltean/podium/issues/47

## Fix (landed) — `podium.dynamics.disturbances`

`DisturbanceModel` holds the fixed spacecraft geometry (inertia,
orbital rate, aero cd_area + cp, SRP area/Cr/cp) and `.torque(
nadir_body, v_rel_body, rho, sun_body, illuminated)` returns the total
body-frame environmental torque. Torques superpose (independent
physical effects), so the aggregate is the sum of gravity gradient
(#44), aerodynamic (#45), and SRP (#46); any term disables by leaving
its config None (or rho=0 for aero).

## Receipts (tests/test_disturbances.py)

- Superposition: total torque == gravity-gradient + aero + SRP
  exactly.
- Terms disable cleanly (gravity gradient always contributes).
- The three magnitudes are physically reasonable (micro-Nm-class) for
  a representative LEO spacecraft.
- CAPSTONE: a quaternion-feedback controller holds an inertially-fixed
  attitude within 1 degree over a full orbit against the combined
  disturbance — nadir and the relative wind rotate at the orbital
  rate, the Sun is inertially fixed — and the control never exceeds
  saturation. Connects the disturbance DYNAMICS to disturbance
  REJECTION, the proximity-ops attitude-hold problem.

## Deferred

Wiring the disturbance model into the reference-mission attitude and
the 6-DOF PTR truth; a magnetic-dipole torque as a fourth term; a
geometric eclipse predicate feeding `illuminated`; worst-case
disturbance-torque bounds for control sizing.

## Push/merge instructions

Single commit on main: `47 — Disturbance environment + rejection
(#47)`; push; close.
