# 46 — Solar-radiation-pressure torque + eclipse cutoff

GitHub issue: https://github.com/adi-oltean/podium/issues/46

## Fix (landed)

`podium.dynamics.attitude.srp_torque(sun_dir_body, area, cr, r_cp,
pressure=SOLAR_PRESSURE, illuminated=True)` = r_cp x F_srp with
F_srp = -P Cr A s_hat — the third classic environmental attitude
disturbance, dominant at GEO. s_hat is the unit direction TO the Sun
in body coordinates (force anti-sunward, hence the minus); Cr in [1,2]
the reflectivity; illuminated=False models the eclipse cutoff.
`constants.SOLAR_PRESSURE = 4.5606e-6 N/m^2` (1-AU solar irradiance /
c) added as the single source of truth.

## Receipts (tests/test_srp_torque.py)

- tau = r_cp x F_srp EXACTLY.
- Force anti-sunward with magnitude P Cr A.
- Zero when the cp lies on the sun line.
- Eclipse (illuminated=False) zeroes the torque.
- Reflectivity scales the force linearly: a perfect reflector (Cr=2)
  pushes exactly twice as hard as a perfect absorber (Cr=1).

## The environmental-torque trio

Gravity gradient (#44), aerodynamic (#45), and SRP (#46) are the three
dominant environmental attitude disturbances, now all modeled with
analytic-validated receipts. Gravity gradient and aero dominate in
LEO; SRP dominates at GEO.

## Deferred

Distance scaling of P for non-1-AU orbits; a cannonball-plus-panels
multi-surface SRP model; a geometric eclipse (cylindrical/conical
shadow) predicate; combining all three torques in an attitude
propagation with the reference mission.

## Push/merge instructions

Single commit on main: `46 — SRP torque + eclipse cutoff (#46)`;
push; close.
