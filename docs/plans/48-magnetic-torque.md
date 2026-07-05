# 48 — Magnetic disturbance torque + dipole field model

GitHub issue: https://github.com/adi-oltean/podium/issues/48

## Fix (landed)

Completes the canonical set of four environmental attitude
disturbances (gravity gradient #44, aero #45, SRP #46, magnetic #48).

- `attitude.magnetic_torque(dipole_body, b_body) = m x B` — residual
  (or magnetorquer-commanded) dipole crossed with the local field.
- `attitude.dipole_field(r_eci, dipole_moment, axis)` — centered
  dipole B = (mu0 m / 4pi r^3)(3(m_hat.r_hat)r_hat - m_hat), aligned
  (-z) by default; `constants.EARTH_DIPOLE_MOMENT = 7.94e22 A m^2`
  (equatorial surface field ~3.1e-5 T).
- `DisturbanceModel` gains an optional `residual_dipole` term and a
  `b_field_body` argument to `.torque`.

## Receipts (tests/test_magnetic_torque.py)

- tau = m x B EXACTLY; zero when the dipole is field-aligned.
- Dipole field ~3.1e-5 T at the equatorial surface and falls as
  1/r^3; the on-axis (pole) field is exactly 2x the equatorial value
  (the classic dipole factor).
- A ~1 A m^2 residual dipole in LEO gives a micro-Nm-class torque,
  comparable to the other disturbances.
- The aggregator includes the magnetic term.

## Note

magnetic_torque doubles as the magnetorquer ACTUATION model (a
commanded dipole producing control torque m x B) — useful for future
magnetic-control receipts.

## Deferred

A tilted/IGRF field for latitude structure; magnetorquer-based
detumbling (B-dot) control; combining the four torques in the
reference-mission attitude.

## Push/merge instructions

Single commit on main: `48 — Magnetic disturbance torque + dipole
field (#48)`; push; close.
