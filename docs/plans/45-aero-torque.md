# 45 — Aerodynamic disturbance torque + weathervane stability

GitHub issue: https://github.com/adi-oltean/podium/issues/45

## Fix (landed)

`podium.dynamics.attitude.aerodynamic_torque(v_rel_body, rho, cd_area,
r_cp)` = r_cp x F_drag with F_drag = -1/2 rho (Cd A) |v_rel| v_rel —
the second dominant LEO attitude disturbance after gravity gradient
(#44). v_rel_body is the atmosphere-relative velocity in body
coordinates; r_cp the center-of-pressure offset from the c.m.

## Receipts (tests/test_aero_torque.py)

- tau = r_cp x F_drag EXACTLY, F_drag matching the standard quadratic
  drag law.
- Zero torque when the cp lies on the velocity line (no perpendicular
  moment arm).
- WEATHERVANE stability: cp behind the c.m. (downstream) -> a +alpha
  angle-of-attack gives a torque about +z that OPPOSES it (restores
  toward the flow); the cp-ahead configuration destabilizes (sign
  flips). Passive aerodynamic stability, the aero analogue of the
  gravity-gradient restoring in #44.
- Force consistency: the aero force equals the truth model's drag
  acceleration times mass (podium.dynamics.nonlinear, cd_area = m/bc)
  to 1e-12 relative.

## Deferred

Coupled aero+gravity-gradient attitude propagation with the rotating
atmosphere in the orbit loop; free-molecular vs continuum Cp models;
solar-radiation-pressure torque; wiring environmental torques into the
reference-mission attitude.

## Push/merge instructions

Single commit on main: `45 — Aerodynamic disturbance torque (#45)`;
push; close.
