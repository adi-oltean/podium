# 15 — Sensor models + actuator MIB/execution error

GitHub issue: https://github.com/adi-oltean/podium/issues/15

## Problem

Closed-loop demos measured the true state plus white noise and executed
burns perfectly. Real RPOD flies on bearing/range sensors and imperfect
thrusters.

## Fix (landed)

- `podium.nav.sensors`: RelGnss (white noise + seeded constant bias),
  DockingCamera / Lidar (az/el/range of the target from the chaser,
  proportional range noise, visibility gating). Flight-side `camera_h`
  and `camera_jacobian` (FD-pinned).
- `podium.nav.ekf.update_joseph_nonlinear`: EKF update for nonlinear
  measurements, Joseph covariance, angle-row innovation wrapping.
- Engine actuator model (Scenario fields): per-tick magnitude cap ->
  MIB quantization -> seeded proportional execution error, in that
  order; the burn log records what was actually applied.

Receipts: sensor statistics match budgets (4000-draw checks incl.
recovered bias); camera Jacobian vs FD; camera-only EKF through the
engine converges from 50+ m error to <6 m RMS at ~800 m range; MIB/cap
arithmetic visible in the burn log (sub-MIB burns vanish); the physics
lesson — 2% execution error on a 1.4 m/s open-loop insertion burn
misses by >100 m (3*dv*t drift) while LQR feedback through the same
noisy actuator converges to <5 m; bit-identical replay with the
actuator noise seeded.

## Acceptance Criteria

- [x] Suite green (138 tests, ruff, mypy)
- [x] Roadmap/README/architecture refreshed (docs-upkeep policy)

## Push/merge instructions

Single commit on main: `15 — Sensor models + actuator imperfections
(#15)`; push; close.
