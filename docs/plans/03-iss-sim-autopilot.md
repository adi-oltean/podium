# 03 — Podium autopilot for the SpaceX ISS docking sim

GitHub issue: https://github.com/adi-oltean/podium/issues/3

## Problem

Requested live demo: Podium GNC flying https://iss-sim.spacex.com/ (SpaceX's
public Crew Dragon docking simulator). We cannot host a modified copy of
their app; the honest integration is a **console autopilot**: a script the
user pastes into the browser devtools on the sim page. It reads the sim's
HUD telemetry from the DOM and actuates by clicking the sim's own control
buttons — Podium's control law in the loop, their vehicle model.

Sim interface (recovered from page source):
- Telemetry: `#x-range/#y-range/#z-range .distance` (m), `#pitch/#roll/#yaw
  .error` (deg) and `.rate` (deg/s), `#range .rate` (m), `#rate .rate` (m/s).
- Actuation: `#translate-{left,right,up,down,forward,backward}-button`,
  `#{yaw,pitch,roll}-{left,right,up,down}-button` — each click is a fixed
  velocity/rate increment (pulsed, minimum-impulse-bit actuation).
- Success criteria (from the sim): |lateral| < 0.2 m, |angles| < 0.2 deg,
  closing rate < 0.2 m/s at contact.

## Fix

1. `src/podium/control/docking.py` (static-subset style, contracted):
   - `approach_rate(range_m, ...)` — range-scheduled closing-rate profile
     (fast far out, ramping to a soft terminal rate inside the capture zone).
   - `axis_rate_cmd(pos, v_max, k)` — saturated proportional velocity command
     for lateral position nulling.
   - `pulses_needed(v_cmd, v_est, dv_pulse)` — minimum-impulse-bit allocator:
     integer click count with deadband (round-to-zero inside half a pulse).
2. `tests/test_docking.py` — closed-loop click-dynamics simulation per axis
   (double integrator actuated only by fixed-increment pulses at 1 Hz):
   assert docking criteria met from representative initial conditions, no
   overshoot through the port, monotone terminal approach.
3. `viewer/iss-sim/index.html` — demo page: what it is, how to run it
   (open sim, begin, paste script in console), the autopilot JS (ported
   line-for-line from docking.py, constants identical), troubleshooting.
   JS estimates velocities/rates by finite-differencing DOM telemetry, then
   applies the same pulse allocator. Parity: constants and formulas mirrored;
   fermi-style discipline (Python is the source of truth).
4. Link from the main viewer page and README.

## Tests

`tests/test_docking.py` (new). Browser-side script verified manually against
the live sim (record outcome in the issue).

## Acceptance Criteria

- [ ] docking.py + tests green (pytest, ruff, mypy)
- [ ] Demo page live on Pages
- [ ] Autopilot verified to dock successfully on iss-sim.spacex.com

## Push/merge instructions

Single commit on main: `03 — Podium autopilot for the SpaceX ISS docking sim
(#3)`, push (triggers Pages deploy), close #3 after manual verification.

## Verification steps

Run the suite; open the demo page; follow its instructions on
iss-sim.spacex.com and confirm SUCCESS screen.
