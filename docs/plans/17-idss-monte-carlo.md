# 17 — IDSS docking acceptance + Monte Carlo campaigns

GitHub issue: https://github.com/adi-oltean/podium/issues/17

## Problem

No acceptance criterion tied docking sims to a real interface standard,
and there was no campaign machinery for dispersed runs.

## Fix (landed)

- `podium.sim.idss`: IdssBox (IDD Rev G contact conditions, SI) +
  `check_translation` (closing window, lateral offset/rate about the
  approach axis) and `check_attitude` (misalignment, angular rate) as
  margin dictionaries.
- `podium.sim.monte_carlo`: `run_campaign(n, master_seed, make_case)` —
  one master Generator spawns per-run seeds (engine noise and case
  dispersions share the run seed); structured-array output includes the
  per-run seed so any run replays for post-mortem; `summary()` helper.

Receipts: single-run terminal approach (rate-command profile tapering
to 0.075 m/s) through the engine with 2 cm/3 mm/s sensor noise, 1 mm/s
MIB, and 1% execution error contacts inside all four translational
margins; attitude hold from 10 deg / 0.5 deg/s satisfies the rotational
box well before any plausible contact time; box constants pinned to the
Rev G numbers; 20-run dispersed campaign 100% in-box with the whole
metric table bit-identical across repeats and different under a new
master seed.

Decoupling note: translation flies in the engine, rotation in the
attitude propagator; they are coupled only at the contact instant. The
honest 6-DOF coupling (thrust in body frame, torque allocation) is the
v0.4 MuJoCo/6-DOF layer.

## Acceptance Criteria

- [x] Suite green (146 tests, ruff, mypy)
- [x] Roadmap/README/sim docstring refreshed

## Push/merge instructions

Single commit on main: `17 — IDSS acceptance + Monte Carlo (#17)`;
push; close.
