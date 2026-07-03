# 16 — Attitude dynamics + quaternion-feedback controller

GitHub issue: https://github.com/adi-oltean/podium/issues/16

## Problem

Docking acceptance (IDSS box) needs angular rates and misalignment;
nothing rotational existed beyond the quaternion kernel.

## Fix (landed)

- `podium.dynamics.attitude` (truth): Euler equations with full inertia
  tensor, quaternion kinematics from `core.quat`, coupled RK4 with
  per-step renormalization; kinetic-energy and inertial-angular-momentum
  helpers for receipts.
- `podium.control.attitude` (flight, static subset): classical
  quaternion-feedback regulator, shortest-way error from `quat.error`,
  per-axis saturation, contracts on gains.

Receipts: torque-free intermediate-axis tumble conserves energy and
inertial |L| to 1e-9 relative over 10 min with |q|=1 (and visibly
exhibits the instability); detumble 10 deg/s -> <0.01 deg/s; 20-degree
slew with gains designed for wn=0.1/zeta=0.9 shows the predicted
overshoot bound and 2% settling; saturation held.

## Acceptance Criteria

- [x] Suite green (142 tests, ruff, mypy)
- [x] Roadmap + package docstrings refreshed

## Push/merge instructions

Single commit on main: `16 — Attitude dynamics + quaternion feedback
(#16)`; push; close.
