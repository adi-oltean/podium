# 14 — Relative-nav EKF (Joseph form, static subset)

GitHub issue: https://github.com/adi-oltean/podium/issues/14

## Problem

v0.3 "Full loop" needs navigation: every closed-loop demo so far flies
on perfect (or perfectly-noisy) measurements fed straight to control.
No estimator exists.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/nav/ekf.py` | MISS | predict / Joseph update kernels + CW wrapper |
| `tests/test_ekf.py` | MISS | receipts |
| `src/podium/nav/__init__.py` | PARTIAL | exports |
| roadmap | PARTIAL | check off |

## Fix

Kernels (static-subset style: pure, fixed shapes, closed-form):
- `predict(x, p, phi, q)`: STM propagation of state and covariance.
- `update_joseph(x, p, z, h, r)`: gain via 3x3/6x6 solves, state update,
  covariance in Joseph form P = (I-KH) P (I-KH)' + K R K' (symmetric
  and PSD-preserving under roundoff, unlike the naive form), returns
  innovation and innovation covariance for consistency monitoring.
- `process_noise_wna(dt, q_accel)`: discrete white-noise-acceleration
  Q (per-axis blocks [[dt^3/3, dt^2/2], [dt^2/2, dt]] * q_accel).
- Measurement matrices: `H_POS` (position-only), `H_POSVEL`.

Sandbox wrapper `RelNavEkf(n, dt, q_accel, r_pos, r_vel=None)`: holds
(x, P), CW STM prediction, position or pos+vel updates; `step(z, dv)`
handles the impulsive-burn feedthrough (dv adds to the velocity estimate
at the tick — burns are known inputs, not disturbances).

Receipts:
- Joseph invariants: 500-step random-model sequence keeps P symmetric to
  1e-12 and PD.
- Consistency vs the engine: unforced scenario with seeded measurement
  noise; time-averaged NEES within the chi-square 95% band for 3/6 dof
  and NIS likewise (q_accel absorbs the CW-vs-nonlinear model error).
- Convergence: 100 m / 0.5 m/s initial estimate error collapses below
  measurement noise within the transient.
- Closed loop: LQR on the EKF estimate (position-only measurements,
  sigma = 5 m) stabilizes the nonlinear truth through the engine.

## Acceptance Criteria

- [ ] Suite green (pytest, ruff, mypy)
- [ ] NEES/NIS receipts pass with documented q_accel choice
- [ ] Closed-loop nav+control receipt green
- [ ] Roadmap updated

## Push/merge instructions

Single commit on main: `14 — Relative-nav EKF (#14)`; push; close.
