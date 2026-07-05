# 43 — Analytic torque-free attitude oracle

GitHub issue: https://github.com/adi-oltean/podium/issues/43

## Why analytic instead of tudatpy

The remaining v0.7 6-DOF-oracle item named tudatpy — but tudatpy is
itself a numerical integrator (and a heavy conda-only install).
Torque-free rigid-body rotation has EXACT closed-form solutions, so
they are a stronger, self-contained ground truth: the Euler+quaternion
RK4 is validated against mathematics, not against another integrator.

## Fix (landed) — `tests/test_attitude_analytic.py`

- Asymmetric body (I1<I2<I3): Euler's equations solve to Jacobi
  elliptic functions. With 2T (energy) and L^2 (momentum) from the
  initial state, omega = (a1 cn, a2 sn, a3 dn)(rate*t, m) with exact
  amplitudes/rate/modulus. Measured agreement with the integrator:
  ~2e-13 over 40 s.
- Axisymmetric body: closed-form regular precession. omega_3 constant,
  transverse omega rotates at the body rate lam = n(Ia-It)/It; the
  symmetry axis keeps a constant angle to the conserved inertial
  angular momentum H (nutation cone) and precesses about it at the
  inertial rate |H|/It. Validated against the integrated quaternion:
  nutation cos constant to 1e-9, H conserved to 1e-10, precession rate
  within 1%.
- Torque-free invariants (kinetic energy, |H_inertial|) conserved to
  1e-9 over 50 s — the assumptions the analytic solutions rest on.

## Deferred

A tudatpy (or Basilisk) numerical lane for COUPLED attitude +
translation under gravity-gradient/aero torques — a heavier install,
best via the Docker pattern (Orekit/CompCert/qemu); the analytic
oracle covers the torque-free rotational core exactly.

## Push/merge instructions

Single commit on main: `43 — Analytic torque-free attitude oracle
(#43)`; push; close.
