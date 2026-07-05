# 50 — Certificate-carrying Lyapunov ellipsoid invariants

GitHub issue: https://github.com/adi-oltean/podium/issues/50
Realizes recorded idea R2 (Feron-style credible autocoding).

## Fix (landed)

- `lqr.dlqr_cert(a, b, q, r)` returns (K, P): the gain AND the
  value-function (Riccati) matrix, which is the closed loop's Lyapunov
  certificate.
- `podium.verify.lyapunov`: `EllipsoidInvariant(P)` and
  `verify_lyapunov(a_cl, p)` re-verify, in EXACT `fractions.Fraction`
  arithmetic (no floats in the trusted path, reusing
  `barrier.is_psd`'s all-principal-minors test): P >= 0 (the ellipsoid
  {x'Px<=c} is bounded) and the Lyapunov decrease
  P - A_cl' P A_cl >= 0 (x'Px non-increasing along the closed loop).

## Why exact verification is robust here

For the LQR value matrix, P - A_cl'PA_cl = Q + K'RK, which exceeds Q >= 0
by the full state-cost margin. Rationalizing the float Riccati solution
perturbs the inequality by ~1e-12, far below the Q eigenvalues, so the
exact PSD check passes comfortably — the certificate is robust, not
knife-edge (unlike the hand-built barrier certificate, which sits on
its constraint boundary).

## Receipts (tests/test_lyapunov.py)

- The CW LQR value matrix certifies exactly (P >= 0 and decrease >= 0).
- x'Px is monotonically non-increasing along a 300-step closed-loop
  trajectory and converges toward zero; the sublevel set {x'Px<=c} is
  invariant over 500 steps.
- Negatives: an identity P is rejected, and the OPEN-loop CW plant
  (not contracting) fails the Lyapunov decrease against the closed-loop
  P. The EllipsoidInvariant value x'Px is computed exactly.

## The exact-arithmetic certificate family

With the abort-safety barrier (#20) and the online-solver KKT
certificates (#40/#41), this is the control-Lyapunov member: three
independent safety/stability/optimality certificates all re-verified
over the rationals. The C emitter can render P as a quadratic
PROVE/ACSL obligation on the flight controller (the credible-autocoding
endpoint).

## Deferred

Continuous-time (CARE) Lyapunov certificate; strict-PD margin via
P - eps*I; emitting the quadratic ACSL obligation and discharging it
with Frama-C/WP; robust invariance under bounded disturbances (the
ultimate-bound ellipsoid).

## Push/merge instructions

Single commit on main: `50 — Certificate-carrying Lyapunov ellipsoid
(#50)`; push; close.
