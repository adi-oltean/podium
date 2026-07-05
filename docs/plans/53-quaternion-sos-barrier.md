# 53 — Exact SOS barrier for the quaternion-feedback closed loop

GitHub issue: https://github.com/adi-oltean/podium/issues/53

## Motivation

The SOS thread's remaining open piece (after #51 checker, #52 validated
rounding) was problem-specific: apply the machinery to the polynomial
RPOD closed loop the paper names --- quaternion feedback. This issue
does exactly that.

## The system

Single-axis attitude with quaternion feedback, state (q_s, q_v, w),
torque tau = -kp q_v - kd w:

    q_s' = -1/2 w q_v,   q_v' = 1/2 w q_s,   w' = (1/I)(-kp q_v - kd w).

A genuinely nonlinear (quadratic) polynomial vector field on the
unit-quaternion manifold q_s^2 + q_v^2 = 1, and RPOD-core (attitude
control for docking).

## Result

The Lyapunov barrier V = 2 kp (1 - q_s) + 1/2 I w^2 has

    V-dot = -kd w^2

EXACTLY: the kp q_v w terms from the kinematics and the dynamics cancel,
verified over the rationals for arbitrary rational gains (the same kind
of exact cancellation as the Duffing case, on a real controller). So
-V-dot = kd w^2 is SOS and every sub-level set {V <= c} is an
infinite-horizon attitude-stability invariant of the nonlinear system.

## Receipts (tests/test_quaternion_barrier.py)

- The cross-term cancellation V-dot = -kd w^2 is exact for several gain
  triples.
- -V-dot certifies SOS (basis [w]).
- The certificate synthesized by an UNTRUSTED Clarabel SDP over a richer
  basis [w, q_v, q_s] is validated to an exact rational certificate via
  validate_gram (#52) and re-checked by is_sos --- the full pipeline on
  the real closed loop.
- A simulation confirms the loop converges (q_v, w -> 0), V is monotone,
  and the quaternion norm q_s^2 + q_v^2 is preserved (the certified
  invariant is physical).

## Significance

A real RPOD nonlinear closed loop now runs through the full exact-SOS
pipeline --- synthesize with an untrusted float SDP, then verify the
answer in exact arithmetic --- the same discipline as the barrier
(#20), KKT (#40/#41), Lyapunov (#50), and SOS (#51/#52) certificates.
The remaining open direction is a HIGHER-degree barrier (proving a
non-trivial safe cone, needing the off-diagonal validated synthesis of
#52) for the full three-axis quaternion loop.

## Push/merge instructions

Single commit on main: `53 — Exact SOS barrier for the
quaternion-feedback closed loop (#53)`; push; close.
