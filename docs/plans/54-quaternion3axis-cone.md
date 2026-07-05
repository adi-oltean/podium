# 54 — Three-axis quaternion safe-attitude-cone barrier

GitHub issue: https://github.com/adi-oltean/podium/issues/54

## The frontier

After #51 (checker), #52 (validated rounding), and #53 (single-axis
quaternion barrier), the last open SOS direction was a higher-degree
barrier for the FULL three-axis quaternion loop proving a non-trivial
safe attitude cone, not just global convergence. This issue delivers it
in closed form, so no blind 7-state SDP is needed.

## The system

Full three-axis attitude closed loop, 7 states (q0..q3, w1..w3),
asymmetric inertia, torque tau = -kp q_vec - kd w. Genuinely nonlinear
polynomial dynamics on the unit-quaternion manifold.

## Two exact results

1. **V-dot = -kd ||w||^2.** For V = 2kp(1-q0) + 1/2 w' I w, the Lie
   derivative is exactly -kd(w1^2+w2^2+w3^2): the kp coupling from the
   kinematics cancels, and the gyroscopic w x I w term is workless
   (w . (w x I w) = 0), so it cancels too. Verified over the rationals
   for several asymmetric inertia/gain sets.
2. **Safe-cone Positivstellensatz.** {V <= c} subset {q0 >= q0_min}
   (c = 2kp(1 - q0_min)) via the EXACT identity
   q0 - q0_min = (1/2kp)(c - V) + (1/4kp) w' I w, whose remainder
   (1/4kp) w' I w is SOS. On {V <= c} (c - V >= 0), q0 - q0_min is a
   sum of two nonnegative terms, so q0 >= q0_min.

Together: starting inside {V <= c} the attitude stays in the safe cone
for all time --- an infinite-horizon safe-cone certificate for the full
three-axis nonlinear loop.

## Receipts (tests/test_quaternion3axis_cone.py)

- The gyroscopic + kp cancellation V-dot = -kd||w||^2 is exact for
  three inertia/gain triples.
- -V-dot is SOS.
- The cone-containment remainder equals (1/4kp) w' I w exactly and is
  SOS.
- The remainder's SOS certificate synthesized by an untrusted Clarabel
  SDP is validated to an exact rational certificate (validate_gram, #52).
- A 7-state RK4 simulation confirms the attitude keeps q0 >= q0_min for
  all time, preserves ||q||, and converges.

## Significance

This closes the SOS thread: an exact, infinite-horizon, safe-attitude-
cone barrier for the full three-axis quaternion-feedback closed loop,
synthesized (untrusted SDP) and verified in exact arithmetic. The
paper's last SOS future-work sentence is now a shipped RPOD result.

## Deferred

Cone barriers that are not sub-level sets of the natural Lyapunov
function (needing genuine multi-variable SOS multipliers with
off-diagonal Gram freedom); coupling to the reference-mission attitude.

## Push/merge instructions

Single commit on main: `54 — Three-axis quaternion safe-cone barrier
(#54)`; push; close.
