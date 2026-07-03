# 20 — SOS abort-safety certificates (exact-rational verification)

GitHub issue: https://github.com/adi-oltean/podium/issues/20

## Problem

The reachability gate proves bounded-horizon safety. Abort safety is an
infinite-horizon claim: after a failure, the free drift must stay clear
of the keep-out zone forever, not just for the analyzed window.

## Design decisions (the ones that made it exact)

1. **Time-scaled coordinates**: with u = (r, v/n) and tau = n t, the CW
   matrix is INTEGER — the entire certificate problem is rational, so
   the trusted checker runs in `fractions.Fraction` with zero floats.
2. **Invariant-basis barrier**: B is a combination of CW flow invariants
   (c1 = 4X + 2VY, in-plane amplitude^2, cross-track amplitude^2), so
   dB/dtau = 0 STRUCTURALLY; the checker verifies A'P + PA == 0 exactly
   and conservation gives invariant sublevel sets — all-time safety with
   no discretization or horizon.
3. **RN-plane keep-out**: ||(x,y,z)|| >= ||(x,z)||, so certifying
   radial/cross-track separation covers any along-track drift — the
   e/i-vector heritage argument, machine-checked.
4. **Untrusted/trusted split (R4 pattern)**: tiny SDP (5 barrier scalars
   + 2 S-procedure multipliers) synthesized by cvxpy/Clarabel with slack
   margins; solution rationalized with bounded denominators (margins
   absorb rounding); checker rebuilds every matrix from the certificate
   scalars and integer bases and decides PSD by all principal minors
   (127 exact determinants at 7x7).

## Receipts (all green)

- checker's integer matrix == cw.cw_deriv in scaled coordinates;
- end-to-end: synthesize -> rationalize -> verify_certificate returns
  no violations; certificate carries only Fractions;
- hand-derived certificate B = -c1^2 + 2Ax^2 + 2Az^2 + 2R'^2 with
  lam_u = 2: the identity B - 2 g_koz = (X-u2)^2 + Z^2 + ... closes on
  paper AND under the exact checker;
- tamper detection: sign-flipped coefficient and negative multiplier
  both rejected;
- V-bar-centered hold: synthesis infeasible — correctly, since V-bar
  coast passes through the target (heritage agrees);
- corroboration: X0 corner states through two orbits of dense exact CW
  flow never enter the keep-out radius.

## Deferred

Quartic barriers (tighter X0 geometries), J2-perturbed drift variant,
attitude-loop barriers, emitting the certificate as an audit artifact
alongside the reach gate.

## Push/merge instructions

Single commit on main: `20 — SOS abort-safety certificates (#20)`;
push; close.
