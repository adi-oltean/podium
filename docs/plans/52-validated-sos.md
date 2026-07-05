# 52 — Validated SOS: untrusted float SDP to exact rational certificate

GitHub issue: https://github.com/adi-oltean/podium/issues/52

## Motivation

#51 gave an exact SOS/Positivstellensatz CHECKER of arbitrary degree,
but the demonstrations used exactly-rational certificates (hand-built
or from an exactly-representable identity). Real SOS SYNTHESIS uses a
floating-point SDP, whose solution is not exactly rational. This issue
adds the missing link: turn a float SDP Gram into an exact rational
certificate (the validated-SOS tradition), which is the piece the
paper's higher-degree-barrier future-work flagged as open.

## Fix (landed) — `sos.validate_gram`

`validate_gram(target, basis, gram_float, margin)` round-and-corrects:

1. Rationalize the float Gram (limit_denominator).
2. Symmetrize exactly and add a small rational `margin` to the diagonal
   — the interior-point SDP Gram is strictly positive definite, so this
   preserves PSD with slack after the tiny corrections.
3. Absorb the exact coefficient residual `target - z^T G z` monomial by
   monomial. The correction DECOUPLES: each Gram entry contributes to
   exactly one product monomial, so for each residual monomial pick one
   entry with that product (diagonal preferred) and adjust it by
   residual/weight (weight 1 diagonal, 2 off-diagonal). No general
   linear solve is needed.

Returns the exact Gram, or None if a target monomial has no basis pair
producing it (the basis is too small for an SOS form).

## Receipt (tests/test_sos.py)

An untrusted Clarabel SDP synthesizes an SOS Gram for a quartic that
REQUIRES off-diagonal entries, q = (x1^2 + x1 x2 + x2^2)^2. A naive
rationalization of the float Gram leaves a nonzero residual (so
validation is genuinely necessary, not cosmetic); `validate_gram`
produces an exact Gram that reproduces q identically and passes the
exact `is_sos` check (#51), with the off-diagonal freedom genuinely
used. A too-small basis is rejected.

## Significance

This completes the validated-SOS pipeline: an untrusted float SDP
synthesizes, the shipped certificate is exact rational, and it is
re-verified with no floating point in the trusted path — the same
"synthesize with an untrusted solver, check the answer in exact
arithmetic" discipline as the barrier (#20), KKT (#40/#41), and
Lyapunov (#50) certificates, now for arbitrary-degree SOS.

## Deferred

Applying it to a synthesized higher-degree barrier for the polynomial
quaternion-feedback closed loop (needs the barrier-synthesis SDP set up
for that system); a rigorous a-priori margin bound (currently the
margin is a fixed small rational validated by the exact PSD check
post hoc).

## Push/merge instructions

Single commit on main: `52 — Validated SOS synthesis (#52)`; push;
close.
