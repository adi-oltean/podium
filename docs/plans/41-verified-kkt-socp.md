# 41 — Verified-KKT checker for SOCP (exact cone complementarity)

GitHub issue: https://github.com/adi-oltean/podium/issues/41

## Fix (landed) — `verify_socp` in `podium.verify.kkt`

Extends the exact-rational optimality checker (#40) from QP to
second-order-cone programs — the form the conic Layer-0 guidance
actually solves (keep-out cones, thrust cones). For a standard-form
SOCP

    minimize   c' x
    subject to A x = b,  G x + s = h,  s in K

with K = R+^l x (product of second-order cones), `verify_socp` returns
a `SOCPReport` of EXACT residuals in `fractions.Fraction`:

- stationarity   max|c + A' y + G' z|
- eq_residual    max|A x - b|
- conic_residual max|G x + s - h|
- comp_slack     |s' z|
- s_in_cone / z_in_cone : exact cone membership of the slack and dual

Both K and its dual are self-dual (SOC and nonneg orthant), so z is
checked against the same K.

## The sqrt-free cone test

Second-order-cone membership s0 >= ||s1:|| involves a square root,
which has no exact rational value. `_soc_margin_ok` avoids it: for a
block (s0, s1:), approximate membership to tolerance tol is
`s0 >= -tol AND ||s1:||^2 <= (s0+tol)^2` — pure Fraction comparisons,
so the verdict has no floating-point uncertainty even for the cone.

## Receipts (tests/test_kkt.py)

- `_soc_margin_ok` decided exactly on interior / boundary / exterior /
  tolerance-boundary points (the exterior (1,1,1) is pulled into the
  cone exactly when tol crosses sqrt2 - 1).
- LP-as-SOCP (pure nonneg cone) with a rational optimum: every conic
  residual EXACTLY Fraction(0), certified at tol 0.
- Untrusted ECOS solve (native standard-form x/y/z/s duals) of a
  min-norm-with-thrust-cone SOCP: certifies within 1e-5 with the cone
  active (t = ||u||), objective matches.
- Negating the cone dual (out of K*) or moving the primal off the cone
  fails membership exactly and is rejected.

## Deferred

Feeding an embedded CVXPYgen/QOCOGEN KKT dump straight to the checker;
exponential/PSD cones; an interval variant for problems without a
clean rational optimum.

## Push/merge instructions

Single commit on main: `41 — Verified-KKT checker for SOCP (#41)`;
push; close.
