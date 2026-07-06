# Certified successive convexification

Successive convexification (SCvx/PTR) rounds a nonconvex keep-out
`||r - c|| >= R` by linearizing it into a half-space cut at each
reference. With a slack the cut is soft, so intermediate iterates can
violate the true constraint; only the converged trajectory is (nearly)
feasible. Applied **hard** with a **rational-unit** normal, the cut is a
*sound convex inner-approximation* — every point of the half-space
genuinely lies outside the ball — so **every node of every iterate is
feasible for the true nonconvex problem**, not just the limit.

`podium.verify.scvx_cut` makes that soundness a machine-checked,
exact-rational certificate, in the same discipline as the barrier / KKT
/ Lyapunov / SOS certificates.

## The certificate

For a rational unit normal `n`, `u = r - c`, `h = n·u - R`, the module
verifies the Positivstellensatz identity

```
||u||^2 - R^2  =  ( (n_perp·u)^2 + (n·u - R)^2 )  +  2R (n·u - R)
                   \-------- SOS  s0(r) --------/     \-- s1·h --/
```

with SOS block `s0` (carried by an exact rational Gram) and multiplier
`s1 = 2R >= 0`. On the cut `{h >= 0}` the right side is `SOS + nonneg
>= 0`, so `||u|| >= R`. The Gram can be re-synthesized from an untrusted
float SDP and validated exactly with `sos.validate_gram` — the same
float-synthesis → exact-certificate mechanism as the SOS work.

`snap_rational_unit` returns the nearest rational unit vector to a
reference direction exactly (circle parametrized by `t = tan(θ/2) ∈ Q`),
keeping the whole certificate in the rationals.

## What is demonstrated (`test_certified_scvx.py`)

A 2-D min-energy transfer rounds a spherical keep-out with hard
certified cuts. The receipts:

- **Per-iterate feasibility.** Every node of every SCvx iterate
  satisfies `||r|| >= R` — the true nonconvex constraint — from the
  first iterate on.
- **Sound cuts.** Every cut certifies exactly (SOS block + S-procedure
  identity).
- **Validated synthesis.** One cut's SOS Gram, synthesized by an
  untrusted float SDP, is validated to an exact rational certificate.

## Scope

The guarantee is at the trajectory nodes; continuous-time (between-node)
keep-out uses the dense CTCS cuts in `podium.guidance.scp`. Extending
the certified sound cut to higher-degree (ellipsoidal / superquadric)
keep-outs is the natural next step — the exact SOS machinery already
supports it.
