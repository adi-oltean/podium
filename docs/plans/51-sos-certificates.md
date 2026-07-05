# 51 — Higher-degree exact SOS/Positivstellensatz certificates

GitHub issue: https://github.com/adi-oltean/podium/issues/51

## Motivation

The abort-safety barrier (#20) is the quadratic S-procedure: the
degree-two special case of the sum-of-squares / Positivstellensatz
barrier-certificate approach. Reviewer feedback on the paper's SOS
discussion asked whether the SOS work is innovative; the honest answer
was that the quadratic case is standard. This issue takes the
defensible step: a general-degree exact SOS certificate checker,
demonstrated on a genuinely nonlinear (cubic) system with a quartic
barrier.

## Fix (landed) — `podium.verify.sos`

Polynomials as `dict[exponent-tuple, Fraction]` with exact operations:
`padd`, `psub`, `pscale`, `pmul`, `pdiff` (partial derivative), and
`lie_derivative(V, f)` = sum_i (dV/dx_i) f_i along a polynomial vector
field. `is_sos(p, basis, gram)` certifies p is a sum of squares by
checking the EXACT polynomial identity p = z^T G z (coefficient by
coefficient) AND G >= 0 (the all-principal-minors PSD test reused from
`verify.barrier`). No floating point in the trusted path: a float SOS
solver may synthesize (G, z), but the shipped certificate is the
rational Gram matrix re-verified here.

## Receipts (tests/test_sos.py)

- Diagonal Gram: x1^4 + x2^4 is SOS. Non-diagonal rank-1 Gram:
  (x1^2 + x2^2)^2 is SOS. Both by exact identity + PSD.
- Indefinite x1^2 - x2^2 rejected two ways (indefinite Gram fails PSD;
  a wrong Gram fails the identity).
- Duffing oscillator xdot1 = x2, xdot2 = -x1 - x1^3 - x2 with the
  quartic energy V = 1/2 x1^2 + 1/4 x1^4 + 1/2 x2^2: `lie_derivative`
  yields dV/dt = -x2^2 EXACTLY (the cubic cross terms cancel), and
  -dV/dt = x2^2 is certified SOS -> every sub-level set {V <= c} is an
  infinite-horizon invariant of a NONLINEAR system, with a quartic
  certificate, verified exactly.
- Simulation confirms V is monotone non-increasing along the flow.

## Significance

This is the higher-degree, exact-Positivstellensatz barrier the paper
flagged as the defensible novelty (versus merely integrating the
standard quadratic case). #20's abort-safety barrier is now the
degree-2 special case of this machinery.

## Deferred

Synthesizing higher-degree barriers from an untrusted SDP + a validated
rational rounding (so the exact identity holds with a margin) for RPOD
systems that are not exactly-rational by construction; SOS-multiplier
Positivstellensatz for set-separation (init/unsafe) beyond the Lie
condition; applying it to the polynomial quaternion-feedback closed
loop.

## Push/merge instructions

Single commit on main: `51 — Higher-degree exact SOS certificates
(#51)`; push; close.
