# Trajectory optimization: G-FOLD lineage and the convex RPOD stack

How rpod-lib incorporates convexification, from the G-FOLD heritage to a
layered RPOD formulation. References collected at the end of each section.

## What G-FOLD contributes, and what transfers

G-FOLD (Açıkmeşe & Ploen 2007; Blackmore, Açıkmeşe & Scharf 2010) solved
min-fuel powered descent by **lossless convexification (LCvx)**: the nonconvex
thrust annulus `0 < ρ₁ ≤ ‖T‖ ≤ ρ₂` is relaxed with a slack Γ
(`‖T‖ ≤ Γ`, `ρ₁ ≤ Γ ≤ ρ₂`, minimize ∫Γ), and a Pontryagin/maximum-principle
argument shows the relaxation is exact — the SOCP returns the global optimum
of the nonconvex problem. A change of variables (`u = T/m`, `z = ln m`)
linearizes the mass depletion. Flight-proven onboard (Masten Xombie, 2012-13).

**Transfers cleanly to RPOD** — in fact more cleanly, since CW/TH dynamics are
already linear time-varying:

- Thrust-annulus LCvx for LTV systems (Kunhippurayil et al., Automatica 2021).
- **Lu & Liu (JGCD 2013)** — the key RPOD-LCvx result: exactness of the
  relaxation for rendezvous with a target in an arbitrary orbit, with active
  approach-corridor / plume / velocity constraints, solved as a sequence of
  SOCPs.
- Semi-continuous thrust (`‖u‖ ∈ {0} ∪ [ρ₁, ρ₂]`, i.e. RCS on/off with a
  minimum throttle) convexified without integer variables (Malyuta &
  Açıkmeşe 2020).

**Does not transfer:**

- **Keep-out zones are not losslessly convexifiable** — an open problem
  (Malyuta et al. survey, 2022). Handled by rotating hyperplanes, SCP
  linearization, or (for benchmarking) mixed-integer formulations.
- LCvx guarantees weaken when convex state constraints are active over
  *intervals* (a corridor ridden through final approach) rather than isolated
  instants, and are continuous-time statements — discrete-time repair exists
  (Luo et al. 2024) and validity checks (controllability/normality) must be
  part of the library, not an exercise for the user.

## Layered formulation in rpod-lib

**Layer 0 — one-shot convex (LP/QP/SOCP).** CW or Yamanaka-Ankersen STMs give
*exact* discretization, so dynamics are linear equality constraints with zero
integration error. Costs: axis-wise L1 fuel → LP; `Σ‖Δv‖₂` → SOCP; quadratic
tracking → QP. Convex constraint catalogue:

| Constraint | Form |
|---|---|
| Approach cone (apex, axis d̂, half-angle θ) | `‖(I − d̂d̂ᵀ)(r − r_p)‖ ≤ tanθ · d̂ᵀ(r − r_p)` — exact SOC |
| Keep-out zone | rotating hyperplane `n̂ₖᵀ(rₖ − r_t) ≥ R` (Mueller & Larsson 2008) — linear, conservative |
| Plume inhibition | thrust half-space `n̂ᵀuₖ ≤ 0` near dock (Weiss et al. 2015) |
| Passive safety (Breger & How 2008) | free-drift states after failure at tⱼ are linear in the decision variables ⇒ KOZ avoidance over the safety horizon as linear constraints per failure scenario |
| Thrust annulus | LCvx relaxation with validity checks (Lu & Liu style) |

Layer 0 is globally optimal, needs no initial guess, and covers most far/mid-
field planning. This is the v0.1 target.

**Layer 1 — successive convexification for near-field and 6-DOF docking.**
For coupled translation+attitude, exact KOZ shaping, min-impulse-bit and
range-triggered (state-triggered) constraints:

- **PTR (penalized trust region) as the default engine** — simplest, best
  empirical convergence, flight pedigree via NASA SPLICE's dual-quaternion
  guidance (ran real-time onboard Blue Origin NS-13/NS-17, listen-only).
- **SCvx\* as the guaranteed option** (Oguri 2023) — an augmented-Lagrangian
  outer loop over the same subproblem, with convergence *to a feasible* local
  optimum (plain SCvx can converge with nonzero virtual control).
- Exact FOH discretization via STM integration; virtual control; continuous-
  time constraint satisfaction (CTCS) augmentation for inter-sample cone/KOZ
  violations; time dilation for free final time; state-triggered constraints
  for plume/MIB logic (Malyuta et al. 2019 solved Apollo transposition &
  docking this way — no integers, interactive speeds).

GuSTO is cited but not implemented (control-affine restriction, deprecated
reference code).

**Layer 2 — safety modules.** Breger-How passive-safety constraints in Layers
0/1; abort-safe backward-reachable-set half-spaces (Aguilar Marsillach, Di
Cairano & Weiss 2020/22) as an optional module.

## Solver stack (license-vetted, permissive end-to-end)

| Stage | Choice | License |
|---|---|---|
| Modeling front end | CVXPY, **DPP-parametrized from day one** (makes codegen free later) | Apache-2.0 |
| Prototyping solver | Clarabel (default); QOCO, SCS alternates | Apache-2.0 / BSD-3 / MIT |
| Embedded SOCP | **CVXPYgen → QOCOGEN**: library-free, static-memory custom C, purpose-built for trajectory SOCPs | Apache-2.0 / BSD-3 |
| Embedded QP (tracking MPC) | OSQP v1.0 codegen: malloc-free, division-free, warm-startable | Apache-2.0 |
| SCP outer loop | In-house, dependency-free C (~hundreds of lines): linearize → generated solver → penalty update, fixed iteration caps, deterministic fallback trajectory | ours |

**GPL islands kept out of the dependency tree:** ECOS, CVXOPT, SCPToolbox.jl,
jonnyhyman/G-FOLD. SCPToolbox.jl and the CSM 2022 survey are the *algorithm
references* (read, don't port). The ART repo (Stanford, MIT license) is
legally liftable for RPOD SCP pipeline patterns. **OpenSCvx** (Apache-2.0,
active, JAX-based) is the closest existing SCP core — we interoperate with or
depend on it rather than duplicating; rpod-lib's differentiated value is the
**RPOD domain layer** (CW/YA dynamics, corridors, KOZ, plume, passive/abort
safety, docking phase logic), which nothing open-source currently packages.

The embedded path aligns with the verification story (docs/verification.md):
generated solvers are static-memory, fixed-iteration, `math.h`-only C — the
same dialect the external abstract-interpretation tool is designed to prove.
Real-time reference points: customized IPM SOCP for 3-DOF descent ran in
~0.7 s on a RAD750; SPLICE replaced its IPM with a customized first-order
PIPG solver (matrix-free, auto-coded C) to hit budget — the same escalation
path is open to us.

## Honest caveats

- LCvx guarantees can silently break under interval-active state constraints
  and after discretization — rpod-lib ships validity checks and treats the
  discrete-time repair as part of the implementation.
- PTR has no formal convergence guarantee; SCvx* mode exists for when that
  matters.
- Rotating-hyperplane KOZ is conservative and rate-tuned; SCP-linearized KOZ
  and an optional MICP benchmark keep it honest.
- No public evidence confirms SCP running closed-loop in operational orbital
  RPOD; the certification argument rests on descent heritage (G-FOLD, SPLICE)
  plus fixed-iteration determinism. Visiting-vehicle RPOD flown to date is
  classically glideslope/corridor-based — which is why rpod-lib implements
  those laws first and treats optimization as the upgrade path.

## Key references

Açıkmeşe & Ploen, JGCD 2007 (G-FOLD SOCP) · Blackmore, Açıkmeşe & Scharf,
JGCD 2010 (min-landing-error) · Açıkmeşe & Blackmore, Automatica 2011 (LTV
LCvx) · Lu & Liu, JGCD 2013 (rendezvous LCvx) · Mueller & Larsson 2008
(rotating hyperplane) · Richards, Schouwenaars, How & Feron, JGCD 2002 (MILP
avoidance) · Breger & How, JGCD 2008 (safe trajectories) · Weiss, Baldwin,
Erwin & Kolmanovsky, IEEE TCST 2015 (CW MPC with LOS/plume/soft-dock) ·
Malyuta et al., IEEE CSM 2022, arXiv:2106.09125 (SCP survey) · Mao, Szmuk &
Açıkmeşe, arXiv:1608.05133 (SCvx) · Oguri, arXiv:2304.14564 (SCvx*) · Szmuk &
Açıkmeşe, arXiv:1802.03827 (PTR) · Malyuta et al., arXiv:1906.04857 (docking
via STCs) · Elango et al., arXiv:2404.16826 (CTCS) · Luo, Echigo & Açıkmeşe,
arXiv:2410.09748 (discrete-time LCvx repair) · Dueri et al., JGCD 2017
(RAD750 timing) · NASA SPLICE DQG, NTRS 20240014010.
