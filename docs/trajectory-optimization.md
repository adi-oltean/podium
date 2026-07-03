# Trajectory optimization: G-FOLD lineage and the convex RPOD stack

How Podium incorporates convexification, from the G-FOLD heritage to a
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

## Layered formulation in Podium

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
depend on it rather than duplicating; Podium's differentiated value is the
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
  and after discretization — Podium ships validity checks and treats the
  discrete-time repair as part of the implementation.
- PTR has no formal convergence guarantee; SCvx* mode exists for when that
  matters.
- Rotating-hyperplane KOZ is conservative and rate-tuned; SCP-linearized KOZ
  and an optional MICP benchmark keep it honest.
- No public evidence confirms SCP running closed-loop in operational orbital
  RPOD; the certification argument rests on descent heritage (G-FOLD, SPLICE)
  plus fixed-iteration determinism. Visiting-vehicle RPOD flown to date is
  classically glideslope/corridor-based — which is why Podium implements
  those laws first and treats optimization as the upgrade path.

## Key references

**Lossless convexification (LCvx) and powered-descent heritage**

1. B. Açıkmeşe and S. R. Ploen, "Convex Programming Approach to Powered
   Descent Guidance for Mars Landing," *Journal of Guidance, Control, and
   Dynamics* 30(5):1353–1366, 2007.
   [doi:10.2514/1.27553](https://doi.org/10.2514/1.27553) — the G-FOLD SOCP.
2. L. Blackmore, B. Açıkmeşe, and D. P. Scharf, "Minimum-Landing-Error
   Powered-Descent Guidance for Mars Landing Using Convex Optimization,"
   *JGCD* 33(4):1161–1171, 2010.
   [doi:10.2514/1.47202](https://doi.org/10.2514/1.47202)
3. B. Açıkmeşe and L. Blackmore, "Lossless Convexification of a Class of
   Optimal Control Problems with Non-Convex Control Constraints,"
   *Automatica* 47(2):341–347, 2011.
   [doi:10.1016/j.automatica.2010.10.037](https://doi.org/10.1016/j.automatica.2010.10.037)
   — the general LTV LCvx theorem.
4. P. Lu and X. Liu, "Autonomous Trajectory Planning for Rendezvous and
   Proximity Operations by Conic Optimization," *JGCD* 36(2):375–389, 2013.
   [doi:10.2514/1.58436](https://doi.org/10.2514/1.58436) — LCvx exactness
   for rendezvous with active state constraints.
5. J. Kunhippurayil, M. W. Harris, and O. Jansson, "Lossless Convexification
   of Optimal Control Problems with Annular Control Constraints,"
   *Automatica* 133:109848, 2021.
   [doi:10.1016/j.automatica.2021.109848](https://doi.org/10.1016/j.automatica.2021.109848)
6. K. Luo, T. Echigo, and B. Açıkmeşe, "Lossless Convexification in
   Discrete-Time Optimal Control," 2024.
   [arXiv:2410.09748](https://arxiv.org/abs/2410.09748) — discrete-time
   violation bound and repair.
7. D. Malyuta and B. Açıkmeşe, "Lossless Convexification of Optimal Control
   Problems with Semi-Continuous Inputs," IFAC World Congress, 2020.
   [arXiv:1911.09013](https://arxiv.org/abs/1911.09013) — RCS on/off with
   minimum throttle, no integers.

**Convex RPOD formulations and constraints**

8. M. Tillerson, G. Inalhan, and J. P. How, "Co-ordination and Control of
   Distributed Spacecraft Systems Using Convex Optimization Techniques,"
   *Int. J. Robust and Nonlinear Control* 12(2–3):207–242, 2002.
   [doi:10.1002/rnc.683](https://doi.org/10.1002/rnc.683) — canonical LP
   transcription on relative-motion STMs.
9. E. Mueller and R. Larsson, "Collision Avoidance Maneuver Planning with
   Robust Optimization," ESA GNC Conference, 2008 (also AIAA 2009-2051,
   [doi:10.2514/6.2009-2051](https://doi.org/10.2514/6.2009-2051)) —
   rotating-hyperplane keep-out-zone constraint.
10. A. Richards, T. Schouwenaars, J. P. How, and E. Feron, "Spacecraft
    Trajectory Planning with Avoidance Constraints Using Mixed-Integer
    Linear Programming," *JGCD* 25(4):755–764, 2002.
    [doi:10.2514/2.4943](https://doi.org/10.2514/2.4943)
11. L. S. Breger and J. P. How, "Safe Trajectories for Autonomous
    Rendezvous of Spacecraft," *JGCD* 31(5):1478–1489, 2008.
    [doi:10.2514/1.29590](https://doi.org/10.2514/1.29590) — passive-safety
    constraints, linear per failure scenario.
12. A. Weiss, M. Baldwin, R. S. Erwin, and I. Kolmanovsky, "Model
    Predictive Control for Spacecraft Rendezvous and Docking: Strategies
    for Handling Constraints and Case Studies," *IEEE Trans. Control
    Systems Technology* 23(4):1638–1647, 2015.
    [doi:10.1109/TCST.2014.2379639](https://doi.org/10.1109/TCST.2014.2379639)
    — CW MPC with LOS cone, soft-dock velocity bounds, plume direction.
13. D. Aguilar Marsillach, S. Di Cairano, and A. Weiss, "Abort-Safe
    Spacecraft Rendezvous in Case of Partial Thrust Failure," IEEE CDC 2020
    / *IEEE TCST* 2022.
    [MERL TR2022-142](https://www.merl.com/publications/docs/TR2022-142.pdf)
    — backward-reachable-set abort safety.

**Successive convexification**

14. D. Malyuta, T. P. Reynolds, M. Szmuk, T. Lew, R. Bonalli, M. Pavone,
    and B. Açıkmeşe, "Convex Optimization for Trajectory Generation," *IEEE
    Control Systems Magazine* 42(5):40–113, 2022.
    [arXiv:2106.09125](https://arxiv.org/abs/2106.09125) — the survey; also
    states KOZ non-convexifiability as an open problem.
15. Y. Mao, M. Szmuk, and B. Açıkmeşe, "Successive Convexification of
    Non-Convex Optimal Control Problems and Its Convergence Properties,"
    IEEE CDC 2016. [arXiv:1608.05133](https://arxiv.org/abs/1608.05133)
16. K. Oguri, "Successive Convexification with Feasibility Guarantee via
    Augmented Lagrangian for Non-Convex Optimal Control Problems," IEEE
    CDC 2023. [arXiv:2304.14564](https://arxiv.org/abs/2304.14564) — SCvx*.
17. M. Szmuk and B. Açıkmeşe, "Successive Convexification for 6-DoF Mars
    Rocket Powered Landing with Free-Final-Time," AIAA SciTech 2018.
    [arXiv:1802.03827](https://arxiv.org/abs/1802.03827) — PTR.
18. D. Malyuta, T. Reynolds, M. Szmuk, B. Açıkmeşe, and M. Mesbahi, "Fast
    Trajectory Optimization via Successive Convexification for Spacecraft
    Rendezvous with Integer Constraints," AIAA SciTech 2020.
    [arXiv:1906.04857](https://arxiv.org/abs/1906.04857) — Apollo
    transposition & docking via state-triggered constraints.
19. P. Elango, D. Luo, A. G. Kamath, S. Uzun, T. Kim, and B. Açıkmeşe,
    "Successive Convexification for Trajectory Optimization with
    Continuous-Time Constraint Satisfaction," 2024.
    [arXiv:2404.16826](https://arxiv.org/abs/2404.16826) — CTCS.

**Real-time / onboard**

20. D. Dueri, B. Açıkmeşe, D. P. Scharf, and M. W. Harris, "Customized
    Real-Time Interior-Point Methods for Onboard Powered-Descent Guidance,"
    *JGCD* 40(2):197–212, 2017.
    [doi:10.2514/1.G001480](https://doi.org/10.2514/1.G001480) — ~0.7 s
    SOCP solves on RAD750.
21. NASA SPLICE Dual-Quaternion Guidance flight results, AIAA SciTech 2025.
    [NTRS 20240014010](https://ntrs.nasa.gov/citations/20240014010) —
    onboard SCvx-class guidance (listen-only) on Blue Origin NS-13/NS-17;
    IPM-to-PIPG solver migration.
