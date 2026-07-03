# Verified onboard optimization: research synthesis and Podium ideas

Research date: 2026-07-03. Sources: published work of Behçet Açıkmeşe (UW
ACL) and Pierre-Loïc Garoche (ENAC), plus the certified-optimization
literature connecting them. Verification flags: [V] = source fetched and
read; [S] = corroborated search metadata; [gap] = verified *absence* after
targeted search.

## The strategic picture

Two schools, one empty middle ground:

- **Açıkmeşe's school builds solvers simple enough to verify but never
  proves them.** Customized IPMs (Dueri, JGCD 2017,
  [doi:10.2514/1.G001480](https://arc.aiaa.org/doi/10.2514/1.G001480)) and
  the PIPG line (Automatica 2022,
  [arXiv:2108.10260](https://arxiv.org/abs/2108.10260); xPIPG
  [arXiv:2203.04188](https://arxiv.org/abs/2203.04188); infeasibility
  detection [arXiv:2109.02756](https://arxiv.org/abs/2109.02756)) are
  factorization-free, static-memory, division-light — culminating in
  DQG+PIPG meeting NASA update-rate requirements on SPLICE hardware
  ([arXiv:2508.10439](https://arxiv.org/abs/2508.10439)). Their doctrine is
  verifiability *by construction* plus HWIL demonstration; [gap] no
  per-problem iteration certification for PIPG exists (rates are
  asymptotic), and no ACL paper does code-level formal proof.
- **Garoche's school proves solvers, but only toy/LP-class ones.** Credible
  autocoding (Wang et al., [arXiv:1307.2641](https://arxiv.org/abs/1307.2641);
  convex-optimization autocoding [arXiv:1403.1861](https://arxiv.org/abs/1403.1861);
  Princeton UP book 2019): Lyapunov/ellipsoid invariants as ACSL contracts
  discharged by Frama-C/WP. Cohen/Feron/Garoche verified the *ellipsoid
  method* precisely because it has an a-priori iteration bound
  ([arXiv:2005.12588](https://arxiv.org/abs/2005.12588)); Davy et al. did
  the first code-level IPM proof ([arXiv:1801.03833](https://arxiv.org/abs/1801.03833)).
  All LP-class, real-arithmetic, never flown; the Simulink toolchain
  (CoCoSim/Gene-Auto) is dormant.
- **Convergence is happening now:** Garoche is autocoding Açıkmeşe's
  algorithm class (SCvxPyGen, IEEE CDC 2024; verifiable LP-based SCvx,
  IAGNC 2025, [hal-05175042](https://hal.science/hal-05175042)) with an ESA
  verified-optimization project running 2026–29. Stellato's performance
  verification bounds worst-case residuals after N iterations over
  parametric QP families ([arXiv:2403.03331](https://arxiv.org/abs/2403.03331)).
  Arnström/Axehill certify active-set QP complexity to exact hardware WCET
  ([arXiv:2304.11576](https://arxiv.org/abs/2304.11576)); [gap] nothing
  equivalent for SOCP/ADMM/PIPG or any spacecraft problem.
- **The agencies confirm the vacuum**: NASA/ESA "Call to Action for
  Advanced GNC Algorithm V&V," 2024
  ([NTRS 20240003178](https://ntrs.nasa.gov/citations/20240003178)) — the
  certification methodology for optimization-based GNC "is not currently
  established."

A living open-source library whose design layer emits both the flight code
and its certificates does not exist. Podium's contracts-as-data + static
subset + planned codegen is architecturally the right vehicle.

## Eight ideas, ranked by impact x feasibility

| # | Idea | Impact | Feasibility | First-in-open-source? |
|---|------|--------|-------------|----------------------|
| R1 | **CI reachability regression** on the ARCH rendezvous benchmark: export Podium's FSM+gains as the hybrid model, JuliaReach/CORA re-proves LOS-cone/velocity/passive-abort specs on every guidance/control PR | High | High | Yes ([gap]: no GNC library runs reachability as a CI gate; ARCH-COMP is annual+manual) |
| R2 | **Certificate-carrying contracts**: `EllipsoidInvariant(P, level)` beside `Interval` — synthesis returns the Lyapunov certificate (FP-inflated via validated SDP), sandbox checks x'Px preservation each step, C emitter renders it as a quadratic PROVE/ACSL obligation | High | High | Yes, as a library feature (Feron/Wang chain is dead; Khalife/Garoche is hand-written C) |
| R3 | **@contract -> ACSL backend** in the C emitter: generated flight code checkable by anyone with Frama-C/WP; E-ACSL gives runtime checks in golden-vector CI for free | High | Very high | Yes ([gap]: no Python-contracts->ACSL pipeline exists; CVXPYgen/OSQP codegen emit nothing) |
| R4 | **Static-subset xPIPG + verified KKT checker**: the solver is ideal static-subset material (matrix-vector only, fixed MAX_ITER, anytime); the *formally verified artifact is the ~50-line a-posteriori certificate checker*, not the solver — proof-carrying answers | Very high | Medium | Yes ([gap]: no code-level PIPG verification, no published in-flight certificate checking anywhere) |
| R5 | **Certified MIB pulse allocator**: allocation-error-per-cycle bounded by the quantization cell as a prove() contract, plus an LMI ultimate-bound ellipsoid on the dither/carry closed loop ("corridor holds despite pulse quantization") | High (domain-unique) | Medium | Yes — zero formal-verification prior art for MIB allocation anywhere |
| R6 | **Contract-parametrized iteration certification**: the contract box *is* the parametric problem family; certify "N iterations => KKT residual <= eps" over it (Ranjan-Stellato-style), regenerate in CI when contracts change; MAX_ITER becomes a derived, certified quantity | Very high (publishable) | Medium-low | Yes for SOCP/PIPG/spacecraft families |
| R7 | **FP-bound-derived CI tolerances**: FPTaylor/Daisy bounds on straight-line kernels stored as contract metadata; Python<->C golden-vector gate asserts the analyzer-derived tolerance | Medium | High | Effectively yes |
| R8 | **Proof-carrying SCP**: per-iteration acceptance certificates + a contractual guarantee of "certificate-backed trajectory OR pre-verified fallback whose safety R1 re-proves" | Medium-high | Medium | Contested — track Garoche's SCvxPyGen/IAGNC 2025 closely |

Sequencing against the roadmap: R7+R3 as groundwork (cheap; v0.2–v0.5
emitter), R1 at v0.2 (benchmark export already planned — the novelty is
gating on it), R2 with the v0.3 control work, R5 at v0.3 actuators, R4 at
the v0.5 solver stage, R6/R8 as the research-grade layer. The composite
endgame — one release where CI simultaneously re-proves reachability (R1),
re-derives certified MAX_ITER (R6), and re-discharges ACSL obligations
(R3) — is precisely the artifact the NASA/ESA Call to Action says doesn't
exist.

## Caveats

A few AIAA/IEEE landing pages 403'd (flagged [S]); Dueri's iteration-bound
specifics are attributed rather than fetched; read Garoche's IAGNC 2025
verifiable-SCvx paper in full before committing to R8's novelty claim.
