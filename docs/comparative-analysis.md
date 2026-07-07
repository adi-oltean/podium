# Comparative analysis of existing open-source frameworks

Survey date: July 2026 (updated 2026-07-03). Full framework-by-framework
details in the table below; conclusions first.

Update 2026-07-03 (verified against repos/registries):

- **Orekit MR !1133 is still open** (milestone 14.0, under active review):
  a `propagation.relative` package (CW + YA propagators, two-impulse
  rendezvous) plus relative maneuvers. When merged, Orekit gains its first
  native RPOD package — Java-side, CW/YA only, no ROE, no optimization.
- **Basilisk latest is v2.10.2** (2026-05-08); RPOD-relevant additions:
  `hillFrameRelativeControl`, differential-drag formation scenario,
  `stochasticAtmDensity`, docked-two-craft constraint effector (2.9).
- **No open implementation of relative-orbital-elements (ROE) dynamics**
  (Koenig/D'Amico STMs, e/i-vector passive safety) exists in any surveyed
  tool — see the roadmap's ROE module (issue #5).
- SCPToolbox.jl remains dormant; **OpenSCvx** (Python/JAX, CTCS built-in)
  is the active SCP core. **Brahe** (Rust/Python, MIT, JOSS) is a new
  pedagogy-first entrant. `dyn4space` does not exist; SaRA is
  human-robot-safety software, not RPOD (both occasionally mis-cited).
- AFRL `safe-autonomy-sims` remains the RL docking reference and models
  continuous thrust only — no minimum-impulse-bit actuation.

## The gap Podium fills

No maintained open-source tool combines all of:

1. **validated relative-motion dynamics** (CW + Tschauner-Hempel/Yamanaka-
   Ankersen + nonlinear truth cross-checks) — *no surveyed tool ships an open
   Yamanaka-Ankersen implementation at all*;
2. **docking contact/capture mechanics** coupled to those dynamics;
3. **sensor/actuator models** relevant to RPOD (relative GNSS, docking
   camera, lidar, thruster MIB);
4. **pythonic GNC algorithm hooks** (the layer between sim and flight code
   that vendors sell closed);
5. **Monte Carlo** as a first-class citizen;
6. **a deliberate path to flight C code** with formal validation.

State of practice remains STK/FreeFlyer/MATLAB-Simulink cobbled together, or
Basilisk plus custom code. The open niches Podium targets: RPOD guidance
primitives as library objects (glideslope, hops, safety ellipses, KOZ/passive-
abort, plume constraints), analytic relative-motion rigor, and the
verification-ready codegen path — under a permissive license (the incumbents
are fractured across ISC/NOSA/AGPL/"government-owned").

## Decisions: build on / interop / inspiration

### Build on (dependencies, all optional extras)

| Dependency | License | Role |
|---|---|---|
| **Basilisk** (`pip install bsk`) | ISC | Optional high-fidelity dynamics backend: coupled 6-DOF, rich sensor/actuator models, Hill-frame modules, mature Monte Carlo. Kept a *backend*, never the core — our Python→C story must not depend on its C++ internals |
| **MuJoCo/MJX** | Apache-2.0 | Optional contact/docking-capture backend, already wrapped by `podium.sim.contact`: Podium owns the RPOD-specific probe-drogue geometry, capture criteria, IDSS tie-in, and tests, while MuJoCo supplies the contact solve. The field converged on it (Basilisk ≥2.8, SmallSatSim) |
| **ANISE** | MPL-2.0 | Ephemerides/frames/time (SPICE-compatible, flight-proven, pip wheels). spiceypy (MIT) as an alternate adapter behind the same abstraction |
| **CVXPYgen / TinyMPC** | MIT | The Python→embedded-C path for optimization-based guidance |
| **heyoka.py** | MPL-2.0 | Optional truth-propagation/variational backend (free STMs, event detection) |

### Interop (adapters and CI oracles, never dependencies)

- **Orekit** (orekit-jpype) and **GMAT** — translational-propagation
  cross-validation oracles in CI; **tudatpy** for torque-coupled 6-DOF.
  Interchange via **CCSDS OEM/OPM**.
- **cFS** and **F´** — downstream flight-software consumers of the generated
  C GNC modules (software-bus/UDP bridge, NOS3-style). **NASA 42** as an
  optional socket-IPC cross-check.
- **Gymnasium env conventions** (compat shim for the RL community),
  **SPEED+/SPNv2 datasets** for vision-in-the-loop navigation validation.

### Design inspiration only

- **poliastro** (MIT, archived) — the strongest architectural donor: a
  units-safe object shell over a jitted numerical core written in a
  C-translatable Python subset; pluggable plotters; CZML export.
- **Trick** — Monte Carlo generation pattern (dispersed input files, execute
  anywhere) and the variable-server concept.
- **nyx-space** (AGPL — inspiration only) — trait-composition of dynamics,
  NIS/NEES filter validation, built-in multithreaded MC.
- **Basilisk** — message-bus architecture; C FSW modules portable to flight.
- **Astrobee** — docking-phase state machines. **SPHERES** — the
  guest-scientist "same code in sim and on hardware" API.
- **pyrpod** (GPL) — plume-impingement scope definition; we reimplement.
- **safe-autonomy-sims / run-time-assurance** (non-OSS license) — RTA
  (CBF/simplex) wrapper concepts for unverifiable controllers.

## Competitive landscape (18-month watch list)

- **Orekit MR !1133** (opened 2026-03): CW + Yamanaka-Ankersen + two-impulse
  rendezvous guidance — the incoming competitor, but Java-side and JVM-bound.
- **Basilisk's MuJoCo + ROS2 trajectory** — well-funded; strengthens the case
  for treating Basilisk as our backend rather than competing on raw dynamics.
- **SmallSatSim** (Apache-2.0, MuJoCo/MJX, 2026) — closest direct neighbor:
  contact-rich docking with GPU-parallel MC, but local-inertial dynamics, no
  validated orbital relative motion, no classical GNC layer, no codegen.
  Collaboration candidate more than competitor.

Differentiation strategy: analytic relative-motion rigor + RPOD-native
guidance API + verified flight codegen. Do **not** re-fight general-purpose
simulation frameworks.

## Survey table

| Framework | Language / license | Status (mid-2026) | RPOD-relevant strengths | Why not build the core on it |
|---|---|---|---|---|
| Basilisk (AVS Lab) | C/C++/Python, ISC | v2.10.2, very active | Coupled 6-DOF, formation/Hill FSW modules, sensors/actuators, bit-repeatable MC, MuJoCo contact path | No analytic CW/YA GNC layer; C++ friction; no codegen story |
| NASA 42 | C, NOSA-ish (ambiguous) | Alive, single maintainer, no releases | Multi-body attitude, contact forces, multi-SC prox ops | License ambiguity; edit-C-and-recompile extensibility |
| NASA Trick | C++, NOSA 1.3 | v25.1 active | Executive/MC/variable-server patterns | No dynamics included; NOSA is not OSI; C++-centric |
| NASA cFS/cFE | C, Apache-2.0 | v7.0.1 active | Flight heritage; the natural home for generated C GNC apps | It *is* flight software, not a sim |
| GMAT | C++, Apache-2.0 | R2026a | NASA-V&V'd force models, optimal control, OD | Attitude kinematic-only; RPOD absent; Python API unsuited to tight loops |
| Orekit | Java, Apache-2.0 | v13.1.6, most active astro lib | Best-in-class propagation/OD, CCSDS I/O; rendezvous MR pending | JVM boundary; RPOD not yet in a release |
| tudat/tudatpy | C++/Python, BSD-3 | v1.0 active | Torque-integrated rotational propagation; best Python guidance hooks | conda-only; no RPOD primitives; per-step GIL cost |
| poliastro/hapsira | Python, MIT | **dead** (archived 2023 / stalled 2024) | Design donor (units shell / jitted core split) | Unmaintained; no relative motion at all |
| nyx-space | Rust, **AGPL-3.0** | Very active, flight-proven | Trait composition, MC, OD validation | AGPL fatal for a permissive library |
| ANISE | Rust, MPL-2.0 | v0.10.3 active | SPICE-compatible frames/ephemerides, fast, thread-safe | (build on it — scope is frames/time only) |
| SPICE/spiceypy | C non-OSI / MIT wrapper | spiceypy active, CSPICE frozen | Kernel I/O, geometry | Non-OSI core; no dynamics |
| pyrpod | Python, GPL-3.0 | Small, unpackaged | Only open plume-impingement tool | GPL; prescribed profiles, not validated dynamics |
| Astrobee | C++/ROS1, Apache-2.0 | Maintained for ISS ops | Docking state machines, sim-to-hardware | ROS1 EOL; station-interior free-flyer, not orbital |
| safe-autonomy-sims (AFRL) | Python, "government owned" | Slowing | CW docking/inspection RL envs, RTA concepts | Not open source; CW-only, no contact |
| SmallSatSim (USC/JPL/ETH) | Python/JAX, Apache-2.0 | New (2026), tiny | Contact-rich docking, GPU MC, fault models | Local-inertial only; no orbital rel-motion fidelity, sensors, classical GNC |
| BSK-RL | Python, MIT | Active | Proof that pip-packaged layers over Basilisk work | Scheduling/inspection focus, not docking |
| TinyMPC / CVXPYgen | C / Python, MIT | Active (ICRA'24/'26) | Embedded MPC codegen incl. SOC glideslope constraints | (build on — guidance backend, not a sim) |
| heyoka.py | C++/Python, MPL-2.0 | v7.11 active | Taylor integration, free variational STMs | Scope is propagation, not RPOD |
| KSP tools (kRPC, kspdg) | mixed | kspdg active as AIAA challenge | Community visibility | Game physics; not validation-grade |

## Trajectory-optimization and verification landscape (update 2026-07-06)

The survey above is of *simulation frameworks*. Two adjacent landscapes matter
for what distinguishes Podium — the convex/sequential-convex guidance layer and
the exact-rational verification layer — and are surveyed at length in the tool
paper's related work. Neither is an RPOD library, and none couples the three
things Podium does (sequential-convex planners + exact-rational per-solve
certificates re-run in CI + a restricted Python-to-C flight-code path checked by
golden vectors, Frama-C/EVA, and CompCert).

### Sequential-convex / trajectory-optimization codebases

| Tool | Language / license | RPOD relevance | Why it is not the core |
|---|---|---|---|
| **SCPToolbox.jl** (UW-ACL) | Julia, MIT | Open SCvx / GuSTO / PTR / LCvx toolbox; ships planar-rendezvous and Apollo transposition-and-docking examples | Algorithms only — no certificate layer, no verified-code path, not RPOD-packaged; the closest open SCP-for-rendezvous prior art |
| **GuSTO.jl** (Stanford ASL) | Julia | SCP with convergence guarantees; free-flyer berthing hardware demo | Same: planner only |
| **OpenSCvx / ct-scvx** | Python/JAX | Continuous-time constraint-satisfaction SCvx; in-house C generation | No exact-rational certificate; floating-point |
| **successive\_rendezvous** (Malyuta et al.) | Python/CVXPY | Released 6-DoF Apollo transposition-and-docking via SCvx | Research code, planner only |
| **PIPG / factorization-free SCP** (Chari, Açıkmeşe) | Python | First-order rendezvous SCP solver | Solver/planner, no certificates or verified C |
| **STROM** (Kang, Yang) | Julia/CUDA | Certifiably-optimal nonconvex trajectory optimization via sparse moment-SOS SDP (chain-block) | Certificate is **floating-point** GPU SDP, not exact-rational; not RPOD |

### Verification / certificate tooling

| Tool | Language / license | Role | Why it is not the core |
|---|---|---|---|
| **RealCertify** | Maple | Exact rational SOS non-negativity certificates | General polynomials, not online, not RPOD, not in CI |
| **ValidSDP** | Coq | Machine-checked SDP/SOS positivity certificates | Proof-assistant checker for polynomial inequalities; not a GNC stack |
| **VSDP** | MATLAB/INTLAB | Rigorous interval bounds / infeasibility certificates for conic programs | Interval, not exact rational; not RPOD |
| **Credible autocoding** (Feron et al.); **verified MPC solver** (Davy, Cohen–Feron–Garoche) | C + Frama-C/ACSL | Autocode/verify convex-optimization solver code | Not RPOD, not exact-rational per-solve certificates, no CompCert stack |
| **Copilot / CopilotVerifier** (NASA/Galois) | Haskell DSL → C | Restricted-source → C monitors, Frama-C/CompCert-compatible, CI, bisimulation proof | Runtime monitors, not guidance; not RPOD |
| **CVXPYgen** | Python → C | Embedded C from CVXPY convex problems | Codegen only, no formal verification or certificates |
| **CvxLean** | Lean 4 | Verified convex-optimization problem reductions | Formal reductions, not an RPOD certificate stack |
| **EIQP** (Wu, Xiao, Braatz) | C | Execution-time-certified, infeasibility-detecting QP solver | Solver with a WCET bound, no per-solve exact certificate |

Podium's contribution relative to both landscapes is the *integration*: an
open, RPOD-focused library that couples a sequential-convex planner, exact-rational
certificates re-run in continuous integration (tolerance-free for the barrier,
Lyapunov, sum-of-squares, and optimality-gap certificates, and an exact
suboptimality bound for KKT), and a restricted Python-to-C flight-code path
checked by golden vectors, Frama-C/EVA, and CompCert.
See the tool paper (`docs/paper/`) for the full comparison and citations.
