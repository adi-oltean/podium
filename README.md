# Podium

[![Paper DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21225267.svg)](https://doi.org/10.5281/zenodo.21225267)

**Physics-precise RPOD (rendezvous, proximity operations & docking) GNC library, simulation sandbox, and visualization — focused on LEO/MEO.**

Podium is an open-source framework for *developing, testing, and formally validating* guidance, navigation & control algorithms for spacecraft rendezvous and docking. It is Python-first for iteration speed, but the algorithm core is written in a deliberately restricted style (the *static subset*) so that flight algorithms translate mechanically to C, are proven runtime-error-free by an integrated sound static analyzer (Frama-C/EVA), and compile through a formally-verified compiler (CompCert).

> Status: v0.1–v0.7 complete, v0.8 in progress. The verified-flight-code path runs end to end and is checked by **eight CI lanes**. Highlights: the CW / **Yamanaka-Ankersen** / **relative-orbital-elements** kernels (Koenig Keplerian/J2/J2+drag STMs, e/i-vector passive safety — no other shipped open implementation exists); the nonlinear truth model (J2 + differential drag + seeded stochastic atmosphere) **cross-validated against Orekit**; **convex Layer-0 guidance** (DPP-compiled impulsive planning on CW/YA/ROE, LCvx finite-burn with losslessness audit, constraint library) and **SCP Layer-1** docking (PTR with continuous-time cuts, 6-DOF attitude-coupled planning with a body-fixed thruster); a Joseph-form **relative-nav EKF**; rigid-body **attitude** with the full **environmental-torque suite** (gravity gradient, aerodynamic, SRP, magnetic — analytically validated); a **C99 emitter** (20 flight kernels, ACSL-annotated) whose output is **bit-exact vs Python** for scalar kernels (correctly-rounded via CORE-MATH; matmul agrees up to reassociation), bit-identical on aarch64, checked through **CompCert** and **proven alarm-free by Frama-C/EVA**; **exact-rational certificates** for abort-safety barriers and online-solver KKT re-verification (QP suboptimality bounds; SOCP conic-KKT re-check); a **reference mission** shipped per-tag as an evidence-gated audit bundle; and the **ARCH-COMP rendezvous benchmark as a CI reachability gate** re-proving closed-loop safety on every relevant commit. See [`docs/verification.md`](docs/verification.md).

## Why another space simulator?

Mature tools exist — Basilisk, NASA 42, Trick, Orekit, GMAT — and Podium does not re-fight their battles (see the [comparative analysis](docs/comparative-analysis.md)). What none of them offer together:

1. **RPOD as the first-class problem.** Relative motion (CW/Tschauner-Hempel), approach corridors, keep-out zones, plume impingement, passive abort safety, and docking contact — not an afterthought bolted onto an orbit propagator.
2. **A verification-ready algorithm core.** GNC algorithms written as pure, statically-shaped, bounded-loop step functions with machine-readable contracts — the style that abstract interpretation (Astrée-class tools) can actually prove things about, and that translates line-for-line to embedded C.
3. **Convex trajectory optimization built in.** Direct LP/SOCP transcription on the exact CW/YA/ROE discretizations (DPP-compiled, Clarabel), lossless-convexification finite-burn planning shipped with validity audits rather than assumptions, a constraint library (approach cone, rotating-hyperplane KOZ, plume, Breger-How passive safety), successive convexification (PTR with continuous-time cuts) for the nonconvex docking layer, and 6-DOF attitude-coupled planning with a body-fixed thruster — plus an exact-rational KKT checker that bounds the online solver's suboptimality (when a valid dual point exists) with no trust in its floating-point.
4. **A sandbox you can trust.** Deterministic fixed-step simulation (bit-identical replays, enforced by test), truth/flight separation, seeded noise and stochastic atmosphere, STL-robustness spec oracles, MuJoCo probe-drogue contact, and cross-validation against Orekit (translational) and exact analytic solutions (attitude).
5. **Verification as a regression, not a ceremony.** Eight CI lanes re-prove safety on every relevant commit: reachability (JuliaReach), exact-rational barrier certificates, sound static analysis (Frama-C/EVA), golden vectors through CompCert and on aarch64, and Orekit cross-validation — with every tagged release shipping a byte-deterministic, evidence-gated audit bundle.

## Layout

```
src/podium/
  core/        Verifiable algorithm core (static subset): CW dynamics & STM,
               Yamanaka-Ankersen STM (elliptic orbits), relative orbital
               elements (Koenig Keplerian/J2/J2+drag STMs, LVLH maps,
               control matrix), quaternion kernel, fixed-step integrators
  dynamics/    Truth models: nonlinear dual-ECI relative motion, J2 + drag
               (+ seeded stochastic atmosphere), eccentric-valid ROE map,
               rigid-body attitude (Euler+quaternion RK4) with the
               environmental-torque suite (gravity gradient, aero, SRP,
               magnetic) and a disturbance aggregator
  guidance/    Glideslope, convex Layer-0 (impulsive DPP planners on
               CW/YA/ROE, LCvx finite-burn + audits, constraint library,
               MIB bridge), SCP Layer-1 PTR docking (CTCS cuts, STL
               timed-reach), 6-DOF attitude-coupled PTR (body-fixed
               thruster), tumbling-target study, ARCH benchmark model
  control/     LQR (discrete Riccati + continuous CARE synthesis),
               quaternion-feedback attitude, pulsed docking control,
               push-only thruster allocation (min-propellant LP)
  nav/         Relative-navigation EKF (Joseph form, linear + bearing/range
               updates, sequential scalar update), sensor models
               (relative GNSS, camera, lidar)
  sim/         Deterministic engine (actuator MIB/execution-error truth),
               spec registry (STL robust semantics), IDSS contact-box
               checkers, MuJoCo probe-drogue contact, seeded Monte Carlo,
               the end-to-end reference mission + audit bundle
  verify/      Contracts (ranges/invariants → ACSL), exact-rational
               barrier certificates AND KKT suboptimality-bound
               checkers (QP + SOCP), correctly-rounded transcendental oracle
  emit/        C99 emitter for the static subset (bounded loops, matmul,
               ACSL rendering), the CORE-MATH correctly-rounded option,
               the EVA driver generator, and the cFS app generator;
               golden-vector equivalence vs Python (bit-exact scalars/sqrt,
               matmul up to reassociation)
tools/         reach/ (JuliaReach CI gate), eva_gate.py (Frama-C/EVA),
               tier2_build_run.sh (aarch64/qemu), build_audit_bundle.py,
               deploy_viewer.py, UI Playwright suites
tests/         pytest receipts (truth-model validations, closed-loop
               flights, statistical consistency, exact-arithmetic audits)
examples/      Runnable scenarios (V-bar approach); cfs_nav_app/ — a
               Core Flight System app running the verified kernels
third_party/   Vendored CORE-MATH (correctly-rounded sin/cos, MIT)
docs/          Architecture, comparative analysis, verification approach,
               numerical reproducibility, visualization, per-issue plans
```

**Live demos:**
- [V-bar approach viewer](https://adi-oltean.github.io/podium/) — a 1 km glideslope approach propagated through the nonlinear truth model, in a self-contained page (fermi-style: no build system, no external requests).
- [3-D viewer](https://adi-oltean.github.io/podium/3d/) — the same approach in an interactive three.js LVLH scene: approach corridor, keep-out sphere, burn glyphs, follow/orbit cameras (vendored three.js, zero build).
- [ISS-sim autopilot](https://adi-oltean.github.io/podium/iss-sim/) — Podium's pulsed docking-control laws (`podium.control.docking`) flying [SpaceX's ISS docking simulator](https://iss-sim.spacex.com/) via a paste-in-console autopilot.


## Quick start

```bash
pip install -e ".[dev,viz]"
pytest
python examples/vbar_approach.py
```

```python
import numpy as np
from podium.core import cw

n = cw.mean_motion(3.986004418e14, 6_778_137.0)   # 400 km target orbit
x0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])  # 1 km behind on V-bar
dv1, dv2 = cw.two_impulse(x0, np.zeros(6), n, 1500.0)
```

## Design principles

- **LVLH conventions:** x radial (zenith), y along-track, z cross-track; SI units everywhere; quaternions scalar-first, body→reference.
- **Truth vs. flight separation:** truth models may use anything in SciPy; flight algorithms live in `podium.core` under the static-subset rules ([`docs/verification.md`](docs/verification.md)) and are exercised in the sim through the same step-function interface they will have in C.
- **Determinism:** fixed-step master clock, seeded noise, no wall-clock or platform dependence in results.
- **Contracts as data:** core functions declare input ranges/invariants via `podium.verify` (with explicit `DEFAULT_RANGES` operating assumptions recorded in the generated driver for any remaining gaps); checked at runtime in the sandbox, and rendered as ACSL preconditions on the emitted C that Frama-C/EVA discharges in CI.

## Documentation

| Document | Contents |
|---|---|
| [`docs/architecture.md`](docs/architecture.md) | Frames, units, module boundaries, sim engine design |
| [`docs/comparative-analysis.md`](docs/comparative-analysis.md) | Survey of existing simulators; build-on vs interop decisions |
| [`docs/paper/podium-paper.pdf`](docs/paper/podium-paper.pdf) | The tool paper describing the library (DOI [10.5281/zenodo.21225267](https://doi.org/10.5281/zenodo.21225267)) |
| [`docs/exact-arithmetic-certificates/note.pdf`](docs/exact-arithmetic-certificates/note.pdf) | Technical note: constructions, proofs, and prior art for the exact-rational certificates (DOI [10.5281/zenodo.21247381](https://doi.org/10.5281/zenodo.21247381)) |
| [`docs/optimality-gap-certificates.md`](docs/optimality-gap-certificates.md) | Index of the optimality-gap results mapped to code and tests |
| [`docs/verification.md`](docs/verification.md) | The ten shipped verification modalities, static-subset rules, contract→ACSL mapping, layered assurance |
| [`docs/numerical-reproducibility.md`](docs/numerical-reproducibility.md) | Golden-vector methodology, equality classes, and cross-ISA bit-exactness conditions |
| [`docs/visualization.md`](docs/visualization.md) | Rendering architecture (patterns adopted from fermi) |
| [`docs/plans/`](docs/plans/) | One design/receipt plan per numbered issue |

## Paper and citation

The library is described in the tool paper:

> Adi Oltean. *Podium: An Open-Source Library for Rendezvous and Docking Guidance, Navigation, and Control with Integrated Formal Verification.* 2026. [`docs/paper/podium-paper.pdf`](docs/paper/podium-paper.pdf). DOI: [10.5281/zenodo.21225267](https://doi.org/10.5281/zenodo.21225267).

The mathematics behind the exact-rational certificates — the constructions, proofs, and prior-art positioning for the barrier, KKT, control-Lyapunov, sum-of-squares, and optimality-gap certificates — is collected in a companion technical note, itself citable:

> Adi Oltean. *Exact-Rational Certificates in Podium: Constructions, Proofs, and Prior Art.* Technical note, 2026. [`docs/exact-arithmetic-certificates/note.pdf`](docs/exact-arithmetic-certificates/note.pdf). DOI: [10.5281/zenodo.21247381](https://doi.org/10.5281/zenodo.21247381).

The optimality-gap results are also indexed to their code and tests in [`docs/optimality-gap-certificates.md`](docs/optimality-gap-certificates.md).

To cite the software itself, use [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button); the code is archived on Zenodo (software DOI [10.5281/zenodo.21225268](https://doi.org/10.5281/zenodo.21225268)), with a version DOI per tagged release.

## License

[MIT](LICENSE) — maximally permissive and simple. Note for algorithm contributions: some convexification methods in this space are patent-encumbered upstream (e.g., G-FOLD, US 8,489,260 — Caltech/JPL); an MIT license neither grants nor affects third-party patent rights, so implementations here follow the published papers.
