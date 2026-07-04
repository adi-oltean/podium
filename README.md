# Podium

**Physics-precise RPOD (rendezvous, proximity operations & docking) GNC library, simulation sandbox, and visualization — focused on LEO/MEO.**

Podium is an open-source framework for *developing, testing, and formally validating* guidance, navigation & control algorithms for spacecraft rendezvous and docking. It is Python-first for iteration speed, but the algorithm core is written in a deliberately restricted style (the *static subset*) so that flight algorithms translate mechanically to C and can be proven safe by an external abstract-interpretation tool.

> Status: alpha — v0.1 and v0.2 complete, v0.3 in progress. Implemented and tested: the CW kernel, **Yamanaka-Ankersen STM** and **relative-orbital-elements kernel** (Koenig Keplerian/J2/J2+drag STMs, e/i-vector passive safety — no other open implementations exist), the nonlinear truth model (J2 + differential drag + seeded stochastic atmosphere), the deterministic sim engine with an STL-semantics spec registry, **convex Layer-0 guidance** (DPP-compiled impulsive planning on CW/YA/ROE, LCvx finite-burn with losslessness audit and primer certificate, approach-cone/KOZ/plume/passive-safety constraints, MIB quantization bridge), LQR (discrete + continuous CARE), pulsed docking control, a Joseph-form **relative-nav EKF**, and the **ARCH-COMP rendezvous benchmark wired as a CI reachability gate** — every guidance/control/dynamics commit re-proves LOS-cone, velocity, and abort-avoidance safety with JuliaReach, for both the published reference controller and Podium-synthesized gains. See [`docs/roadmap.md`](docs/roadmap.md).

## Why another space simulator?

Mature tools exist — Basilisk, NASA 42, Trick, Orekit, GMAT — and Podium does not re-fight their battles (see the [comparative analysis](docs/comparative-analysis.md)). What none of them offer together:

1. **RPOD as the first-class problem.** Relative motion (CW/Tschauner-Hempel), approach corridors, keep-out zones, plume impingement, passive abort safety, and docking contact — not an afterthought bolted onto an orbit propagator.
2. **A verification-ready algorithm core.** GNC algorithms written as pure, statically-shaped, bounded-loop step functions with machine-readable contracts — the style that abstract interpretation (Astrée-class tools) can actually prove things about, and that translates line-for-line to embedded C.
3. **Convex trajectory optimization built in.** Direct LP/SOCP transcription on the exact CW/YA/ROE discretizations (DPP-compiled, Clarabel), lossless-convexification finite-burn planning shipped with validity audits rather than assumptions, and a constraint library (approach cone, rotating-hyperplane KOZ, plume, Breger-How passive-safety scenarios) — successive convexification for the remaining nonconvexities is the v0.4 layer.
4. **A sandbox you can trust.** Deterministic fixed-step simulation (bit-identical replays, enforced by test), truth/flight separation, seeded noise and stochastic atmosphere, STL-robustness spec oracles, and cross-validation hooks against established stacks.
5. **Verification as a regression, not a ceremony.** The ARCH-COMP rendezvous benchmark runs as a CI gate: reachability analysis re-proves closed-loop safety properties on every relevant commit — including for controllers synthesized by Podium's own LQR machinery.

## Layout

```
src/podium/
  core/        Verifiable algorithm core (static subset): CW dynamics & STM,
               Yamanaka-Ankersen STM (elliptic orbits), relative orbital
               elements (Koenig Keplerian/J2/J2+drag STMs, LVLH maps,
               control matrix), quaternion kernel, fixed-step integrators
  dynamics/    Truth models: nonlinear dual-ECI relative motion, J2 + drag
               (+ seeded stochastic atmosphere), eccentric-valid ROE map,
               rigid-body attitude (planned)
  guidance/    Glideslope, convex Layer-0 (impulsive DPP planners on
               CW/YA/ROE, LCvx finite-burn + audits, constraint library,
               MIB bridge), passive-safety metrics, ARCH benchmark model
  control/     LQR (discrete Riccati + continuous CARE synthesis,
               flight-side gain application), pulsed docking control,
               CW/YA ZOH discretizations
  nav/         Relative-navigation EKF (Joseph form, linear + bearing/range
               EKF updates), sensor models (relative GNSS, camera, lidar)
  sim/         Deterministic engine (with actuator MIB/execution-error
               truth), spec registry (STL robust semantics), IDSS
               contact-box checkers, seeded Monte Carlo, analysis plots
  viz/         (viewer/ on Pages): canvas 2-D viewer, three.js 3-D viewer,
               ISS-sim autopilot page — all zero-build, self-contained
  verify/      Contracts (input ranges, invariants), exact-rational barrier
               certificates, + export to the external validation tool
  emit/        C99 emitter for the static subset (ACSL + analyzer
               annotations; bit-exact tier-1 golden vectors vs Python)
tools/reach/   JuliaReach reachability regression (CI gate)
tests/         pytest receipts (truth-model validations, closed-loop
               flights, statistical consistency, audits)
examples/      Runnable scenarios (V-bar approach); cfs_nav_app/ — a
               Core Flight System app running the verified kernels
docs/          Architecture, comparative analysis, trajectory optimization,
               verification approach, roadmap, per-issue plans
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
- **Contracts as data:** every core function declares input ranges/invariants via `podium.verify`; checked at runtime in the sandbox, emitted as annotations for the external prover in the C translation.

## Documentation

| Document | Contents |
|---|---|
| [`docs/architecture.md`](docs/architecture.md) | Frames, units, module boundaries, sim engine design |
| [`docs/comparative-analysis.md`](docs/comparative-analysis.md) | Survey of existing simulators; build-on vs interop decisions |
| [`docs/verification.md`](docs/verification.md) | Static subset rules, contract→annotation mapping, validation flow |
| [`docs/visualization.md`](docs/visualization.md) | Rendering architecture (patterns adopted from fermi) |
| [`docs/roadmap.md`](docs/roadmap.md) | Milestones and per-release status (v0.1/v0.2 complete, v0.3 current) |
| [`docs/plans/`](docs/plans/) | One design/receipt plan per numbered issue |

## License

[MIT](LICENSE) — maximally permissive and simple. Note for algorithm contributions: some convexification methods in this space are patent-encumbered upstream (e.g., G-FOLD, US 8,489,260 — Caltech/JPL); an MIT license neither grants nor affects third-party patent rights, so implementations here follow the published papers.
