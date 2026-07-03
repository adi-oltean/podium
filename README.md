# Podium

**Physics-precise RPOD (rendezvous, proximity operations & docking) GNC library, simulation sandbox, and visualization — focused on LEO/MEO.**

Podium is an open-source framework for *developing, testing, and formally validating* guidance, navigation & control algorithms for spacecraft rendezvous and docking. It is Python-first for iteration speed, but the algorithm core is written in a deliberately restricted style (the *static subset*) so that flight algorithms translate mechanically to C and can be proven safe by an external abstract-interpretation tool.

> Status: pre-alpha scaffold. The CW kernel, **Yamanaka-Ankersen STM** (elliptic relative motion — no other open implementation exists), quaternion kernel, glideslope guidance, and LQR synthesis/application are implemented and tested; the rest of the tree is designed but stubbed. See [`docs/roadmap.md`](docs/roadmap.md).

## Why another space simulator?

Mature tools exist — Basilisk, NASA 42, Trick, Orekit, GMAT — and Podium does not re-fight their battles (see the [comparative analysis](docs/comparative-analysis.md)). What none of them offer together:

1. **RPOD as the first-class problem.** Relative motion (CW/Tschauner-Hempel), approach corridors, keep-out zones, plume impingement, passive abort safety, and docking contact — not an afterthought bolted onto an orbit propagator.
2. **A verification-ready algorithm core.** GNC algorithms written as pure, statically-shaped, bounded-loop step functions with machine-readable contracts — the style that abstract interpretation (Astrée-class tools) can actually prove things about, and that translates line-for-line to embedded C.
3. **Convex trajectory optimization built in.** Direct SOCP transcription of relative dynamics plus successive convexification for the nonconvex constraints (keep-out zones, plume), following the G-FOLD / SCvx lineage (planned; see the roadmap).
4. **A sandbox you can trust.** Deterministic fixed-step simulation (bit-identical replays), truth/flight separation, seeded Monte Carlo, and cross-validation hooks against established stacks.

## Layout

```
src/podium/
  core/        Verifiable algorithm core (static subset): CW dynamics & STM,
               Yamanaka-Ankersen STM (elliptic orbits), quaternion kernel,
               fixed-step integrators
  dynamics/    Truth models: Tschauner-Hempel, nonlinear relative motion,
               J2 + drag, rigid-body attitude
  guidance/    Glideslope, multi-impulse targeting, convex/SCP trajectory
               optimization, passive-safety checks
  control/     LQR (offline synthesis / flight-side gain application),
               attitude control, thruster allocation
  nav/         Relative-navigation EKF, sensor models (RGPS, camera, lidar)
  sim/         Deterministic engine, events (contact/KOZ/abort), Monte Carlo
  viz/         Analysis plots + interactive three.js viewer
  verify/      Contracts (input ranges, invariants) + export to the external
               abstract-interpretation validation tool
tests/         pytest + hypothesis property tests
examples/      Runnable scenarios (V-bar approach, ...)
docs/          Architecture, comparative analysis, trajectory optimization,
               verification approach, roadmap
```

**Live demos:**
- [V-bar approach viewer](https://adi-oltean.github.io/podium/) — a 1 km glideslope approach propagated through the nonlinear truth model, in a self-contained page (fermi-style: no build system, no external requests).
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
| [`docs/roadmap.md`](docs/roadmap.md) | Milestones toward v0.1 |

## License

[MIT](LICENSE) — maximally permissive and simple. Note for algorithm contributions: some convexification methods in this space are patent-encumbered upstream (e.g., G-FOLD, US 8,489,260 — Caltech/JPL); an MIT license neither grants nor affects third-party patent rights, so implementations here follow the published papers.
