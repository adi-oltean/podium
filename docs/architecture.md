# Architecture

## Scope

Initial focus: **RPOD in LEO/MEO** — far-range rendezvous handoff (~10 km) through proximity operations, final approach, and docking/berthing contact. Deep-space and cislunar RPOD are out of scope until v1.

## Frames and conventions

| Item | Convention |
|---|---|
| Relative frame | Target-centered LVLH: **x radial (zenith), y along-track, z cross-track** (right-handed; z along orbit normal) |
| Inertial frame | ECI (J2000/GCRF) for truth propagation |
| Attitude | Quaternions, **scalar-first** `[w, x, y, z]`, rotate body→reference |
| Units | SI throughout (m, s, kg, rad). No mixed-unit APIs, ever |
| Time | Mission elapsed seconds (float64) inside the sim; epoch handling only at the truth-model boundary |

## The two-layer rule

The single most important boundary in the codebase:

```
┌──────────────────────────────────────────────────────────────┐
│  SANDBOX (full Python)                                       │
│  truth dynamics · sensor/actuator error models · Monte Carlo │
│  gain synthesis · cvxpy prototyping · plotting · viz export  │
│        │                        ▲                            │
│        │ sensor outputs         │ actuator commands          │
│        ▼                        │                            │
│  ┌──────────────────────────────────────────┐                │
│  │  FLIGHT CORE (static subset, podium.core)  │ → C + external │
│  │  pure step functions, fixed shapes,      │   abstract-    │
│  │  bounded loops, contracts                │   interpretation│
│  └──────────────────────────────────────────┘   validation   │
└──────────────────────────────────────────────────────────────┘
```

Everything inside the flight core obeys the static subset (see `verification.md`); everything outside may use the full scientific-Python stack. The sim engine calls flight-core functions through the exact interface they will have after C translation: `step(state_in, inputs, params) -> (state_out, outputs)` with fixed-shape arrays.

Some algorithms straddle the line deliberately — e.g. LQR: the Riccati recursion (synthesis) is sandbox-side and produces a constant gain; only `apply_gain` (a saturated matrix-vector product) is flight-side. The same split applies to trajectory optimization: cvxpy prototyping sandbox-side, with the online solves covered by the certificate layer rather than the verified-C path (embedded fixed-iteration solver codegen is a partial, offline prototype).

## Simulation engine

- **Fixed-step master clock** at the GNC rate (default 10 Hz). Truth dynamics integrate at an integer multiple (default 10× → 100 Hz RK4). No adaptive stepping in the loop — determinism and replayability outrank integrator elegance; adaptive SciPy integrators are used only offline to bound fixed-step truncation error.
- **Determinism:** identical config + seed ⇒ bit-identical trajectories. All randomness flows from one seeded generator; no wall-clock, no dict-ordering dependence, no platform-dependent reductions.
- **Events** are evaluated on the master clock, with linear interpolation for threshold crossings via `Trace.crossing_times` (between-step bisection refinement is not yet implemented): docking-interface contact, keep-out-zone entry, approach-corridor departure, abort triggers.
- **Monte Carlo:** dispersion campaigns are (config, seed-list) pairs; results are structured arrays written to `.npz`. Because single runs are deterministic, any MC outlier replays exactly for debugging.

## Truth models (sandbox side)

- Nonlinear two-body relative motion in the target LVLH frame (exact, no linearization), with **J2** and **exponential-atmosphere drag** — the dominant LEO perturbations for RPOD timescales. Differential drag matters for dissimilar chaser/target ballistic coefficients. An optional seeded mean-reverting (OU) density perturbation covers storm-class dispersion while preserving bit-identical replay (exact discretization on a fixed time grid).
- Tschauner-Hempel (Yamanaka-Ankersen STM) for eccentric targets; CW for near-circular.
- Relative orbital elements (quasi-nonsingular ROE) with Koenig closed-form Keplerian/J2 STMs for perturbed multi-orbit relative motion, plus the e/i-vector passive-safety geometry (`podium.guidance.safety`). ROE STMs propagate mean elements; the osculating/mean distinction is handled truth-side and documented in the tests.
- Rigid-body attitude with reaction-wheel and thruster torques; thruster minimum-impulse-bit and rise/tail-off shaping.
- Sensor models with error budgets: relative GNSS, docking camera (bearing + fiducial pose), lidar. (Gyro / star-tracker models are not yet implemented.)

Fidelity hierarchy: every linearized model used by guidance is validated in CI against the nonlinear truth model, and the truth model itself is cross-validated against external references (Orekit test vectors for translational motion and analytic oracles for attitude) — errors quantified, not assumed.

## GNC dataflow per tick

```
sensors(truth) → nav filter (EKF, core) → guidance (waypoint/optimal traj, core)
    → control (tracking + attitude, core) → allocation (thruster map, core)
    → actuators(truth) → dynamics(truth) → …
```

Each core block is a pure function of `(state, measurement/reference, params)`. Blocks are individually replaceable — swapping a glideslope for an SCP-generated trajectory changes one constructor argument in the scenario config.

Implemented so far along this dataflow: the deterministic engine (`podium.sim.engine`) with actuator MIB/cap/execution-error truth; the sensor truth models (`podium.nav.sensors`: relative GNSS, camera, lidar); the relative-nav EKF (`podium.nav.ekf`: Joseph-form update, white-noise-acceleration process model, commanded burns fed through prediction as known inputs; NEES/NIS consistency enforced by test); LQR/glideslope/pulsed-docking control and push-only thruster allocation; the convex Layer-0 planners (impulsive and LCvx finite-burn) and the SCP Layer-1 / 6-DOF attitude-coupled docking planners; and rigid-body attitude with the environmental-torque suite. The continuous-thrust engine interface and between-step event bisection remain open.

## Mission phases (typical LEO profile)

1. **Far-range rendezvous** (~30–10 km): ground/absolute-nav targeting, CW two-impulse transfers, hop sequences.
2. **Close rendezvous** (10 km–500 m): relative GNSS nav, safety ellipse / football orbits, hold points.
3. **Final approach** (500 m–10 m): V-bar or R-bar approach, glideslope or optimized trajectory, approach corridor + KOZ constraints, camera/lidar nav.
4. **Docking** (10 m–contact): sub-cm/s closing rates, tight lateral/angular corridors, plume-inhibited zones, capture envelope of the docking mechanism (IDSS-like numbers as defaults).

Passive abort safety is a cross-cutting constraint: from any point of a nominal approach, thrust-free coasting (or a single canned abort burn) must not intersect the keep-out zone for N orbits.
