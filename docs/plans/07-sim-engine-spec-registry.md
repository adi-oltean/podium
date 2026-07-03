# 07 — Sim engine + spec registry v0

GitHub issue: https://github.com/adi-oltean/podium/issues/7
(filed after a network outage; number to be confirmed)

## Problem

Podium has no closed loop: every test hand-rolls its own propagation, no
engine hosts the GNC step functions against the truth model, determinism
is a convention rather than an enforced engine property, and there is no
event/monitor layer. The spec registry must be co-designed with the
engine's event system (same design decision) to avoid retrofitting the
one-spec-many-consumers architecture later.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/sim/engine.py` | MISS | Scenario, Trace, run() loop, viewer export |
| `src/podium/sim/spec.py` | MISS | registry v0: named predicates, robust margins |
| `tests/test_sim.py` | MISS | determinism, closed-loop receipts, spec semantics |
| `docs/roadmap.md`, `docs/architecture.md` | PARTIAL | mark items, note interfaces |
| `pyproject.toml` | PARTIAL | `slow` pytest marker registration |

## Fix

Engine (`podium.sim.engine`):
- Fixed-step master GNC clock (`dt_gnc`, default 0.1 s architectural /
  1-2 s in tests); truth propagates between ticks with `truth_substeps`
  RK4 steps of the dual-ECI truth model.
- Flight-block interface v0: `controller(t, measured_x_rel) -> dv_lvlh`
  — impulsive Δv at ticks (matches glideslope, LQR-as-impulse, pulsed
  docking). Continuous-thrust interface deferred.
- Measurement = true relative state + seeded Gaussian noise (one
  `np.random.default_rng(seed)`; the ONLY randomness in the engine).
- Trace: times, relative states, target ECI, burn log; named channels
  (range, range-rate, speed, components); `crossing_times()` helper with
  linear interpolation; `to_viewer_json()` matching the live viewer's
  schema (closes part of the export-API roadmap item).
- Determinism receipt: identical scenario+seed ⇒ `np.array_equal` traces.

Spec registry (`podium.sim.spec`):
- `Spec` = named requirement over a trace channel with STL robust
  semantics for the base fragment (PUS-12-shaped): `always_above/below/
  between`, `eventually_below/above`, `final_between`, all with optional
  time windows; margin > 0 ⇔ satisfied, magnitude = robustness.
- `evaluate(specs, channels)` -> margins dict; engine attaches margins to
  the Trace; pytest oracles assert margins.
- rtamt adapter for full STL is a follow-up (network outage during
  implementation; also keeps the core dependency-free).

Receipts:
- bit-identical replay (same seed), reproducible-noise checks;
- glideslope flown closed-loop through the engine against the nonlinear
  truth (arrival + corridor/terminal specs as margins);
- LQR gains synthesized on the CW model stabilize the nonlinear truth
  through the engine (impulse-equivalent actuation);
- spec semantics pinned on synthetic traces; crossing-time
  self-consistency.

## Acceptance Criteria

- [ ] Engine + registry implemented; suite green (pytest, ruff, mypy)
- [ ] Determinism enforced by test
- [ ] Closed-loop glideslope + LQR receipts green
- [ ] Viewer JSON export schema-compatible
- [ ] `slow` marker registered and applied to the heavy truth tests

## Push/merge instructions

Single commit on main: `07 — Sim engine + spec registry v0 (#7)`; push;
close the issue (file it first if the outage persisted).

## Verification steps

Full suite; run the glideslope scenario and load the exported JSON in the
live viewer.
