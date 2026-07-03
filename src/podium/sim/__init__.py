"""Deterministic simulation engine.

Design (see docs/architecture.md):

- **Fixed-step master clock.** The GNC stack runs at a fixed rate (default
  10 Hz) against truth dynamics integrated at an integer multiple of that
  rate. Determinism is a hard requirement: identical inputs give bit-identical
  trajectories, which makes Monte Carlo campaigns and regression tests exact.
- **Truth vs. flight separation.** Truth models (nonlinear dynamics, sensor
  noise, actuator imperfections) use full Python/SciPy. Flight algorithms are
  called through the same pure step-function interface they will have in C.
- **Events.** Contact/capture detection, keep-out-zone violation, abort
  triggers — evaluated on the master clock with interval refinement.
- **Monte Carlo.** Seeded dispersion campaigns over initial state, sensor and
  actuator errors; results as structured arrays for batch analysis.
- **Specs.** Named requirements (STL robust semantics, PUS-12-shaped base
  fragment) evaluated over trace channels; margins double as pytest oracles.
- **Plots.** `podium.sim.plots` (import explicitly; matplotlib optional):
  trajectory plane view, channel time series, dv timeline.
- **Acceptance + campaigns.** `podium.sim.idss` (IDSS Rev G contact box
  checkers) and `podium.sim.monte_carlo` (seeded campaigns, structured
  arrays, exact reproducibility).
- **Contact.** `podium.sim.contact` (import explicitly; mujoco optional):
  probe-drogue capture truth model with envelope sweeps.
"""

from podium.sim import idss, monte_carlo, spec  # noqa: F401
from podium.sim.engine import (  # noqa: F401
    Controller,
    Scenario,
    Trace,
    circular_target,
    mean_motion_of,
    run,
)
