# 13 — v0.1 stragglers: stochastic atmosphere + matplotlib plots

GitHub issue: https://github.com/adi-oltean/podium/issues/13

## Problem

Two roadmap items kept v0.1 open: Monte Carlo without density dispersion
understates the dominant LEO uncertainty, and traces had no analysis
plots beyond the web viewer.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/dynamics/nonlinear.py` | PARTIAL | DensityPerturbation; time-dependent density plumbed through accel/deriv |
| `src/podium/sim/plots.py` | MISS | trajectory / channels / dv figures |
| `tests/test_stochastic_atmo.py`, `tests/test_plots.py` | MISS | receipts |
| `pyproject.toml` | PARTIAL | matplotlib into dev extras (viz extra exists) |

## Fix

- `DensityPerturbation`: OU log-density factor, EXACT discretization
  (p+ = e^{-dt/tau} p + sigma sqrt(1-phi^2) xi) precomputed on a fixed
  grid at construction from one seeded Generator, linearly interpolated,
  clamped outside the grid — density(t) is a deterministic function of
  time, so RK4 and bit-identical replay survive. Default calibration
  sigma_log=0.35 (=> +2sigma ~ 2.0x, inside the +50-125% storm band),
  tau = 6 h. `DragConfig.density(h, t)`; t threaded through
  perturb_accel/total_accel/_deriv (defaults keep old call sites valid).
  Known approximation: the LVLH transform's frame-precession term uses
  t=0 density; drag contributes ~nothing out-of-plane, so the effect is
  nil.
- `podium.sim.plots`: object-API figures (no pyplot, no backend state):
  plot_trajectory (y/x plane + burn markers), plot_channels (stacked
  series), plot_dv (stems + cumulative). Not imported by podium.sim —
  matplotlib stays optional.

Receipts: seeded determinism; OU stationary std and lag-tau
autocorrelation pinned; +2sigma calibration inside the band; sigma=0
bit-exactly the baseline; clamping (fixed a real int()-truncation bug
where small negative t extrapolated); 3-orbit decay ratio bounded by the
realized factor extremes with bit-identical replay across fresh objects;
plot data arrays match trace channels; PNGs save non-empty headless.

## Acceptance Criteria

- [x] All receipts green (suite, ruff, mypy)
- [x] v0.1 marked COMPLETE in the roadmap

## Push/merge instructions

Single commit on main: `13 — v0.1 complete: stochastic atmosphere +
analysis plots (#13)`; push; close.

## Verification steps

Full suite; render the three figures from a demo trace and eyeball.
