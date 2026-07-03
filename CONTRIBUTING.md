# Contributing

Contributions are welcome — bug reports, dynamics/GNC algorithm
implementations, validation cases against flight data or other simulators,
and documentation all help.

## Ground rules

- **`src/podium/core` follows the static subset** documented in
  [`docs/verification.md`](docs/verification.md): pure functions, fixed
  shapes, compile-time loop bounds, contracts on every scalar input. PRs into
  `core` that break these rules will be asked to move the code to the sandbox
  layer instead.
- **Physics claims need receipts.** New dynamics or guidance implementations
  must include tests against an independent reference: a closed-form solution,
  a published numerical example (cite it), or cross-validation against the
  nonlinear truth model with quantified error.
- **Determinism is non-negotiable** in `sim`: no wall-clock, no unseeded
  randomness, no platform-dependent reductions.
- SI units, LVLH x-radial/y-along-track/z-cross-track, scalar-first
  quaternions — no exceptions, no per-module conventions.

## Workflow

```bash
pip install -e ".[dev]"
ruff check src tests examples
mypy
pytest
```

All three must pass. Property-based tests (hypothesis) are encouraged for
kernel math. Keep functions small; the verification story depends on it.

Multi-orbit truth-model propagations are marked `slow`; use
`pytest -m "not slow"` for the fast lane (<5 s) during development. CI
runs the full suite.

## Documentation upkeep

Code changes are not done until the docs that describe them are current.
These files describe the *state* of the project and must be refreshed
whenever the code they summarize changes:

| Document | Refresh when... | Staleness risk |
|---|---|---|
| `README.md` (status line, pillar list, layout tree, doc table) | any module lands, a milestone completes, a live demo changes | high — it is the front door |
| `docs/roadmap.md` | every issue closes (check the box, record measured results honestly, note deferrals with issue numbers) | high |
| `docs/architecture.md` | module boundaries, dataflow, truth-model fidelity, or conventions change | medium |
| `docs/verification.md` | static-subset rules, contract semantics, or the verification pipeline change | medium |
| `docs/comparative-analysis.md` | periodically (quarterly): competitor releases, ecosystem-gap claims ("no other open implementation") must be re-verified before repeating them | medium, time-driven |
| `docs/trajectory-optimization.md` | guidance-layer capabilities or solver choices change | medium |
| `docs/visualization.md` | viewer architecture changes | low |
| `docs/plans/NN-*.md` | written at issue start, acceptance boxes checked at close; historical afterwards (do not retro-edit closed plans except to fix errors) | low after close |
| Package docstrings (`src/podium/*/__init__.py`) | a module in that package lands or changes scope — these are the API-level docs | medium |
| `pyproject.toml` extras + `CONTRIBUTING.md` commands | dependencies or dev workflow change | low |
| `viewer/` demo pages | exported-data schema or control laws they embed change | low |

Rule of thumb for PRs: touch code in `src/podium/X/`, then check the
README layout tree, the roadmap item that covers X, `X/__init__.py`, and
`docs/architecture.md` in that order.

## License

By contributing you agree your contributions are licensed under MIT.
