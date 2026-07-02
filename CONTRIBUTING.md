# Contributing

Contributions are welcome — bug reports, dynamics/GNC algorithm
implementations, validation cases against flight data or other simulators,
and documentation all help.

## Ground rules

- **`src/rpod/core` follows the static subset** documented in
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

## License

By contributing you agree your contributions are licensed under Apache-2.0.
