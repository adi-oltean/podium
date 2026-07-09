# Contributing

Podium is a verification-first library: its value is that the flight-relevant code
carries machine-checkable evidence. Contributions are welcome — bug reports,
dynamics/GNC implementations, validation cases, documentation — but the bar is
deliberately high. A change that weakens the verification story will not be merged,
however useful the feature. Read this before opening a PR.

## Non-negotiables

Hard requirements. A PR that violates any of them is rejected or asked to move the
offending code out of the trusted path.

1. **The trusted checkers stay float-free.** Everything under `src/podium/verify/`
   runs in exact `fractions.Fraction` arithmetic — no `float`, no NumPy, no tolerance
   comparisons in the trusted path. Solver output is rationalized *before* it reaches
   the checker, which then verifies exactly or **refuses**. A checker must never
   accept on a tolerance.
2. **`src/podium/core` obeys the static subset**
   ([`docs/verification.md`](docs/verification.md)): pure step functions, statically
   shaped arrays, compile-time-bounded loops, no dynamic allocation, no recursion,
   machine-readable contracts on every input. Code that cannot meet this belongs in
   the sandbox layer, not `core`.
3. **Determinism is absolute** in `sim` and the audit path: no wall-clock, no unseeded
   randomness, no platform-dependent reductions. Same source + same seed → same bytes.
4. **Physics and optimality claims need receipts.** New dynamics, guidance, or
   certificates ship with a test against an independent reference: a closed-form
   solution, a cited published example, cross-validation against the nonlinear truth
   model with quantified error, or an exact certificate. "It looks right" is not
   evidence.
5. **Emitted kernels stay verifiable.** Changes to `emit/` or to an emitted kernel
   must keep the C in the analyzable subset, pass Frama-C/EVA alarm-free, and match
   the Python reference on the golden vectors — bit-exact for the scalar/sqrt class,
   documented tolerance classes otherwise (see
   [`docs/numerical-reproducibility.md`](docs/numerical-reproducibility.md)).
6. **Conventions are fixed, repo-wide:** SI units; LVLH x-radial / y-along-track /
   z-cross-track; scalar-first quaternions. No per-module exceptions.

## The gate — every PR, zero warnings

```bash
pip install -e ".[dev]"
ruff check src tests examples   # zero findings
mypy                            # clean
pytest                          # green, coverage may not drop
```

All must pass with **zero** warnings or errors — not "mostly." CI enforces `ruff`
and the full `pytest` suite on every push; a PR does **not** merge over red CI, and CI
is never bypassed with `--no-verify` on shared branches. New code carries tests;
kernel math carries property-based (`hypothesis`) tests and, for emitted kernels,
golden vectors.

Multi-orbit propagations are marked `slow`; `pytest -m "not slow"` is the fast local
lane (<5 s), but the full suite must pass before you push.

## Dependencies

The core imports only NumPy; the trusted checkers import only the standard library.
Optional stacks (cvxpy, clarabel, ecos, mujoco, Julia/JVM tooling) are
**lazy-imported** so the trusted path stays importable and float-free without them.
Adding a runtime dependency to `core`/`verify` requires explicit justification and
maintainer sign-off; otherwise put it behind a lazy import in the sandbox layer.

## Scope discipline

Truth models may use anything in SciPy; flight algorithms live in `podium.core` under
the subset rules and are exercised through the same step-function interface they will
have in C. Keep functions small — the verification story depends on it. Do not widen a
PR beyond its stated purpose.

## Documentation upkeep

A change is not done until the docs that describe it are current:

| Document | Refresh when… |
|---|---|
| `README.md` (status, pillar list, layout tree, doc table) | any module lands, a demo changes, a count moves |
| `docs/architecture.md` | module boundaries, dataflow, fidelity, or conventions change |
| `docs/verification.md` | subset rules, contract semantics, or the pipeline change |
| `docs/numerical-reproducibility.md` | golden-vector methodology, equality classes, or cross-ISA conditions change |
| `docs/comparative-analysis.md` | quarterly — re-verify ecosystem-gap claims before repeating them |
| Package docstrings (`src/podium/*/__init__.py`) | a module in that package lands or changes scope |
| `pyproject.toml` extras + `CONTRIBUTING.md` | dependencies or the dev workflow change |

Rule of thumb: touch code in `src/podium/X/`, then refresh the README layout/doc
table, `X/__init__.py`, and `docs/architecture.md`.

## Commits

Focused commits, clear messages, and sign-off (`git commit -s`, Developer Certificate
of Origin). By contributing you agree your contributions are licensed under MIT.
