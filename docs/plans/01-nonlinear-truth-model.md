# 01 — Nonlinear truth model with J2 and drag

GitHub issue: https://github.com/adi-oltean/rpod-lib/issues/1

## Problem

Podium's linearized guidance models (CW, Yamanaka-Ankersen) currently have
no in-library truth reference. The nonlinear elliptic relative dynamics used
to validate the YA STM live inside `tests/test_ya.py` and cover two-body
only — no J2, no drag, no reusable module, and no quantified statement of
where the linear models stop being valid (roadmap v0.1 items 2 and CI
validity envelopes).

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/constants.py` | MISS | SI constants (mu, Re, J2, omega_earth) |
| `src/podium/dynamics/nonlinear.py` | MISS | the truth model |
| `tests/test_nonlinear.py` | MISS | physics validation |
| `tests/test_validity_envelopes.py` | MISS | CW/YA error envelopes in CI |
| `docs/roadmap.md` | PARTIAL | check off truth-model item |

## Fix

New sandbox-side module `podium.dynamics.nonlinear` (full Python allowed —
this is the truth layer, not flight code):

- **Formulation:** integrate *both* spacecraft in ECI with identical force
  models and difference into the target LVLH frame. This is exact (no
  relative-dynamics approximations), extends to any force model, and the
  float64 differencing error (~1e-9 m at LEO radii) is far below truth-model
  needs.
- **Forces:** two-body + J2 (WGS84) + exponential-atmosphere drag with a
  co-rotating atmosphere and per-spacecraft ballistic coefficients
  (differential drag matters for dissimilar chaser/target).
- **LVLH transform:** rotation from (r, v); frame angular velocity
  `omega = [(r/h)(a_pert . z_hat), 0, h/r^2]` — the x-component (orbit-plane
  precession under out-of-plane perturbation) is included exactly, not
  dropped.
- `elements_to_rv` for scenario setup; `propagate_relative` convenience
  wrapper (fixed-step RK4, deterministic).

Validation tests (physics receipts):

1. Two-body relative propagation through the full ECI pipeline matches the
   YA STM to linearization error (cross-checks transforms end-to-end).
2. LVLH velocity consistency: central-difference of the position history
   equals the reported relative velocity with J2+drag active (validates the
   omega term including its x-component).
3. J2 secular RAAN drift over 10 orbits matches the analytic rate
   -1.5 n J2 (Re/p)^2 cos(i) within short-period noise.
4. Drag: circular-orbit semi-major-axis decay matches the analytic
   per-orbit estimate; differential BC produces the expected along-track
   drift sign.
5. Two-body energy conservation bounds RK4 truncation error.

Validity envelopes (CI): for separations {0.1, 1, 10} km, e in {0, 0.05,
0.2}, one orbit: assert YA error stays below a quadratic-in-separation
bound; assert CW at e=0 matches YA; document CW degradation with e.

## Tests

New: `tests/test_nonlinear.py`, `tests/test_validity_envelopes.py`.
Existing `tests/test_ya.py` stays (independent implementation — the two
truth models now cross-check each other).

## Acceptance Criteria

- [x] Code change implemented
- [x] Tests pass (pytest, ruff, mypy strict) — 48 tests green
- [x] Roadmap updated

## Push/merge instructions

Single commit on `main`: `01 — Nonlinear truth model with J2 and drag (#1)`,
push to origin, close #1 with a comment linking the commit.

## Verification steps

`./.venv/bin/pytest` green; envelope test prints/asserts the documented
bounds; example still runs.
