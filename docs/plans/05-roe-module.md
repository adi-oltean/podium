# 05 — ROE module: quasi-nonsingular state, Koenig STMs, e/i-vector passive safety

GitHub issue: https://github.com/adi-oltean/podium/issues/5

## Problem

Podium's relative-motion layer is Cartesian-only (CW + Yamanaka-Ankersen).
For perturbed, multi-orbit, km-scale relative motion — and for the
flight-proven passive-safety formulation — the field standard is relative
orbital elements: secular J2/drag effects are absorbed analytically,
linearization error is lower at large separations, and passive safety
becomes an algebraic phase condition on the relative e/i vectors rather
than a propagation screen. No maintained open-source implementation of the
Koenig closed-form STMs exists (verified July 2026), so this module is
both a capability gap and a differentiator, mirroring the YA situation.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/core/roe.py` | MISS | state defs, maps, STMs, control-input matrix |
| `src/podium/guidance/safety.py` | MISS | e/i-separation + ellipse screening checks |
| `tests/test_roe.py` | MISS | truth-model receipts + structural identities |
| `docs/architecture.md` | PARTIAL | add ROE to the dynamics section |
| `docs/roadmap.md` | DONE | v0.2 items added 2026-07-03 |

## Fix

`podium.core.roe` (static subset, contracts on every scalar):

1. Quasi-nonsingular ROE state (da, dλ, dex, dey, dix, diy) with
   conversions: osculating/mean orbital elements ↔ ROE; ROE ↔ LVLH
   Cartesian (first-order maps, near-circular bridge per Gaias-Lovera for
   cross-validation against the CW layer).
2. Closed-form STMs (Koenig, Guffanti & D'Amico, JGCD 40(7) 2017,
   doi:10.2514/1.G002409): Keplerian, J2, and J2 + differential-drag
   (density-model-free variant first), arbitrary eccentricity, fixed-shape
   6x6 matrices — bounded trigonometry only, no iteration.
3. Impulsive control-input matrix (delta-v in RTN → delta-ROE).
4. `podium.guidance.safety`: relative e/i-vector separation metric and
   minimum-separation ellipse checks (AAS 23-155 formulation) as pure
   functions returning margins, usable as guidance constraints and sim
   monitors.

Validation (physics receipts, per house rules):

- STM vs truth: propagate the nonlinear ECI truth model (J2 on/off,
  differential BC on/off) and compare ROE STM predictions over 1-30
  orbits; tolerance derived from the linearization error scaling in
  Koenig Table/analysis, with the quadratic-scaling discriminator test.
- Cross-model: ROE↔LVLH maps against CW/YA propagation in the
  near-circular unperturbed limit (Gaias-Lovera equivalence).
- Structural: composition, identity at dt=0, invertibility; J2 secular
  rates in ROE space vs the analytic RAAN/argument drift already tested
  in the truth model.
- Safety checks: e/i separation vs brute-force minimum-distance
  propagation screening over dispersed along-track errors (target the
  published ~99% agreement; document disagreement cases).

## Tests

`tests/test_roe.py` (new; hypothesis for map round-trips), extension of
`tests/test_validity_envelopes.py` with an ROE column (CW vs YA vs ROE
error over separation/eccentricity/duration/perturbation grid).

## Acceptance Criteria

- [ ] Code implemented under static-subset rules with contracts
- [ ] Receipts green (pytest, ruff, mypy) incl. truth-model comparisons
      and quadratic-scaling discriminators
- [ ] Validity-envelope table extended and documented
- [ ] architecture.md updated

## Push/merge instructions

Single commit on main: `05 — ROE module (#5)`, push, close #5 with a
summary of measured envelopes.

## Verification steps

Full suite green; envelope table shows the expected regime ordering
(ROE ≥ YA ≥ CW under J2/drag at multi-orbit horizons); e/i safety checks
agree with propagation screening on the dispersion grid.
