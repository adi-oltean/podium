# 32 — Orekit cross-validation lane for the truth model

GitHub issue: https://github.com/adi-oltean/podium/issues/32

## Fix (landed)

`tests/test_orekit_validation.py` (skips cleanly without JVM/data;
slow-marked) + `.github/workflows/validate.yml` (weekly + dynamics
paths + dispatch; setup-java temurin 21, cached orekit-data.zip).
Local environment: a JRE + orekit-data.zip under `./.orekit-cache`
(or set `JAVA_HOME` / `OREKIT_DATA_ZIP`), orekit-jpype in
the validate extra. numpy note: orekit-jpype pins numpy<2.3 as
metadata; dev venv keeps numpy 2.5 (works fine at runtime, and the
main CI lane never installs orekit-jpype).

## Three receipts, three failure classes

1. Two-body vs Orekit's ANALYTIC Keplerian propagator: <5 cm position
   after 2 orbits (11,200 km of arc) at RK4 dt=5 s — pure integrator
   drift, zero model ambiguity.
2. J2 vs Holmes-Featherstone degree 2: measured residual 1.13 m over
   2 orbits (bound 25 m, >20x margin; a J2 sign error diverges by
   kilometers).
3. Drag DELTAS — (J2+drag) minus (J2) end positions — agree within
   10%: validates drag magnitude, direction, and the co-rotating-
   atmosphere convention (both stacks co-rotate; the drag scenario
   must produce >50 m of effect or the test fails as vacuous).

## The finding worth remembering

The first naive J2 comparison differed by 185 m. Not a bug in either
stack: podium's ECI convention is "inertial, z = CURRENT pole", while
EME2000's z is the J2000 pole — 0.13 deg away in 2026 through
precession — which tilts the J2 symmetry axis (0.0024 rad times ~60 km
of J2 secular displacement = the observed 185 m). Feeding Orekit the
state in the True-of-Date frame (z = pole of date) collapses the
residual to 1.13 m; the remainder is nutation motion during the
interval plus the EGM-vs-IERS J2 in the 8th digit. Documented in the
test docstrings — this is exactly the class of convention error the
lane exists to catch, demonstrated on itself.

## Deferred

Viewer ECI frame-blending (split to its own roadmap line); SPICE
ephemeris comparisons (spiceypy already in the validate extra);
eccentric-orbit and higher-degree gravity cases.

## Push/merge instructions

Single commit on main: `32 — Orekit cross-validation lane (#32)`;
push; validate.yml triggers on its own path; close after green.
