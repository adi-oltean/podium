# Roadmap

Focus: LEO/MEO RPOD. Each milestone is shippable and validated before the
next starts. (Revised 2026-07-03 after a literature/practice sweep; the
headline addition is the ROE module, #5.)

## v0.1 — "CW sandbox" (current)

- [x] CW kernel: dynamics, closed-form STM, two-impulse targeting
- [x] Quaternion kernel (scalar-first, static subset)
- [x] Fixed-step integrators (RK4, Euler)
- [x] Glideslope guidance (Hablani), LQR synthesis/application split
- [x] Contract layer (`podium.verify`): interval contracts + `prove()`
- [x] Tschauner-Hempel / **Yamanaka-Ankersen STM** (no other open
      implementation exists — immediate differentiator); validated against
      nonlinear elliptic relative dynamics up to e = 0.35 and the CW limit
- [x] Nonlinear relative-motion truth model (dual ECI propagation, exact
      LVLH transform incl. frame-precession term) + J2 + exponential drag
      with differential ballistic coefficients (#1)
- [x] Quantified CW/YA validity envelopes in CI: YA position error
      < C(e)·sep²/a per orbit with C ≈ 40 (e ≤ 0.05) / 200 (e = 0.2),
      quadratic scaling asserted; CW-vs-YA degradation ratios documented
- [x] Sim engine (#7): fixed-step master clock with truth substeps,
      scenario config, seeded measurement noise (bit-identical replay
      enforced by test), impulsive flight-block interface v0, burn log,
      channel extraction, threshold-crossing events; closed-loop
      glideslope and CW-synthesized LQR flown against the nonlinear truth
      through the engine. Continuous-thrust interface and bisection event
      refinement deferred
- [ ] Stochastic atmosphere option for the truth model: seeded
      mean-reverting density perturbation on the exponential baseline,
      envelope calibrated to observed storm excursions (+50–125% at
      200–400 km, May 2024 class events) — deterministic replay preserved
- [x] v0 web viewer on GitHub Pages (#2): canvas playback of the V-bar
      approach, follow camera, burn timeline, log-range scrubber
- [x] `Trace.to_viewer_json()` export API (schema-compatible with the live
      viewer)
- [ ] matplotlib analysis plots

## v0.2 — "Relative motion complete + convex guidance" (Layer 0)

- [x] **ROE module (#5)**: quasi-nonsingular relative orbital elements,
      ROE↔elements and ROE↔LVLH maps, Koenig closed-form Keplerian and J2
      STMs (every J2 entry pinned by a finite-difference Jacobian of the
      exact secular flow), impulsive control-input matrix (validated by
      finite-difference impulses in the ECI truth model). Follow-up (#6):
      J2+differential-drag STM variant, eccentric-valid LVLH map,
      analytic minimum-separation ellipse screening
- [x] **ROE-native passive safety** (`podium.guidance.safety`): e/i-vector
      separation angle and RN-plane minimum-separation screening
      (along-track-drift-independent lower bound), validated against
      brute-force scans; Breger-How Cartesian constraints remain the
      planned terminal-phase complement
- [ ] Direct transcription on CW/YA/ROE STMs: LP (L1 fuel), SOCP (L2), QP
      (tracking) via CVXPY, DPP-parametrized from day one
- [ ] Constraint library: approach cone, rotating-hyperplane KOZ, plume
      half-space, passive-safety scenarios
- [ ] LCvx thrust-annulus option with validity checks — including the
      discrete-time lossless-convexification conditions (continuous-time
      guarantees do not survive discretization unconditionally)
- [ ] Clarabel default solver; QOCO alternate
- [ ] **ARCH-COMP rendezvous benchmark** as an executable example with
      model export for reachability tools (CORA/JuliaReach), designed for
      later use as a CI regression gate
- [x] Spec registry v0 (#7, `podium.sim.spec`): named requirements with
      STL robust semantics over trace channels (PUS-12-shaped base
      fragment: always/eventually/final with windows), margins as pytest
      oracles, co-designed with the engine. Follow-ups: rtamt backend for
      full STL, robustness-guided falsification lane in CI

## v0.3 — "Full loop"

- [ ] Relative-nav EKF (fixed-dimension, Joseph form, static subset)
- [ ] Sensor models: relative GNSS, docking camera, lidar; actuator MIB/rise
- [ ] Attitude dynamics + quaternion-feedback controller; thruster
      allocation with explicit minimum-impulse-bit handling
- [ ] Docking acceptance tests against the IDSS IDD Rev G contact-condition
      box (closing 0.05–0.10 m/s, lateral rate 0.04 m/s, angular rates
      0.20 deg/s, lateral offset 0.10 m, angular misalignment 4 deg)
- [ ] Monte Carlo campaigns (seeded, structured-array output)
- [ ] three.js interactive viewer (fermi patterns; docs/visualization.md)

## v0.4 — "SCP docking" (Layer 1)

- [ ] PTR successive convexification for 6-DOF approach+docking; SCvx* mode
- [ ] Continuous-time constraint satisfaction (CTCS) so corridor/KOZ hold
      between nodes, not only at them; state-triggered constraints (plume,
      min-impulse-bit)
- [ ] Temporal-logic mission constraints in the SCP stack via smooth
      robustness encodings, with mixed-integer reference solutions as
      offline validation
- [ ] Evaluate OpenSCvx as the SCP core vs. in-house loop
- [ ] Contact/capture via MuJoCo backend; capture-envelope MC analysis
- [ ] Tumbling-target terminal guidance (rotating corridor, variable-
      horizon endpoint) — scoped study

## v0.5 — "Flight path"

- [ ] C emitter for the static subset + contract→annotation rendering
      (analyzer annotations and ACSL), binary64 throughout, generated C
      kept within the CompCert-compilable subset
- [ ] Golden-vector Python↔C equivalence harness, two tiers: bit-exact on
      host (pinned FP semantics, evaluation-order-matched references,
      correctly-rounded transcendentals) and ULP-bounded on target;
      differential fuzzing from contract ranges
- [ ] Open abstract-interpretation gate in CI: sound float-interval
      analysis as the primary value gate plus a memory/index gate;
      reproducible audit evidence
- [ ] CVXPYgen→QOCOGEN embedded solver generation for Layer-0 problems
- [ ] cFS/F´ integration example (generated GNC app on a software bus)

## Cross-cutting, every release

- Cross-validation oracles: Orekit (orekit-jpype) and/or GMAT propagation
  comparisons in CI; tudatpy for 6-DOF once attitude lands
- Determinism tests (bit-identical replay)
- Truth-model credibility documented along the NASA-STD-7009B assessment
  dimensions (verification, validation, input pedigree, uncertainty)
- Generated-code workflow aligned with NPR 7150.2D SWE-146 and
  ECSS-E-ST-40C Rev.1 (2025) expectations for autogenerated software
- No new dependency unless license-vetted permissive (see
  docs/comparative-analysis.md)
