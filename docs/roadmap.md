# Roadmap

Focus: LEO/MEO RPOD. Each milestone is shippable and validated before the
next starts. (Revised 2026-07-03 after a literature/practice sweep; the
headline addition is the ROE module, #5.)

## v0.1 — "CW sandbox" (COMPLETE, 2026-07-03)

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
- [x] Stochastic atmosphere option (#13, `DensityPerturbation`): seeded
      mean-reverting (OU) log-density factor on the exponential baseline,
      exactly discretized on a fixed grid (deterministic in time — RK4 +
      bit-identical replay preserved and enforced by test); default
      calibration puts +2σ at ~2.0× inside the +50–125% storm band;
      σ=0 reduces exactly to the baseline; OU statistics pinned
- [x] v0 web viewer on GitHub Pages (#2): canvas playback of the V-bar
      approach, follow camera, burn timeline, log-range scrubber
- [x] `Trace.to_viewer_json()` export API (schema-compatible with the live
      viewer)
- [x] matplotlib analysis plots (#13, `podium.sim.plots`): trajectory
      plane view with burns, channel time-series grid, dv timeline —
      object-API figures (headless/thread-safe), matplotlib optional

## v0.2 — "Relative motion complete + convex guidance" (Layer 0, current)

- [x] **ROE module (#5)**: quasi-nonsingular relative orbital elements,
      ROE↔elements and ROE↔LVLH maps, Koenig closed-form Keplerian and J2
      STMs (every J2 entry pinned by a finite-difference Jacobian of the
      exact secular flow), impulsive control-input matrix (validated by
      finite-difference impulses in the ECI truth model)
- [x] ROE follow-ups (#6): density-model-free **J2+drag STM** (augmented
      7-state with constant relative sma-decay rate; pinned by an FD
      Jacobian of the exact augmented flow with closed-form rate
      integrals, and validated against the differential-BC truth over 8
      orbits — quadratic mean-longitude runaway within 25%);
      **eccentric-valid ROE→LVLH map** (Jacobian of the exact nonlinear
      chain, sandbox side; first-order-exact at e=0.3 by quadratic-
      scaling test); **analytic min-RN-separation** (quartic stationary
      points on the unit circle, matches 200k-point scans to 1e-4
      relative; the bounded scan remains the static-subset variant)
- [x] **ROE-native passive safety** (`podium.guidance.safety`): e/i-vector
      separation angle and RN-plane minimum-separation screening
      (along-track-drift-independent lower bound), validated against
      brute-force scans; Breger-How Cartesian constraints remain the
      planned terminal-phase complement
- [x] Direct transcription on exact STMs (#8): `RendezvousPlanner`
      (CW/YA, Cartesian) and `RoePlanner` (Keplerian/J2 ROE dynamics +
      control matrix), L1/LP and L2/SOCP fuel objectives, per-burn caps,
      DPP-compiled once with STMs/boundaries/normals as Parameters;
      QP tracking objective deferred
- [x] Constraint library v0 (#8): approach cone (exact SOC),
      rotating-hyperplane KOZ (bounded two-pass reference refinement;
      hyperplane implies true distance), plume half-space (arrival burn
      exempt by design, active constraints recorded on the Plan)
- [x] Breger-How passive-safety scenarios (#9): per-failure-node
      free-drift KOZ avoidance as linear constraints (normal folded with
      the drift STM into one DPP parameter row per scenario sample);
      convex e/i safe-set terminal for the ROE planner (alignment cones +
      minimum magnitudes, exact rn_margin scan as the post-solve receipt);
      QP tracking objective for MPC-style re-solves
- [x] LCvx thrust-annulus (#9, `FiniteBurnPlanner`): classical slack
      relaxation on the exact ZOH CW discretization, shipped with the
      validity checks rather than the assumption — controllability
      precondition plus a per-solve losslessness audit against the
      discrete-time theory bound (non-tight nodes <= state dimension);
      the audit provably flags degenerate excess-capacity problems where
      the relaxation goes loose (tested both ways). YA/ROE finite-burn
      dynamics deferred (#12)
- [x] Clarabel default solver (QOCO alternate deferred to the embedded
      path)
- [x] **ARCH-COMP rendezvous benchmark as a CI reachability gate** (#10):
      executable Podium model (`guidance.arch`, abort mode verified to be
      planar CW at GEO mean motion against `core.cw`), machine-readable
      hybrid-automaton export, JuliaReach proof of the LOS-cone, velocity-
      octagon, and abort-avoidance properties (SRNA01 + SRA01, ~14 s reach
      time, PROVEN locally and gated in `.github/workflows/reach.yml` on
      guidance/control/dynamics changes + weekly). Follow-up: STL-property
      checking (CORA lane)
- [x] **Podium-synthesized gains through the gate** (#11): continuous
      LQR via a CARE solver (`control.lqr.care`/`clqr`, Hamiltonian
      stable-subspace method, residual pinned at machine precision);
      switched-controller gains synthesized on the `core.cw` planar plant
      and PROVEN by JuliaReach for both scenarios alongside the reference
      controller — the full synthesize→export→prove workflow, gated in CI
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
