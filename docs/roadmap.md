# Roadmap

Focus: LEO/MEO RPOD. Each milestone is shippable and validated before the
next starts.

## v0.1 — "CW sandbox" (current)

- [x] CW kernel: dynamics, closed-form STM, two-impulse targeting
- [x] Quaternion kernel (scalar-first, static subset)
- [x] Fixed-step integrators (RK4, Euler)
- [x] Glideslope guidance (Hablani), LQR synthesis/application split
- [x] Contract layer (`podium.verify`): interval contracts + `prove()`
- [x] Tschauner-Hempel / **Yamanaka-Ankersen STM** (no other open
      implementation exists — immediate differentiator); validated against
      nonlinear elliptic relative dynamics up to e = 0.35 and the CW limit
- [ ] Nonlinear relative-motion truth model + J2 + exponential drag;
      quantified CW/TH validity envelopes in CI
- [ ] Sim engine: fixed-step master clock, event detection, scenario config
- [ ] `sim.to_viewer_json()` + matplotlib analysis plots

## v0.2 — "Convex guidance" (Layer 0)

- [ ] Direct transcription on CW/YA STMs: LP (L1 fuel), SOCP (L2), QP
      (tracking) via CVXPY, DPP-parametrized from day one
- [ ] Constraint library: approach cone, rotating-hyperplane KOZ, plume
      half-space, Breger-How passive-safety scenarios
- [ ] LCvx thrust-annulus option with validity checks (Lu & Liu)
- [ ] Clarabel default solver; QOCO alternate
- [ ] **ARCH-COMP rendezvous benchmark** as an executable example with model
      export for reachability tools (CORA/JuliaReach)

## v0.3 — "Full loop"

- [ ] Relative-nav EKF (fixed-dimension, Joseph form, static subset)
- [ ] Sensor models: relative GNSS, docking camera, lidar; actuator MIB/rise
- [ ] Attitude dynamics + quaternion-feedback controller; thruster allocation
- [ ] Monte Carlo campaigns (seeded, structured-array output)
- [ ] three.js interactive viewer (fermi patterns; docs/visualization.md)

## v0.4 — "SCP docking" (Layer 1)

- [ ] PTR successive convexification for 6-DOF approach+docking; SCvx* mode
- [ ] State-triggered constraints (plume, min-impulse-bit); CTCS augmentation
- [ ] Evaluate OpenSCvx as the SCP core vs. in-house loop
- [ ] Contact/capture via MuJoCo backend; capture-envelope MC analysis

## v0.5 — "Flight path"

- [ ] C emitter for the static subset + contract→annotation rendering
- [ ] Golden-vector Python↔C equivalence harness
- [ ] External abstract-interpretation validation gate in CI
- [ ] CVXPYgen→QOCOGEN embedded solver generation for Layer-0 problems
- [ ] cFS/F´ integration example (generated GNC app on a software bus)

## Cross-cutting, every release

- Cross-validation oracles: Orekit (orekit-jpype) and/or GMAT propagation
  comparisons in CI; tudatpy for 6-DOF once attitude lands
- Determinism tests (bit-identical replay)
- No new dependency unless license-vetted permissive (see
  docs/comparative-analysis.md)
