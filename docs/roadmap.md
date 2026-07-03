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
      the relaxation goes loose (tested both ways)
- [x] Layer-0 residuals (#12): **eccentric finite-burn** (YA ZOH
      discretization with composition-identity receipts; per-interval
      time-varying maps and lifted-reachability controllability check);
      **primer normality certificate** from the solver duals (primer/dt
      is scale-free: O(0.1) normal vs ~1e-7 degenerate, with the
      Gamma-stationarity identity primer==dt at interior-slack nodes
      pinning the dual convention); **min-time pre-solve** (bounded
      bisection; LCvx at 1.15x t_min passes the audit — the coast-arc
      antidote); **dense passive-safety verification margins** reported
      on every plan (exposes inter-sample dips instead of hiding them);
      **MIB bridge** (`quantize_plan`: sigma-delta thruster-click
      quantization, residual <= half a click per axis, flown through the
      engine against the nonlinear truth)
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

## v0.3 — "Full loop" (COMPLETE, 2026-07-03)

- [x] Relative-nav EKF (#14, `podium.nav.ekf`): fixed 6-state, Joseph-form
      update (symmetry/PD preserved under roundoff, tested over 500-step
      random sequences), white-noise-acceleration process model, CW STM
      prediction with commanded burns fed through as known inputs.
      Receipts: NIS within the chi-square band against the engine's
      seeded truth; convergence from 100 m initial error; LQR flying on
      EKF estimates from 5 m position-only measurements stabilizes the
      nonlinear truth closed-loop
- [x] Sensor models + actuator imperfections (#15, `podium.nav.sensors`):
      relative GNSS (white noise + seeded constant bias), docking camera
      and lidar (az/el bearing + proportional range noise, visibility-
      gated; Jacobian pinned by FD; EKF gains a nonlinear Joseph update
      with bearing wrap-around); engine actuator model (MIB quantization,
      per-tick cap, seeded proportional execution error — receipts show
      open-loop pulse plans miss by hundreds of meters under 2% execution
      error while feedback through the same actuator converges).
      Thruster rise/tail-off shaping deferred to the 6-DOF layer (v0.4)
- [x] Attitude dynamics + quaternion-feedback controller (#16):
      rigid-body Euler equations + quaternion kinematics (RK4,
      renormalized; torque-free energy and inertial angular momentum
      conserved to 1e-9 through an intermediate-axis tumble);
      shortest-way quaternion-feedback regulator with saturation —
      detumble receipt and a 20-degree slew pinned to the second-order
      design prediction (overshoot and settling of wn=0.1, zeta=0.9).
      Thruster torque allocation deferred to the 6-DOF layer (v0.4)
- [x] Docking acceptance vs the IDSS IDD Rev G box (#17,
      `podium.sim.idss`): terminal rate-command approach flown through
      the engine with sensor noise + MIB + execution error contacts
      inside the translation box (closing 0.05–0.10 m/s, lateral rate
      0.04 m/s, offset 0.10 m); quaternion-feedback hold satisfies the
      rotational box (0.20 deg/s, 4 deg) at contact time. Translation
      and rotation decoupled until the 6-DOF engine (v0.4), coupled at
      the contact instant and documented as such
- [x] Monte Carlo campaigns (#17, `podium.sim.monte_carlo`): master-seed
      spawned per-run seeds, structured-array output with per-run seeds
      for post-mortem replay; 20-run dispersed docking campaign 100%
      in-box with bit-identical campaign reproducibility asserted
- [x] three.js interactive viewer (#18, `viewer/3d/`): zero-build page
      with vendored three.js r172 (MIT); LVLH scene (station + docking
      ring, KOZ wireframe, approach-corridor cone, V-bar grid),
      preallocated-buffer trail, burn glyphs, follow + free-orbit
      cameras, play/scrub/HUD sharing the 2-D viewer's DATA schema;
      Playwright receipts (WebGL renders non-blank, playback/scrub,
      physics end-state, zero console errors). Deferred: ECI/LVLH frame
      blending (needs target-ECI export), log-time map, attitude
      keyframes

## v0.4 — "SCP docking" (Layer 1, current)

- [x] PTR/SCvx* SCP core (#21, `podium.guidance.scp`): penalized trust
      region over the exact-STM transcription; TRUE nonconvex keep-out
      sphere (re-linearized per iteration), virtual buffers, SCvx*-style
      penalty ramp on infeasibility stall, trust-region expansion on
      validated feasible iterates, flat-valley-aware convergence (fuel
      stationarity). Receipts: reduces to Layer-0 on convex problems;
      passage problems converged with virtual buffers at zero and cost
      within 0.1% of the hyperplane heuristic while satisfying the true
      sphere; penalty ramp from 1e-6 demonstrated; deterministic; plan
      flown against the nonlinear truth. 6-DOF (attitude-coupled) PTR
      remains open with the contact layer
- [x] CTCS, exact-flow form (#21): coast arcs follow the exact STM flow,
      so intermediate-time positions are LINEAR in the decision
      variables — continuous-time KOZ violations become exact linear
      cuts (times persistent, directions re-linearized each iteration;
      stale fixed cuts provably wedge the loop and are avoided by
      design), iterated to a clean independent 1000-sample dense check.
      Receipt: a coarse grid whose coast dips inside the sphere between
      nodes is caught and cut. Integral-augmentation CTCS (for future
      non-coast dynamics) and state-triggered constraints stay open
- [ ] 6-DOF attitude-coupled PTR (body-frame thrust, contact attitude,
      angular corridors) — the consciously deferred remainder of the
      PTR and contact items above; pairs with the thruster torque
      allocation deferred from v0.3
- [x] Temporal-logic mission constraints in SCP (#22,
      `EventuallyBoxSpec`): timed-window reach specs via smooth
      robustness with a SOUND encoding split — node box margins enter
      the subproblem exactly as hypograph variables (they are concave),
      and only the convex LSE smooth-max is tangent-linearized (a convex
      function dominates its tangent, so tangent >= eps + ln(K)/tau
      implies true robustness >= eps; linearizing the concave margins
      directly is a relaxation the optimizer provably exploits —
      observed as a period-2 oscillation before the fix). Receipts: spec
      bites (negative robustness without it), converged plans meet TRUE
      non-smooth robustness, KOZ stays clean, and the engine-flown trace
      passes a spec-registry eventually-check — guidance, truth, and
      monitoring agreeing on one temporal property. Deferred: MIP
      reference validation (no MIP dependency), richer fragments
      (until/nested), robustness-maximization objective
- [x] Evaluate OpenSCvx vs in-house (#21, decision recorded): in-house —
      the Layer-0 cvxpy transcription/receipt infrastructure carries
      directly, problem sizes are tiny, and the exact-flow cut mechanism
      is specific to our coast-arc structure. Revisit at 6-DOF where
      OpenSCvx's discretization machinery earns its keep
- [x] Contact/capture via MuJoCo (#24, `podium.sim.contact`):
      probe-drogue model generated programmatically as MJCF (funnel from
      convex plates + throat sleeve — box-edge contact normals stall
      probes, measured and designed out; ~29-deg cone after the 47-deg
      first cut reflected too much axial momentum; optional sustained
      docking thrust per probe-drogue practice). Receipts: IDSS-box
      corner conditions ALL physically capture (the acceptance box tied
      to a capture mechanism), beyond-mouth misses, 2 m/s bounces, peak
      force monotone in closing rate, deterministic, envelope boundary
      bracketed. Deferred: contact attitude (6-DOF), latch model,
      compliant drogue, dispersed MC campaign via monte_carlo
- [x] Tumbling-target terminal guidance — scoped study (#23,
      `podium.guidance.tumbling` + docs/plans/23-tumbling-study.md):
      known-tumble port capture stays CONVEX (deterministic port
      kinematics → terminal boundary state + per-node rotating-corridor
      cones on the exact CW STM). Findings pinned by tests: fuel-vs-rate
      is phase-confounded (fix the arrival phase for clean envelopes);
      the naive co-rotation cost intuition fails at low rates (matching
      a slow port beats nulling all motion — CW drift supplies velocity
      free); measured envelope closes at ~3 deg/s under a 0.35 m/s burn
      cap (~1 deg/s flyable at 2.1 m/s). Engine flight arrives on the
      independently-recomputed rotating port within 1 m / 1 cm/s.
      Follow-on: uncertain tumble (estimator + robust corridors), 3-D
      nutation, plume-vs-rotating-body — where SCP re-enters
- [x] **Infinite-horizon abort-safety certificates** (#20,
      `podium.verify.barrier`): barrier functions over the CW flow
      invariants in time-scaled coordinates (integer dynamics matrix →
      fully rational problem), synthesized by an untrusted SDP
      (cvxpy/Clarabel S-procedure) and re-verified by a TRUSTED checker
      in exact fractions arithmetic (all-principal-minors PSD test, Lie
      derivative exactly zero — no floats in the trusted path). Certifies
      "passive abort drift keeps RN separation outside the keep-out
      radius for ALL time", along-track-independent (the machine-checked
      e/i-separation heritage argument). Receipts: end-to-end synth →
      rationalize → exact verify; hand-derived certificate whose algebra
      closes on paper; tamper detection; V-bar hold correctly infeasible
      (not passively safe); dense-propagation corroboration. Follow-ups:
      quartic barriers for tighter sets, J2-perturbed variant, attitude
      closed loop

## v0.5 — "Flight path"

- [x] C emitter v0 (#25, `podium.emit.cemit`): AST-based C99 emitter for
      the static subset — the supported-subset checker REJECTS anything
      outside it, making the emitter the operational StaticPy
      definition. Covers scalar/fixed-array kernels, constant-index
      subscripts (1-D/2-D), whitelisted math.*, branches, cross-kernel
      calls lowered through temporaries; contracts render as ACSL
      requires clauses + analyzer [spec] blocks. First kernels: the
      quaternion family + CW (mean_motion, cw_deriv, stm). Remaining:
      bounded-for-loop support to cover roe/ya/integrators/EKF (v0.6),
      CompCert-subset audit
- [x] Golden vectors, tier 1 (#25): emitted C under pinned FP semantics
      (-O2 -ffp-contract=off, SSE2 binary64) reproduces Python
      BIT-FOR-BIT over 2000 seeded vectors per kernel including branch
      paths — for every arithmetic+sqrt kernel (sqrt is IEEE
      correctly-rounded). Measured exception that motivates the
      CORE-MATH item: sin/cos differ between this interpreter's libm
      and system glibc on ~0.03% of stm values, bounded ≤4 ulp after
      propagation — asserted as such, not hidden. Tier 2 (ULP-bounded
      on target) open
- [ ] Open abstract-interpretation gate in CI: sound float-interval
      analysis as the primary value gate plus a memory/index gate;
      reproducible audit evidence
- [x] CVXPYgen embedded solver generation (#26,
      `podium.emit.solvergen`): fixed-grid Layer-0 rendezvous with live
      boundary parameters generated to a self-contained C tree (ECOS
      backend), compiled with plain gcc (build recipe encoded: gnu99,
      -fcommon, demo-source exclusions, SuiteSparse includes) and run
      with zero Python — the binary reproduces the Clarabel optimum to
      1e-5. cvxpygen deliberately NOT a dev dependency (its import pulls
      a Julia sidecar via pdaqp; generation is a local/offline step).
      QOCOGEN alternate + verified-KKT checker are the v0.6 items
- [ ] cFS/F´ integration example (generated GNC app on a software bus)

## v0.6 — "Certified reference mission"

The layers exist; v0.6 composes them into one auditable whole.

- [x] Emitter v1 (#27): bounded `for range(N)` loops (compile-time
      bounds, module-constant resolution), module-constant inlining,
      tuple/augmented assignment, loop-var subscripts, np.eye
      allocation. 17 kernels now emit + verify, incl. the bounded
      Newton Kepler solve and the first CONTRACTED kernels through the
      ACSL path (roe maps/STM). Tier-1 policy recalibrated from
      measurement: strict bit-exact for arithmetic+sqrt kernels;
      trig-bearing kernels bounded by output-vector scale (cancelling
      entries diverge by ulps of their INTERMEDIATES — measured and
      documented), <=1% incidence. Remaining for the full core: matmul
      lowering (EKF), integrators (function-typed params), CompCert
      audit
- [ ] Sound value gate in CI: Frama-C/EVA over the emitted-and-annotated
      C (float intervals from the ACSL contracts), memory/index gate,
      reproducible audit report artifacts
- [ ] Correctly-rounded transcendentals option (CORE-MATH) closing the
      measured tier-1 sin/cos gap; tier-2 ULP-bounded golden vectors on
      a cross-compiled target (qemu-aarch64)
- [ ] CVXPYgen/QOCOGEN embedded generation of a Layer-0 problem with the
      verified-KKT-checker pattern (certificate checked by exact/interval
      arithmetic, R4-style)
- [ ] 6-DOF attitude-coupled PTR + contact attitude (carried from v0.4);
      thruster torque allocation
- [ ] End-to-end reference mission: far-range ROE phasing → corridor
      approach (SCP plan, EKF nav, imperfect actuators) → IDSS contact →
      MuJoCo capture, as one seeded scenario with a release-grade audit
      bundle (spec margins, MC table, reach verdicts, barrier
      certificate, golden-vector attestation) published per tag
- [ ] Orekit cross-validation lane in CI (orekit-jpype) for the truth
      model; three.js viewer frame-blending using the target-ECI export
- [ ] cFS or F´ integration example: the generated GNC app on a software
      bus, fed by the reference-mission scenario

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
