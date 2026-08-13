# TLA+ applicability assessment: Podium's discrete mode logic

Podium's verification stack is extensive on continuous-state
guarantees and deliberately limited on discrete-temporal ones. This
note reads the shipped stack through the safety/liveness taxonomy of
Lamport's *Proving the Correctness of Multiprocess Programs* (IEEE
TSE, 1977), inventories the state machines present in the code, and
identifies where TLA+ — the modern descendant of that proof style —
would add guarantees no current modality expresses. It is an
assessment; nothing below is built. The execution plan derived from
it is [add-tla-specs.md](add-tla-specs.md).

## Lamport's two classes, applied to this stack

Lamport's decomposition: a **safety** property says something bad
never happens and is proved by an invariant — an assertion shown
inductive over every program step; a **liveness** property says
something good eventually happens and needs a different proof shape
(well-foundedness, fairness). The distinction matters because
Podium's proof-backed artifacts and its trace-checked obligations sit
on opposite sides of it.

A third class: **bounded liveness** ("reaches X within T") is
formally a safety property of the clock-augmented state. Podium
already uses exactly this reduction — the ARCH model carries the
clock *as a state variable* "so abort transitions are time-triggered
guards" (`src/podium/guidance/arch.py:13-14`, integrator row
`_B_CLOCK` at `:74`), and every eventually-obligation in the repo
carries a deadline (`t_end=` windows, `max_iters`, the 300-minute
`HORIZON`, the contact-sim duration). The reduction is sound; the gap
is in how the resulting properties are checked.

## Classification of the shipped artifacts

| Artifact | Class | Checked by | Scope |
|---|---|---|---|
| Reachability: LOS cone, velocity octagon, abort keep-out (`tools/reach/arch_rendezvous.jl:69-86`) | Safety | flowpipe inclusion / disjointness | full initial set, symbolic |
| Abort barrier (`podium.verify.barrier`, C1/C2 at `barrier.py:262-273`) | Safety, **infinite horizon** | exact rationals | whole X0 ellipsoid |
| Lyapunov ellipsoid (`verify/lyapunov.py:77-109`) | Safety (invariance) | exact rationals | all x, symbolic |
| SOS / quaternion-cone certificates (`podium.verify.sos`) | Safety | exact rationals | all x, symbolic |
| KKT / optimality-gap certificates (`podium.verify.kkt`, `verify/bracket.py`) | Safety-shaped predicate on one solve | exact rationals | per solver output |
| Contracts → ACSL, EVA zero alarms (`podium.verify.contracts`, `eva.yml`) | Safety (RTE-freedom) | abstract interpretation | contracted box, exhaustive |
| Golden vectors, CompCert, cross-ISA | Safety-shaped equivalence | bit comparison | tested inputs |
| STL `always_*` specs (`sim/spec.py:52-55`), e.g. `koz_far_phase` | Safety | trace monitor | **single trace** |
| STL `eventually_*` specs (`sim/spec.py:56-59`): `reach_hold`, `contact` (`sim/mission.py:217-222`) | Bounded liveness | trace monitor | **single trace** |
| Mission capture + seat dwell (`sim/contact.py:174-187`) | Bounded liveness | MuJoCo rollout | single trace (seed 7; seeds 1–3 in `tests/test_mission.py:51-56`) |
| Solver `status == "converged"` (`guidance/scp.py:328-331`, `sixdof.py:281-291`) | Bounded liveness per instance | test assertions | single instances |
| Closed-loop stability (`tests/test_arch.py:34-37`) | Attractivity (liveness-adjacent) | float eigenvalues | one matrix, no certificate |

**Every artifact backed by a proof is a safety property.** Every
liveness obligation — dock, reach the hold, enter the attempt mode,
converge — is checked on finitely many concrete simulated traces.
Three specific gaps:

- **No liveness operator survives into any proof.** The STL registry
  is a flat five-kind fragment with no `until`, no nesting, no boolean
  combinators (deferred in `docs/plans/22-stl-in-scp.md:45-49`); the
  words *fairness* and *termination* appear nowhere; loop termination
  is by construction (static-subset rule 5), which is sound
  engineering but provides no way to state system-level progress.
- **The reachability lane can pass vacuously.** The Julia gate checks
  LOS and velocity only on flowpipe segments whose location matches
  the attempt mode (`arch_rendezvous.jl:71-76`); if mode 2 were never
  entered, every check passes and CI prints PROVEN. The trace side
  documents the same hazard — modes that never occur report `+inf`,
  "vacuously satisfied, visibly so" (`arch.py:217`) — and the only
  backstop is `assert modes[-1] == 2` / `== 3` on five corner traces
  (`tests/test_arch.py:62,75,128`). Mode entry is a liveness
  property; no current lane can state it.
- **The invariant channel is declared but unpopulated.** `prove(cond,
  label)` (`verify/contracts.py:98-106`) has zero call sites in the
  repo, and its lowering to ACSL obligations is planned, not emitted
  ([verification.md](verification.md)). The exact Lyapunov checker
  certifies only the non-strict decrease `P − A_cl'PA_cl ⪰ 0` —
  invariance, i.e. safety — never attractivity
  (`verify/lyapunov.py:77-84`); asymptotic convergence rests on a
  floating-point eigenvalue test.

## State-machine inventory

None of the machines is an enum. Static-subset rule 6 — "Mode logic
is an explicit enum FSM whose modes and guards match the hybrid model
used for reachability analysis" ([verification.md](verification.md))
— is prescriptive; no code in `podium.core` implements mode logic
yet, and the machines that exist live in the sandbox layers:

| Machine | Where | States | Structure |
|---|---|---|---|
| ARCH rendezvous FSM | `guidance/arch.py:184-197` | int 1/2/3 (`MODE_NAMES` `:76`) | urgent guards, abort dominates; strictly one-way 1→2, 1→3, 2→3; mode 3 absorbing |
| Mission phase machine | `sim/mission.py:127,153-154` | str `"A"` → `"BC"` | one-way time trigger at `t_phase_a` |
| Corridor hold gate | `sim/mission.py:163-164` | closing on/off | stateless, re-evaluated each tick, cycles; keyed on the **estimate** |
| Replan queue / burn index | `sim/mission.py:128,142-152` | consumed-once list, monotone index | monotone by convention, unchecked |
| Capture latch | `sim/contact.py:171-187` | `entered` latch + dwell counter | latch one-way; dwell counter resets on unseat |
| Solver status protocols | `scp.py:277-331`, `sixdof.py:183-291`, `convex.py` | strings (`converged`, `max_iters`, `stalled_feasible`, …) | anytime contract: bounded loop, best iterate always returned |
| Viewer autopilot (illustrative) | `viewer/iss-sim/index.html:177` | waiting → calibrating → active → disengaged | the one FSM with a retry self-loop (calibration) |

Two defects of the spec-code-drift class are present in the ARCH FSM:

- **The attempt-box guard is written twice.** Once as Python booleans
  (`_in_attempt_box`, `arch.py:91-99`) and once as eight exported
  halfspaces (`arch.py:265-270`). Nothing checks that they describe
  the same octagon; a change to one silently diverges from the other.
  The export-schema test compares mode *names*, which are themselves
  re-typed as string literals in `export_model` (`arch.py:279-291`)
  rather than drawn from `MODE_NAMES`.
- **The FSM semantics exist only as a comment.** "transitions
  (urgent, abort dominates)" (`arch.py:187`) plus guard ordering *is*
  the discrete specification; it exists nowhere as a checkable
  artifact.

## Two structural exclusions in the current guarantees

Two deliberate engineering exclusions bound the current guarantees.
Both concern discrete coordination logic — the domain model checking
addresses — so the proposal below is complementary to the existing
lanes rather than overlapping with them.

### Exclusion 1: the verified path assumes sequential execution

There is no concurrency anywhere in the ACSL/Frama-C surface, by
construction. The emitted `kernels.c` contains at file scope only
`#include <math.h>`, prototypes, and pure function definitions
(`emit_module`, `src/podium/emit/cemit.py:891-902`); no emitter code
path can produce a global, `static` state, `volatile`, or a thread,
and the subset tripwire bans goto/switch/union/VLAs/dynamic allocation
(`tests/test_cemit.py:395-406`). The EVA harness is a single entry
point calling each kernel check in sequence (`emit/evagen.py:141-145`),
analyzed under EVA's default *sequential* semantics
(`tools/eva_gate.py:38-41`). This is a precondition of the method:
interval-domain analysis is sound on exactly the Astrée/ATV code
class, one synchronous loop of pure step functions
([verification.md](verification.md)).

The consequence: every guarantee on the verified path — zero alarms,
CompCert semantics preservation, bit-exact golden vectors — is
conditional on the kernels being called sequentially, one loop, no
shared mutable state. That condition holds trivially today, and stops
holding by itself the moment the kernels are embedded in real
multi-rate flight software. cFS runs each app as a task on a shared
software bus, and the example app is already there — a task blocking
on the bus with module state in a `static` struct
(`examples/cfs_nav_app/podium_nav_app.c:19,75-83`) — *outside* the
proof boundary: the EVA gate analyzes only `kernels.c` +
`eva_driver.c` (`tools/eva_gate.py:32-40`). Measurements arriving out
of order, mode commands racing telemetry, multiple apps sharing the
bus: ACSL as used here has no vocabulary for interleavings, and EVA is
not run under a concurrent semantics, so the sequencing assumption
becomes a system-level obligation that no current lane can state.
Proving that the flight-software integration maintains the
sequential-call contract the kernel proofs assume — under all message
orders and task interleavings — is a standard application of TLA+.

### Exclusion 2: event-triggered logic is deliberately outside the verified path

The paper's limitations section is explicit: "variable-horizon,
adaptive, and event-triggered logic is deliberately excluded from this
path and lives in the untrusted planner and controller layers." What
that sentence covers in the code:

- **Variable-horizon** — problem sizes fixed only at run time: each
  replan builds a fresh planner over the remaining time to the phase
  boundary (`sim/mission.py:131`); the static subset requires shapes
  known at import time (rule 2) and compile-time loop bounds (rule 5),
  so this cannot enter the C path.
- **Adaptive** — logic that retunes itself online: the SCvx*-style
  penalty/trust-region updates (`guidance/scp.py:304-310`), the
  accept/reject step with trust-region halving
  (`guidance/sixdof.py:261-278`).
- **Event-triggered** — the decisions about *when*: replan firing
  (`sim/mission.py:142-144`), the phase transition (`:153-154`), the
  corridor hold (`:163-164`), the abort and attempt-box guards
  (`guidance/arch.py:187-191`), the solver-status fallbacks.

The exclusion is appropriate for the C path — none of this fits the
interval-analyzable subset. It implies, however, that the decisions
with the highest safety leverage — when to abort, whether to keep
closing, which controller is in command — are exactly the mode/phase
machines inventoried above, and that they sit in the layer with the
least formal coverage: no ACSL, no EVA, no certificate, no exhaustive
check of any kind — only trace tests at fixed seeds. The certificate
layer re-checks the planners' *outputs* offline, never their
switching logic; the reachability lane verifies a hand-exported
abstraction of one FSM, with the duplicated-guard defect noted above.

Together, the two exclusions define the proposal's scope: the current
guarantees stop at a boundary the paper itself names; the artifacts
beyond that boundary are small, explicit state machines characterized
above; and model-checking them requires a `.tla` file and a CI lane,
with no change to the flight path, no new trusted arithmetic, and no
overlap with an existing modality.

## Candidate applications, ranked

TLA+ model-checks discrete transition systems exhaustively and is the
only tool under consideration that states and checks liveness under
fairness. The continuous side is unchanged: JuliaReach retains the
flowpipes and the barrier the infinite horizon; TLA+ covers the
discrete skeleton, with guards abstracted to nondeterministic "may
fire" events.

1. **ARCH FSM discrete abstraction** (highest value, smallest model).
   Modes, the clock, and abstracted guards give a TLC state space of
   trivial size. Safety invariants that today exist as comments
   become checked formulas: abort dominance (mode 3 wins any
   simultaneous enabling), mode monotonicity (order never decreases),
   mode-3 absorption. The principal gain is the liveness property:
   with fairness on the abort deadline,
   `<>[](mode ∈ {attempt, aborting})` closes the vacuity gap — the
   property the reachability lane structurally cannot pose, currently
   backstopped by five corner simulations.
2. **Mission phase machine with a truth/estimate split.** Model
   `phase`, the corridor gate, and *two* copies of the lateral state —
   truth and estimate — related by an assumed error bound. The safety
   invariant "never closing while misaligned (in truth)" is provable
   only under an explicit estimator-error assumption, which is the
   accurate form of the guarantee: the recorded `q_accel` bug — a
   near-open-loop filter accumulated ~0.1 m of lateral bias and "the
   ESTIMATE-keyed corridor gate never saw the true offset"
   (`sim/mission.py:119-125`) — is the counterexample TLC produces
   when the assumption is dropped. Additional low-cost invariants:
   replan-queue and burn-index monotonicity, phase one-wayness.
3. **Single-sourcing the guards** (a code fix TLA+ motivates but does
   not require). Derive `_in_attempt_box`, the exported halfspaces,
   and any future `.tla` guard from one datastructure, and draw mode
   names from `MODE_NAMES` everywhere. This removes the duplicated-
   guard defect regardless of whether a spec is ever written.
4. **Solver status protocol.** "Callers handle every status"
   (static-subset rule 5) is an unchecked convention over strings.
   The fix is an enum plus exhaustive-branch discipline (and
   eventually the prescribed C shape `status = step(...)`); TLA+ is
   marginal here — worth modeling only if statuses start driving mode
   transitions.
5. **The flight loop, where the two exclusions meet.** The cFS
   example app is a bare receive-dispatch loop
   (`examples/cfs_nav_app/podium_nav_app.c:70-85`); the sequential
   assumption (Exclusion 1) and the untrusted event-triggered layer
   (Exclusion 2) intersect here as soon as that loop grows real
   logic. Stale measurements, mode commands racing telemetry, and
   run-time-assurance switching (the documented route for any learned
   component) are standard model-checking subjects. Writing the
   specification before that logic exists minimizes its cost.

## Properties expressible in no current lane

The list above is organized by artifact and effort. The same content
organized by property — each inexpressible in every existing
modality:

1. **Mode entry.** `<>[](mode ∈ {attempt, aborting})` — the attempt
   mode is eventually reached and held, or abort takes over. The
   reachability gate only filters flowpipe segments by location, so
   an unentered mode 2 yields PROVEN (`arch_rendezvous.jl:71-76`);
   the trace side prints `+inf` "vacuously satisfied"; the STL `F`
   operators judge one trace. The only current check is
   `assert modes[-1] == 2` on five corner simulations. TLC checks the
   property over every run of the abstraction; fairness on the abort
   deadline is the assumption that makes it provable.
2. **Guard priority and simultaneous enabling.** Abort dominance,
   mode monotonicity, and mode-3 absorption exist today as an `elif`
   ordering plus the comment at `arch.py:187`. These are functional
   properties outside EVA's scope — a swapped guard priority is
   alarm-free wrong behavior, faithfully compiled by CompCert and
   bit-exactly reproduced by the golden vectors. The Python
   simulation additionally checks guards on a fixed dt grid
   ("switching times are quantized to dt," `arch.py:173-175`), so
   simultaneous enabling — abort deadline and attempt-box crossing in
   the same instant — is quantized away; TLC enumerates it
   exhaustively.
3. **Conditional safety across the truth/estimate boundary.** "Never
   closing while misaligned in truth, given |estimate − truth| ≤ b"
   is a two-variable invariant under an explicit assumption. Nothing
   in the current stack relates the estimate to the truth formally —
   the EKF returns innovations for monitors no caller runs, and
   simulation samples a few noise realizations. In TLA+ the
   estimator-error bound is a stated assumption; removing it makes
   TLC produce the recorded `q_accel` failure (the corridor gate
   blind to a real 0.1 m bias) as a counterexample trace. A bug found
   once by measurement becomes a property class checked on every
   relevant change.
4. **The sequencing obligation of Exclusion 1.** "OnMeas is atomic
   with respect to the filter state, and each app's messages are
   consumed one at a time in a total order" is the unstated contract
   every kernel proof rests on. ACSL as used has no vocabulary for
   interleavings; EVA and CompCert are sequential-semantics tools. A
   TLA+ model of tasks and the software bus proves the integration
   maintains the contract under all message orders — or exhibits the
   interleaving that breaks it. No existing lane expresses this
   property.
5. **General liveness under fairness.** "The corridor hold cannot
   livelock forever (given the estimator converges)," "every queued
   replan eventually fires," "a capture attempt eventually resolves
   as captured or bounced," "abort remains reachable" (checked as the
   invariant that the abort guard is never permanently disabled — the
   direct possibility form is branching-time, outside TLC's
   linear-time checking). All are currently observations about
   particular seeds;
   none is expressible as a certificate, a flowpipe, or an interval
   analysis, because all quantify over the futures of all runs under
   fairness assumptions. This is the half of Lamport's taxonomy for
   which the current stack has no prover.
6. **Protocol completeness for the status machine.** "Every caller
   has a transition for every solver status" — today a convention
   (rule 5) over strings. In a model, an unhandled status appears as
   a deadlock TLC reports. Low value at present (candidate 4 above);
   higher if statuses begin driving mode transitions.
7. **A drift-resistant specification artifact.** A capability rather
   than a property: today the FSM's semantics exist as a comment and
   its guard exists in two hand-synced copies (`_in_attempt_box` vs.
   the exported halfspaces — nothing checks agreement). A `.tla`
   module plus annotation extraction makes the discrete design a
   machine-checked, CI-gated artifact — the traceability rule 6
   demands, currently unenforced.

What TLA+ does not add: anything on the continuous side. It does not
check flowpipes, floating point, or exact rationals; the barrier
keeps the infinite horizon, JuliaReach the reachable sets, EVA the
RTE-freedom. TLA+ covers the one layer those tools abstract away —
the discrete logic that decides when to abort, whether to keep
closing, and which controller is in command — with the same
exhaustive, CI-replayed treatment as the rest of the stack. The cost
is a `.tla` file and a TLC lane; the flight path and the trusted
arithmetic are unchanged.

## Adoption path: source-annotation traceability

The failure mode of a standalone `.tla` file is drift — the spec and
the code diverge silently, and the mapping between spec variables and
code fields exists only in one developer's head. The
annotation-driven methodology (developed in the companion
`tla-annotation` work, applied previously to a Raft implementation)
addresses this: structured markers in the source declare the mapping
— variable mappings (`@tla{var: mode}`), action sites
(`TLA_ACTION("Abort")` on the guard block at `arch.py:187-191`),
named invariants and temporal properties — and an extraction tool
cross-references them against the `.tla` module in CI, reporting
unmapped variables and stale annotations. The markers translate
directly to Python as structured comments; the extraction layer is
language-agnostic. The convention is specified normatively in
[add-tla-specs.md](add-tla-specs.md).

Four reasons this fits Podium:

- **It matches the existing contract pipeline's shape.** `@contract`
  is already "declare once in Python, consume as runtime check + ACSL
  + harness"; `@tla{...}` applies the same pattern to discrete state,
  and it gives the currently unused `prove()` channel a consumer: an
  invariant declared at a program point becomes a named TLA+
  invariant checked by TLC, ahead of the planned ACSL lowering.
- **It operationalizes rule 6.** The rule demands traceability
  between the mode FSM and the verified hybrid model; annotations
  make that traceability extractable and CI-gated, and would have
  flagged the duplicated attempt-box guard as an unmapped-variable
  finding.
- **The spec effort is small and assisted.** Both candidate machines
  are already fully characterized (states, guards, one-wayness,
  absorption); generating the initial `.tla` from annotated source
  and having the developer review it — rather than writing specs from
  scratch — is the intended workflow.
- **The integration layer is generated, so its annotations can be
  too.** The cFS example app is emitted by `podium.emit.cfsapp`
  ("generated … DO NOT EDIT", `examples/cfs_nav_app/podium_nav_app.c:1`).
  The markers for the sequencing obligation of Exclusion 1 —
  `TLA_ATOMIC_BEGIN` around the OnMeas update→predict chain
  (`podium_nav_app.c:33-38`), action names on the receive-dispatch
  loop — can therefore be emitted by the generator rather than
  hand-maintained: drift-free by construction, the declare-once
  pattern of the contract pipeline extended to the app/bus layer
  where the sequential-call contract lives.

## Suggested sequence (not yet built)

The execution plan, with work packages, the normative annotation
convention, and acceptance criteria, is
[add-tla-specs.md](add-tla-specs.md). In outline:

1. Pilot: `tla/ArchRendezvous.tla` for candidate 1, TLC-checked in a
   CI lane alongside `reach.yml`; single-source the guards while
   touching `arch.py` (candidate 3).
2. `tla/Mission.tla` with the truth/estimate split (candidate 2).
3. Python annotation markers + extraction gate, applied to
   `guidance/arch.py` and `sim/mission.py`.
4. `tla/NavApp.tla`: model the cFS app, the software bus, and message
   orders against the *existing* example app, stating the
   sequential-call contract of Exclusion 1 (OnMeas atomicity,
   in-order consumption) as invariants plus an ordering assumption —
   actionable now, since the app and its `static` state already
   exist, and emitting its annotations from `podium.emit.cfsapp` at
   the same time.
5. Extend `NavApp.tla` as flight logic grows — mode commands, stale-
   measurement policy, run-time-assurance switching (candidate 5) —
   rather than writing a spec after the fact.
