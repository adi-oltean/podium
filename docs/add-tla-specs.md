# Plan: machine-checked TLA+ specifications for Podium's discrete mode logic

Status: implemented (2026-08-12); execution notes and deviations in
Section 10. Prerequisite reading: [verification.md](verification.md);
the gap analysis motivating this plan is in
[tla-potential.md](tla-potential.md). This document is self-contained:
every convention it introduces is specified normatively below, and no
external material is required to execute it.

## 1. Objective

Add a model-checking lane to Podium's verification stack covering the
discrete transition systems — mode logic, phase logic, and the flight
application's message-handling protocol — that the existing lanes do
not address. Deliverables are TLA+ specifications checked exhaustively
by TLC in continuous integration, plus a source-annotation convention
and extraction tool that keep the specifications traceable to the code.

The properties targeted are those inexpressible in the current
modalities:

| # | Property class | Example | Why no current lane can state it |
|---|---|---|---|
| P1 | Liveness under fairness | attempt mode is eventually entered, or abort takes over | reachability filters flowpipe segments by location and passes vacuously if a location is unreachable (`tools/reach/arch_rendezvous.jl:71-76`); STL `F` operators judge one trace |
| P2 | Transition-priority invariants | abort dominates simultaneous enabling; mode order is monotone; the abort mode is absorbing | encoded today as `elif` ordering plus a comment (`src/podium/guidance/arch.py:187`); Frama-C/EVA proves RTE-freedom, not functional correctness; simulation quantizes guard times to the dt grid (`arch.py:173-175`) |
| P3 | Conditional safety across the truth/estimate boundary | no closing while laterally misaligned *in truth*, given a stated estimator-error bound | no artifact relates estimate to truth; the failure class is on record (`src/podium/sim/mission.py:119-125`) |
| P4 | Sequencing obligations of the integration layer | filter update→predict is atomic; bus messages are consumed one at a time in order | ACSL as used has no vocabulary for interleavings; EVA and CompCert operate under sequential semantics |

Out of scope (covered by existing lanes, unchanged): continuous
dynamics and set reachability (JuliaReach), infinite-horizon abort
safety (exact-rational barrier), RTE-freedom of the emitted C
(Frama-C/EVA), compilation correctness (CompCert), numerical
equivalence (golden vectors). The TLA+ models abstract continuous
guards as nondeterministic events; no TLA+ artifact replaces or
duplicates a continuous-domain proof.

Assessed and deliberately deferred (from
[tla-potential.md](tla-potential.md), with reasons):

- *Solver status protocol* (candidate 4 there): the useful change —
  status enums plus exhaustive-branch handling in callers — is a code
  change independent of TLA+; modeling is deferred until statuses
  drive mode transitions.
- *Capture latch* (`sim/contact.py:171-187`): its verdict logic
  classifies a physics rollout post hoc and holds no control
  authority; revisit if it acquires any.
- *Possibility properties* ("from any state, abort remains
  reachable"): branching-time, not expressible in TLC's linear-time
  checking; if needed, reformulate as the invariant that the abort
  guard is never permanently disabled.

## 2. Constraints

C1. No new runtime dependencies. Annotations are comments; the TLC
    model checker (Java) is a CI-only dependency, pinned by version
    and checksum.

C2. The audit bundle's byte-determinism is not affected. `kernels.c`,
    `eva_driver.c`, and `bundle.json` (the `SHA256SUMS` manifest set,
    `tools/build_audit_bundle.py:75-78`) are not modified by any work
    package. WP4 modifies the generated example application only.

C3. The reachability lane's inputs are not altered. WP1 refactors
    guard definitions in `arch.py` to a single source; the exported
    hybrid-automaton JSON must remain value-identical (asserted by
    test, Section 7).

C4. Specifications state their environment assumptions explicitly
    (fairness conditions, estimator-error bounds, message-delivery
    assumptions). A property proved under an assumption is reported
    with that assumption; no assumption is left implicit.

## 3. Tooling

- TLA+ tools: `tla2tools.jar`, release pinned (v1.8.0 or later);
  record the release SHA-256 in the workflow on first download.
- Java 17 (Temurin) in CI; TLC invoked headless:
  `java -cp tla2tools.jar tlc2.TLC -workers auto -deadlock
  -config <model>.cfg <Module>.tla`.
  Note: the `-deadlock` flag *disables* TLC's deadlock reporting; it
  is required here because the models terminate by stuttering at the
  horizon (Section 5).
- No TLA+ Toolbox/IDE requirement; plain-text `.tla`/`.cfg` files.
- Extraction tool: `tools/tla_extract.py`, standard library only.

## 4. Annotation convention (normative)

Traceability between source and specification is maintained through
structured comments. The convention is language-neutral; Python uses
`# @tla{...}`, C uses `/* @tla{...} */`. An annotation applies to the
statement on the same line, or to the statement immediately following
when placed on its own line.

Grammar:

```
annotation := "@tla{" pair ("," pair)* "}"
pair       := key ": " value
key        := "module" | "spec" | "var" | "const" | "action"
            | "action-begin" | "action-end"
            | "atomic-begin" | "atomic-end"
            | "invariant" | "property" | "note"
```

Semantics:

| Key | Placement | Meaning |
|---|---|---|
| `module`, `spec` | once per source file | binds the file to a TLA+ module and its `.tla` path; enables extraction |
| `var` | on a variable/field declaration or first assignment | the code entity implements the named TLA+ `VARIABLE`; value `none` with a `note` marks state deliberately not modeled |
| `const` | on a constant/enum | maps to a TLA+ `CONSTANT` or model value set |
| `action` | on the statement(s) implementing one atomic spec transition | the code site implements the named action (a top-level operator in the module's `Next`) |
| `action-begin` / `action-end` | around a multi-statement block | block form of `action` |
| `atomic-begin` / `atomic-end` | around a region that the spec treats as one step | documents the atomicity boundary and what enforces it (`note` names the mechanism, e.g. "single task, single pipe") |
| `invariant` / `property` | at the code site the named formula constrains | reference by *name only*; the formula exists once, in the `.tla` file, never duplicated in a comment |
| `note` | any | free-text qualifier |

Design rule: the `.tla` file is the single source of formulas;
annotations carry names, not formulas, so a formula edit cannot
silently diverge from a comment.

Where a runtime check is also desired at an `invariant` site, the
code may call the existing `podium.verify.contracts.prove(cond,
label)` with `label` equal to the TLA+ invariant name; extraction
then cross-checks the label against the spec, giving that
currently-unused channel its first consumer.

Extraction tool (`tools/tla_extract.py`) contract:

- Inputs: source files carrying a `module` binding; the referenced
  `.tla` and `.cfg` files.
- Parses from the spec: `VARIABLES`, `CONSTANTS`, top-level operator
  names, and the `INVARIANT`/`PROPERTY` entries of the config.
- Reports, with file:line:
  - E1 *stale reference* (error): an annotation names an identifier
    absent from the bound spec/config.
  - E2 *unmapped variable* (error): a spec `VARIABLE` with no `var`
    annotation in any bound source file.
  - E3 *unanchored action* (warning initially; promoted to error once
    coverage is complete): an operator in `Next` with no `action`
    site.
  - E4 *unbalanced block* (error): `action-begin`/`atomic-begin`
    without its matching end.
- Exit code 0 iff no errors; runs in the CI lane (Section 6).

## 5. Work packages

### WP1 — `tla/ArchRendezvous.tla`: the ARCH rendezvous mode logic

Subject: the three-mode switched controller of the ARCH-COMP
rendezvous benchmark (`src/podium/guidance/arch.py:184-197`): modes
approaching → attempt → aborting, urgent guards, abort dominant, abort
absorbing. The continuous question (does the trajectory enter the
attempt octagon; does it respect the cone/velocity/keep-out sets)
remains with the reachability lane; the specification models the
discrete skeleton with box entry as a nondeterministic environment
event and time as a discrete clock in minutes.

Module skeleton (starting point; fairness placement is a review item):

```tla
---------------------- MODULE ArchRendezvous ----------------------
(* Discrete abstraction of src/podium/guidance/arch.py mode logic.
   Box entry is an environment event; continuous semantics are owned
   by the reachability lane (tools/reach). Time unit: minutes. *)
EXTENDS Naturals

CONSTANTS AbortTime, Horizon
ASSUME AbortTime \in Nat /\ Horizon \in Nat \ {0}
  (* "no abort" (SRNA) is modeled as AbortTime > Horizon *)

VARIABLES mode, clock, inBox
vars == <<mode, clock, inBox>>

Modes   == {"approaching", "attempt", "aborting"}
ModeOrd == [m \in Modes |->
             IF m = "approaching" THEN 1
             ELSE IF m = "attempt" THEN 2 ELSE 3]

TypeOK ==
  /\ mode \in Modes
  /\ clock \in 0..Horizon
  /\ inBox \in BOOLEAN

Init == mode = "approaching" /\ clock = 0 /\ inBox = FALSE

AbortEnabled == clock >= AbortTime

Abort ==                        \* urgent, dominant
  /\ AbortEnabled /\ mode # "aborting"
  /\ mode' = "aborting" /\ UNCHANGED <<clock, inBox>>

EnterAttempt ==                 \* urgent; dominated by Abort
  /\ mode = "approaching" /\ inBox /\ ~AbortEnabled
  /\ mode' = "attempt" /\ UNCHANGED <<clock, inBox>>

EnvEnterBox ==                  \* environment event (latched; see note)
  /\ mode = "approaching" /\ ~inBox
  /\ inBox' = TRUE /\ UNCHANGED <<mode, clock>>

Tick ==                         \* time advances only when nothing is urgent
  /\ clock < Horizon
  /\ ~(AbortEnabled /\ mode # "aborting")
  /\ ~(mode = "approaching" /\ inBox /\ ~AbortEnabled)
  /\ clock' = clock + 1 /\ UNCHANGED <<mode, inBox>>

Done == clock = Horizon /\ UNCHANGED vars

Next == Abort \/ EnterAttempt \/ EnvEnterBox \/ Tick \/ Done

Spec == Init /\ [][Next]_vars
             /\ WF_vars(Abort) /\ WF_vars(EnterAttempt) /\ WF_vars(Tick)

(* --- checked formulas ------------------------------------------- *)
AbortByDeadline   == clock > AbortTime => mode = "aborting"
AttemptRequiresBox == mode = "attempt" => inBox
Monotone  == [][ModeOrd[mode'] >= ModeOrd[mode]]_vars
Absorbing == [][mode = "aborting" => mode' = "aborting"]_vars
AbortTaken   == <>[](mode = "aborting")            \* SRA config only
EntryFromBox == (<> inBox) => <>(mode = "attempt")
===================================================================
```

Modeling notes for the executor:

- Urgency is encoded in `Tick`'s guard: time cannot advance while an
  urgent transition is enabled. This turns the bounded-liveness claim
  "abort engaged by the deadline" into the *state invariant*
  `AbortByDeadline` — the clock-augmentation idiom the code base
  already uses (`arch.py:13-14`).
- `EnvEnterBox` is latched (one-way) in this baseline, matching the
  hybrid model in which the attempt-mode invariant confines the state
  to the box. A variant with box exit (`inBox' = FALSE` permitted)
  should be checked as a robustness experiment; the liveness
  properties then require strong fairness, and the difference is
  itself informative.
- `EnvEnterBox` carries no fairness: the environment is never assumed
  helpful. `EntryFromBox` is conditional liveness — *if* the
  continuous layer delivers box entry (which reachability/simulation
  evidence supports), the discrete layer converts it to mode entry.

Two configurations:

```
\* tla/ArchRendezvousSRA.cfg          \* tla/ArchRendezvousSRNA.cfg
SPECIFICATION Spec                    SPECIFICATION Spec
CONSTANT AbortTime = 120              CONSTANT AbortTime = 301
CONSTANT Horizon = 300                CONSTANT Horizon = 300
INVARIANT TypeOK                      INVARIANT TypeOK
INVARIANT AbortByDeadline             INVARIANT AttemptRequiresBox
INVARIANT AttemptRequiresBox          PROPERTY Monotone
PROPERTY Monotone                     PROPERTY Absorbing
PROPERTY Absorbing                    PROPERTY EntryFromBox
PROPERTY AbortTaken
PROPERTY EntryFromBox
```

Code changes in `src/podium/guidance/arch.py` (single-sourcing the
guards, motivated independently of TLA+):

1. Define the attempt octagon once, as module-level data
   (`ATTEMPT_BOX = ((idx, coef, bound), ...)` — eight halfspaces).
   Derive `_in_attempt_box` by evaluating that data and
   `export_model`'s `attempt_box` list from the same data. Today the
   two are independent transcriptions (`arch.py:91-99` vs
   `:265-270`) with no consistency check.
2. Draw mode names in `export_model` from `MODE_NAMES` (`:76`) instead
   of re-typed string literals (`:279-291`).

Annotations to add (Python form):

| Site | Annotation |
|---|---|
| module docstring end | `# @tla{module: ArchRendezvous, spec: tla/ArchRendezvous.tla}` |
| `MODE_NAMES` (`:76`) | `# @tla{const: Modes}` |
| `mode = 1` in `simulate` | `# @tla{var: mode}` |
| `t = k * dt` in `simulate` | `# @tla{var: clock}` |
| `_in_attempt_box` def | `# @tla{var: inBox, note: guard proxy; continuous semantics owned by the reachability lane}` |
| abort branch (`:188-189`) | `# @tla{action: Abort}` and `# @tla{invariant: AbortByDeadline}` |
| attempt branch (`:190-191`) | `# @tla{action: EnterAttempt}` |
| `export_model` transition list | `# @tla{action: Abort, note: exported guard, same transition as simulate}` (and likewise `EnterAttempt`) |

### WP2 — traceability tooling and CI lane

1. Implement `tools/tla_extract.py` per the contract in Section 4.
2. Add `.github/workflows/tla.yml`:

```yaml
name: tla
on:
  push:
    branches: [main]
    paths:
      - "tla/**"
      - "src/podium/guidance/arch.py"
      - "src/podium/sim/mission.py"
      - "src/podium/emit/cfsapp.py"
      - "tools/tla_extract.py"
      - ".github/workflows/tla.yml"
  workflow_dispatch:
jobs:
  tlc:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-java@v4
        with: { distribution: temurin, java-version: "17" }
      - name: Fetch TLC (pinned release)
        run: |
          curl -fsSL -o tla2tools.jar \
            https://github.com/tlaplus/tlaplus/releases/download/<TAG>/tla2tools.jar
          echo "<RECORDED_SHA256>  tla2tools.jar" | sha256sum -c
      - name: Model-check
        run: |
          set -e
          java -cp tla2tools.jar tlc2.TLC -workers auto -deadlock \
            -config tla/ArchRendezvousSRA.cfg tla/ArchRendezvous.tla
          java -cp tla2tools.jar tlc2.TLC -workers auto -deadlock \
            -config tla/ArchRendezvousSRNA.cfg tla/ArchRendezvous.tla
          java -cp tla2tools.jar tlc2.TLC -workers auto -deadlock \
            -config tla/Mission.cfg tla/Mission.tla
          java -cp tla2tools.jar tlc2.TLC -workers auto -deadlock \
            -config tla/NavApp.cfg tla/NavApp.tla
      - name: Annotation extraction gate
        run: python3 tools/tla_extract.py --strict
```

   (Mission and NavApp lines land with WP3/WP4. Record the TLC
   state-space size from each log as a job summary — growth is a
   review signal.)

3. Extraction runs in the same job and fails the lane on E1/E2/E4.

### WP3 — `tla/Mission.tla`: phase logic with a truth/estimate split

Subject: the mission controller's discrete decisions
(`src/podium/sim/mission.py`): the one-way phase transition A → BC
(`:153-154`), the corridor hold gate (`:163-164`), the replan queue
(`:128,142-144`), and the burn index (`:146-152`). The distinguishing
feature is an explicit truth/estimate split: the gate reads the
*estimate*; the safety claim concerns the *truth*. The failure class
is on record — an under-weighted process noise left the filter
near-open-loop, a ~0.1 m lateral bias accumulated, and the
estimate-keyed gate never saw the true offset (`mission.py:119-125`).

State (discretized for TLC; suggested granularity in brackets):

| Variable | Domain | Implements |
|---|---|---|
| `phase` | `{"A", "BC"}` | `state["phase"]` |
| `t` | `0..HorizonTicks` | mission clock |
| `truthOff` | `0..MaxOff` [decimeters] | true lateral offset (environment) |
| `estOff` | `0..MaxOff` | estimated lateral offset |
| `near` | `BOOLEAN` | inside the 40 m gate region |
| `closing` | `BOOLEAN` | commanded closing rate nonzero |
| `replans` | `Seq` or counter | pending replan times |
| `burnIdx` | `Nat` | plan burn cursor |

Actions: `AdvancePhase` (time-triggered, one-way), `Replan` (consumes
the head of `replans`, resets `burnIdx`), `ExecuteBurn` (increments
`burnIdx`), `GateEvaluate` (sets `closing` from `estOff`/`near` per
the code's predicate), `EnvDrift` (environment updates `truthOff`),
`Estimate` (updates `estOff`; conjoined with the gate-local error-bound
assumption below), `Tick`.

Stated assumption (the estimator-error bound at BC near-gate decisions,
a `CONSTANT B`):

```tla
GateRelevant == phase = "BC" /\ near
EstimateOK ==
  (pc = "act" /\ GateRelevant) => AbsDiff(estOff, truthOff) <= B
```

`Estimate` maintains `EstimateOK` by bounding its choice only when
`GateRelevant`; estimates remain adversarial elsewhere, and the environment
never constrains `truthOff` beyond its domain.

Checked formulas:

```tla
CorridorHoldSound ==            \* the central conditional-safety claim
  (phase = "BC" /\ near /\ truthOff > Tol) => ~closing
PhaseMono   == [][phase = "BC" => phase' = "BC"]_vars
BurnMono    == [][burnIdx' >= burnIdx]_vars
ReplanMono  == [][Len(replans') <= Len(replans)]_vars
HoldReleases ==                  \* liveness, under estimator fairness
  [](( ~closing /\ phase = "BC") => <>(closing \/ t = HorizonTicks))
ReplanLive ==                    \* queued replans fire during phase A
  <>(replans = <<>> \/ phase = "BC")
```

`CorridorHoldSound` is expected to hold only for `B <= Tol`-class
bounds; the required relation between `B` and `Tol` is a result of
the work package, stated in the spec's comments. The negative control
(Section 7) removes the bound and must reproduce the recorded failure
shape as a TLC counterexample.

Annotations: `module`/`spec` binding at the top of `mission.fly`;
`var` on `state` dict keys (`phase`, `burn_idx`), the `replan_times`
list, and the gate variables; `action` on the four decision sites;
`invariant: CorridorHoldSound` on the gate lines (`:163-164`).

### WP4 — `tla/NavApp.tla`: the flight application's message protocol

Subject: the generated cFS application
(`examples/cfs_nav_app/podium_nav_app.c`, emitted by
`podium.emit.cfsapp`). Its measurement handler performs a
two-kernel read-modify-write chain on module state
(`update_sequential` then `predict` on `static PODIUM_NAV_App_t g`,
`:33-38`), currently made atomic by structure alone: one task, one
pipe, blocking receive (`:75-83`). The kernel-level proofs (EVA,
CompCert, golden vectors) are established under sequential semantics;
this specification states, for the first time, the integration-level
contract those proofs assume.

Model: a bounded message source, a FIFO bus queue, and the handler as
a two-step program counter (the state between update and predict is
*inconsistent* and must not be observable).

```tla
CONSTANTS NMsgs, QCap
VARIABLES sent, queue, pc, processed
  \* pc \in {"idle", "updated"}: "updated" is the inconsistent
  \* interior of the OnMeas chain (between update and predict)

Recv    == queue # <<>> /\ pc = "idle" /\ ...      \* dequeue head
Update  == pc = "idle"    /\ pc' = "updated" /\ ...
Predict == pc = "updated" /\ pc' = "idle"
           /\ processed' = Append(processed, current)

InOrder     == processed = SubSeqOrdered(sent)   \* prefix-order check
NoPartialPublish == (* telemetry reads occur only when pc = "idle" *)
AllProcessed == <>(Len(processed) = NMsgs)       \* under WF(Recv,...)
```

Baseline (single task, no reader): all three formulas must verify.
Experiment E-NAV (required, documented in the spec file): add a
`Reader` action that may fire when `pc = "updated"`. TLC must produce
the violating interleaving of `NoPartialPublish` — the concrete
witness that the atomicity obligation is real, not hypothetical. This
mirrors the falsification-receipt practice used elsewhere in the
repository (`tests/test_arch.py:1-8`).

Code changes: extend the `podium.emit.cfsapp` generator to emit the
annotations into the generated C (they are comments; no object-code
change):

| Site (generated) | Annotation |
|---|---|
| file header | `/* @tla{module: NavApp, spec: tla/NavApp.tla} */` |
| `static PODIUM_NAV_App_t g;` | `/* @tla{var: pc, note: filter-state consistency; also processed} */` |
| before the update call | `/* @tla{atomic-begin: OnMeas, note: single task, single pipe, blocking receive} */` |
| after the second copy-back | `/* @tla{atomic-end: OnMeas} */` |
| dispatch branch | `/* @tla{action: Recv} */` |

Because the app is generated ("DO NOT EDIT", `podium_nav_app.c:1`),
the annotations are maintained by the generator, not by hand;
regeneration cannot drift from the emitter's declaration. Constraint
C2 is unaffected (the app is not part of the audit manifest).

## 6. CI integration summary

New lane `tla` (WP2): TLC over every `.cfg` + extraction gate; fails
on any invariant/property violation or extraction error. Path-filtered
to the spec files, the annotated sources, and the emitter. Runtime is
expected in seconds (state spaces are small by construction; report
sizes in the job summary). Existing lanes (`ci`, `reach`, `eva`,
`compcert`, `tier2`, `validate`, `release`, `pages`) are unchanged.

## 7. Acceptance criteria and negative controls

Positive criteria:

| # | Criterion |
|---|---|
| A1 | TLC verifies all WP1 formulas in both configurations; state count reported |
| A2 | `arch.py` guard single-sourcing: exported JSON value-identical to pre-change (new regression test compares `export_model()` output against the previous halfspace values exactly); `tests/test_arch.py` and the `reach` lane pass unchanged |
| A3 | TLC verifies WP3 formulas with the estimator bound stated; the bound/tolerance relation is documented in the module |
| A4 | TLC verifies WP4 baseline; experiment E-NAV yields the expected counterexample |
| A5 | `tla_extract.py --strict` exits 0 on the annotated tree; E2 coverage: every spec `VARIABLE` mapped |
| A6 | Constraint C2 verified: audit bundle `SHA256SUMS` entries unchanged by the plan's merges |

Negative controls (each must *fail* when the listed mutation is
applied, then be reverted; these are spec-level mutation receipts,
following the practice of `tools/fault_coverage.py`):

| # | Mutation | Expected failure |
|---|---|---|
| N1 | swap guard priority in `ArchRendezvous.tla` (`EnterAttempt` no longer excludes `AbortEnabled`) | `AbortByDeadline` or `Monotone` violated |
| N2 | remove `Abort`'s urgency from `Tick` | `AbortByDeadline` violated |
| N3 | drop the estimator bound in `Mission.tla` (set `B = MaxOff`) | `CorridorHoldSound` violated; counterexample exhibits estimate/truth divergence at the gate |
| N4 | E-NAV reader enabled | `NoPartialPublish` violated with an explicit interleaving |
| N5 | rename a `VARIABLE` in a spec without touching annotations | extraction E1/E2 failure |

## 8. Milestones

| M | Contents | Depends on |
|---|---|---|
| M1 | WP1 spec + configs; guard single-sourcing in `arch.py`; minimal `tla.yml` (TLC only) | — |
| M2 | WP2 extraction tool; annotations in `arch.py`; lane gates on extraction | M1 |
| M3 | WP3 Mission spec + annotations | M2 |
| M4 | WP4 NavApp spec; `cfsapp` emitter annotations; E-NAV experiment | M2 |

Each milestone is independently mergeable, and each is deliberately
small: the specifications and their state spaces are minimal by
construction. M3 and M4 are independent and may proceed in either
order. Future
flight logic (mode commands, stale-measurement policy, run-time-
assurance switching) extends `NavApp.tla` rather than starting a new
spec; that extension is out of this plan's scope but is the intended
continuation.

## 9. References

- L. Lamport, *Proving the Correctness of Multiprocess Programs*,
  IEEE Trans. Software Eng. SE-3(2), 1977. (Safety/liveness taxonomy.)
- L. Lamport, *Specifying Systems*, Addison-Wesley, 2002. (TLA+,
  TLC, fairness.)
- Repository: [verification.md](verification.md) (existing lanes),
  [tla-potential.md](tla-potential.md) (gap analysis behind this
  plan).

## 10. Execution notes (2026-08-12)

All four work packages are implemented; acceptance criteria A1–A6 and
negative controls N1–N5 all fired as required. Verified state spaces:
ArchRendezvous 725 (SRA) / 903 (SRNA) / 903 (Exit), Mission 16,336,
NavApp 35 (concrete depth 16) / 30 (depth-2 backpressure) distinct
states; E-NAV reaches its required violation after exploring 45 states.
TLC pinned at v1.8.0, SHA-256
`ab323b79802aedc3203b3f9af37c6aca3ed43f4e0225b36f2aa77b26de46c05f`
(recorded in `tla.yml`). Deviations from the plan as drafted, each
with a TLC receipt:

1. **`EntryFromBox` moved out of the SRA configuration.** As drafted
   (Section 5 config listing), TLC falsifies it there: `EnvEnterBox`
   may fire exactly at `clock = AbortTime`, `EnterAttempt` is already
   disabled, and `Abort` takes over — box entry racing the abort
   deadline in the same instant never yields attempt mode. This is
   precisely the simultaneous-enabling class P2 targets (the fixed-dt
   simulation quantizes it away). SRA checks the corrected conditional
   property `EntryBeforeAbort == <>(inBox /\ clock < AbortTime) =>
   <>(mode = "attempt")`; SRNA keeps `EntryFromBox`, where it is
   sound. The counterexample is documented in the module comment.
2. **`AbortDominates` added to both configurations.** N1's mutation
   (drop `~AbortEnabled` from `EnterAttempt`) leaves `AbortByDeadline`
   and `Monotone` intact — urgency stops the clock until `Abort`
   fires, so a swapped guard priority only inserts a transient attempt
   step; TLC confirms the masking. Dominance is expressible only at
   the action level: `[][~(AbortEnabled /\ mode' = "attempt")]_vars`.
   N1 falsifies exactly this formula.
3. **`BurnMono` is conditional.** The drafted unconditional form
   `[][burnIdx' >= burnIdx]_vars` is false — `Replan` resets the
   cursor, as the plan's own action list states. The checked form
   permits the reset exactly on replan steps; `BurnMonoUnconditional`
   is kept in `Mission.tla` (unchecked) for the one-command receipt.
4. **`CorridorHoldSound` is evaluated at the command instant** (the
   post-gate point of the tick cycle). Result: provable iff
   `Tol >= GateTol + B`; the shipped config sits on the boundary
   (`1 + 1 = 2`), and `Tol = 1` is falsified. Truth drift after the
   command, within a tick, is continuous-dynamics territory owned by
   the reachability/simulation lanes. N3 (`B = MaxOff`) reproduces the
   recorded `q_accel` failure shape: `estOff = 0`, `truthOff > Tol`,
   `near`, and the gate commands closing. The estimator assumption is scoped
   to this decision point: the gate cannot use or violate it in phase A or
   outside the near region.
5. **NavApp modeling details.** `pc` has a third state `"have"`
   (post-dispatch, pre-update: consistent but one cycle stale);
   `current` and `dirtyRead` are spec-side variables annotated at
   their most honest code sites. E-NAV is a permanent CI receipt, not
   a one-off: `ReaderEnabled` is a `CONSTANT`, `NavAppReader.cfg`
   must fail, and the lane inverts its exit code and greps for the
   `NoPartialPublish` violation. `NavApp.cfg` and `NavAppReader.cfg`
   use the generated app's concrete `QCap = 16`; the separate positive
   `NavAppBackpressure.cfg` retains `QCap = 2` to exercise a full pipe.
6. **Annotation coverage exceeds the drafted tables.** A5 requires
   every spec `VARIABLE` mapped (E2), so environment-side variables
   (`truthOff`, `sent`, `queue`, `processed`, `t`, `near`, `closing`,
   `replans`, `pc`, …) carry `var` annotations at the code sites that
   most honestly correspond to them, beyond the WP1/WP4 tables.
7. **The `prove()` consumer is deferred.** The `CorridorHoldSound`
   site cannot call `prove()`: the controller cannot observe the
   truth state — that is the point of the property. Runtime checking
   of invariants remains available for conditions computable in situ.
8. **Extractor `--strict` semantics.** Errors E1/E2/E3/E4 always fail
   the tool; `--strict` (the CI entry) additionally requires every
   `tla/*.tla` to be bound by at least one annotated source file. E3
   uses a closed, module-specific allowlist for environment/scheduler
   operators (`Tick`, `Done`, `EnvEnterBox`, `EnvExitBox`, `EnvDrift`,
   `EndActA`, `Send`, `Reader`) that deliberately have no code site;
   every other unanchored action, plus stale or redundant allowlist
   entries, is an error.
9. **A6 verified structurally.** The audit bundle's three inputs are
   untouched (`cemit`, `evagen`) or changed comments/docstrings only
   (`mission.py`, AST-identical modulo docstrings), so the
   `SHA256SUMS` entries are unchanged by construction.
11. **Trace validation added (beyond the plan).**
    `tools/tla_trace_check.py` converts concrete simulation runs to
    their discrete abstractions and has TLC check them as literal spec
    behaviors (generated witness modules constrain `Next` to the trace;
    deadlock detection ON, so divergence deadlocks at the first
    unmatched step). Subjects: `arch` — `simulate()` for SRA/SRNA from
    all five corners vs `ArchRendezvous` (~304 witness states each);
    `mission` — `fly(record_tla_events=True)` vs `Mission` from the
    opt-in discrete events in `extras["tla_events"]` (~19,000 witness states;
    offsets quantized in 0.06 m gate units so `GateTol = 1` is exactly
    the code's predicate; `GateTol`, `B`, and `Tol` are read from
    `Mission.cfg`, while only `MaxOff` is sized from the trace. The
    global and BC near-gate estimator errors are reported, and a gate-region
    out-of-bound trace is rejected rather than weakening the configured
    assumption). Tamper
    receipts (`--receipt`) require TLC to reject a corrupted witness.
    Two consequences for earlier items: (a) `EnvDrift`'s ±1 pace bound
    was removed — no checked property used it, the environment gets
    strictly more adversarial, and literal traces replay (Mission's
    state space first grows 9,144 → 10,576 and later to 16,336 with the
    gate-local estimator assumption, all formulas and receipts N3/R2/R3
    unaffected); (b) supersedes item 9's AST receipt for
    `mission.py`: trace validation opts into event records at the
    decision sites. Default `fly()` calls preserve the prior `extras`
    contract and avoid constructing trace events; no numeric path or
    RNG use is touched, so `bundle.json` remains byte-identical. CI:
    the `trace` job of `tla.yml`.
12. **Box-exit robustness experiment (WP1 modeling note) executed.**
    `ArchRendezvousExit.cfg` (constant `BoxLatched = FALSE`, spec
    `SpecSF`) verifies that all safety/priority properties survive
    revocable box entry and that sustained presence still converts to
    mode entry: `EntrySustainedBox == ([]<>inBox) => <>(mode =
    "attempt")` under strong fairness on `EnterAttempt` (903 distinct
    states). The informative differences, both TLC-falsified as
    receipts: under weak fairness the enter/exit flicker keeps
    `EnterAttempt` enabled infinitely often but never continuously, so
    `EntrySustainedBox` fails; and the unconditional `EntryFromBox`
    fails under ANY fairness once exits are allowed — a single revoked
    entry, never repeated, obliges no fairness assumption. The latched
    baseline (SRA/SRNA) is unchanged: `BoxLatched = TRUE` disables
    `EnvExitBox`, and both configs reverify at their prior state
    counts.
13. **The estimator bound is gate-local, not global.** The concrete seed-7
    trace measures `B = 5` globally but `B = 1` across all 1,878 BC
    near-gate evaluations. `Estimate` therefore constrains its choice only
    when `phase = "BC" /\ near`, exactly where `GateEvaluate` consumes the
    estimate; `EstimateOK` checks that pre-command state. The shipped
    `B = 1`, `Tol = 2` safety boundary is unchanged. TLC verifies all Mission
    properties over 16,336 distinct states, and the 19,220-state seed-7
    witness plus its tamper receipt pass without changing mission code.
