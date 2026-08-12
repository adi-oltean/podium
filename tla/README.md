# tla/ — model-checked specifications of Podium's discrete logic

This directory is the model-checking lane of the verification stack:
TLA+ specifications of the discrete transition systems — mode logic,
phase logic, and the flight application's message protocol — checked
exhaustively by TLC in CI (`.github/workflows/tla.yml`). It covers the
one layer the continuous lanes abstract away: the logic that decides
when to abort, whether to keep closing, and which controller is in
command. The continuous side is unchanged — JuliaReach keeps the
flowpipes, the exact-rational barrier the infinite horizon, EVA the
RTE-freedom; here, continuous guards are nondeterministic environment
events.

The gap analysis behind this lane is
[docs/tla-potential.md](../docs/tla-potential.md); the execution plan,
annotation convention (normative), and acceptance receipts are in
[docs/add-tla-specs.md](../docs/add-tla-specs.md).

## Specifications

| Module | Configs | Subject | Distinct states |
|---|---|---|---|
| `ArchRendezvous.tla` | `ArchRendezvousSRA.cfg` (abort at t=120), `ArchRendezvousSRNA.cfg` (no abort), `ArchRendezvousExit.cfg` (box exit allowed, strong fairness) | The ARCH rendezvous mode FSM of `src/podium/guidance/arch.py`: approaching → attempt → aborting, urgent guards, abort dominant and absorbing | 725 / 903 / 903 |
| `Mission.tla` | `Mission.cfg` | The mission controller's decisions in `src/podium/sim/mission.py` fly(): one-way phase A → BC, the corridor hold gate, replan queue, burn cursor — with an explicit truth/estimate split | 16,336 |
| `NavApp.tla` | `NavApp.cfg` (concrete depth 16), `NavAppBackpressure.cfg` (depth-2 stress), `NavAppReader.cfg` (experiment E-NAV, **must fail**) | The generated cFS app's message protocol (`examples/cfs_nav_app/podium_nav_app.c`): bus FIFO, two-step update→predict handler whose interior state is inconsistent | 35 / 30 (reader: 45 to violation) |

## What is proved (and under which stated assumptions)

**ArchRendezvous** — invariants `AbortByDeadline` (the clock cannot
pass the abort time without the abort engaging; urgency encoded in
`Tick`'s guard) and `AttemptRequiresBox`; action properties `Monotone`,
`Absorbing`, and `AbortDominates` (abort wins simultaneous enabling —
expressible only at the action level, see the module comment); liveness
`AbortTaken` (`<>[](mode = "aborting")`, SRA) and the conditional
mode-entry properties `EntryFromBox` (SRNA) / `EntryBeforeAbort`
(both). Fairness is on the controller and scheduler only; the
environment (box entry) is never assumed helpful. This closes the
vacuity gap of the reachability gate: mode entry is now checked over
every run of the abstraction, not five corner simulations. The
box-exit robustness experiment (`ArchRendezvousExit.cfg`, entry
revocable) grades the fairness needed for entry: weak fairness is
defeated by enter/exit flicker, strong fairness on `EnterAttempt`
recovers entry from *sustained* presence (`EntrySustainedBox`), and
one revoked entry defeats the unconditional form under any fairness —
each direction TLC-checked.

**Mission** — the central conditional-safety claim
`CorridorHoldSound`: never commanding closing while laterally
misaligned *in truth*, given the stated estimator-error bound
`|estOff − truthOff| <= B` immediately before BC near-gate decisions
(invariant `EstimateOK`). Estimates remain adversarial where the gate
cannot act. Result: provable iff `Tol >= GateTol + B`; the shipped config
sits exactly on that boundary. Plus `PhaseMono`, `ReplanMono`, `BurnMono` (monotone except
at replan resets), `HoldReleases`, `ReplanLive`. Dropping the bound
(`B = MaxOff`, negative control N3) makes TLC reproduce the recorded
`q_accel` failure shape: estimate reads aligned, truth is past
tolerance, the gate commands closing.

**NavApp** — the integration-level contract every kernel proof (EVA,
CompCert, golden vectors) silently assumes: `NoPartialPublish` (the
update→predict chain is atomic — nothing observes `g` mid-chain),
`InOrder` (messages consumed one at a time, in order), `AllProcessed`
(under fairness on the source and run loop). In the baseline these
hold by structure alone: one task, one pipe, blocking receive. The
baseline uses the generated app's concrete pipe depth of 16;
`NavAppBackpressure.cfg` separately retains depth 2 to exercise the
temporarily-full queue behavior.

**E-NAV** (`NavAppReader.cfg`) is a permanent falsification receipt,
not a one-off: it enables a second task that reads `g` without the
pipe, and CI *requires* TLC to produce the violating interleaving
(Send → Recv → Update → Reader). That counterexample is the concrete
witness that the atomicity obligation is real — it is what embedding
the kernels in a multi-task flight build would do without a locking or
message-ownership discipline.

## Traceability: the `@tla` annotations

The specs are not standalone files. Structured comments in the source
(`# @tla{...}` in Python, `/* @tla{...} */` in C) bind each file to its
module and map code entities to spec variables, actions, and checked
formulas — by *name only*; the `.tla` file is the single source of
formulas. `tools/tla_extract.py --strict` cross-checks both sides in
CI: stale references (E1), unmapped variables (E2), unanchored actions
outside the explicit environment/scheduler allowlist (E3), and
unbalanced blocks (E4) fail the lane; renaming a spec variable without
touching the code (or vice versa) is a build failure, not silent drift. The
annotations in the generated cFS app are emitted by
`podium.emit.cfsapp`, so regeneration cannot drift from the
declaration.

## Trace validation

`tools/tla_trace_check.py` closes the loop from the other side: instead
of checking all behaviors of the abstraction, it checks that the
*concrete* simulation runs are among them. A real run is converted to
the exact state sequence the spec would traverse — its discrete
abstraction — and a generated witness module constrains `Next` to
follow it, with TLC's deadlock detection ON: a divergence deadlocks at
the first unmatched step (TLC prints the matched prefix), conformance
completes cleanly. This upgrades the code–spec link from name-level
annotations to semantic conformance.

- `arch`: `simulate()` for both scenarios from all five initial corners
  (10 traces, ~304 witness states each) vs `ArchRendezvous.tla`; a
  per-minute sample that appears in no witness state at its clock is
  itself reported as divergence.
- `mission`: `fly(record_tla_events=True)` vs `Mission.tla`, from the
  opt-in discrete events in `extras["tla_events"]` plus the truth trace (~19,000
  witness states for the 4,800-tick mission). Offsets are quantized in
  gate-threshold units (1 unit = 0.06 m) so the spec's `GateTol = 1` is
  exactly the code's gate predicate. `GateTol`, estimator bound `B`, and
  truth tolerance `Tol` are read from `Mission.cfg`; only the finite
  domain size `MaxOff` is fitted to the trace. Global and BC near-gate
  estimator errors are reported; `EstimateOK` enforces the shipped bound
  where the gate can act, and `CorridorHoldSound` checks the resulting
  command. `--selftest` replays a synthetic
  full-length record (no heavy dependencies).
- `--receipt` tampers each witness (a two-tick clock jump) and requires
  TLC to REJECT it — the receipt that the harness can fail.

The `trace` job of `tla.yml` runs both subjects with receipts on every
relevant change; witness modules are generated in a temp directory, and
nothing under `tla/` is written.

## Running locally

One command runs the whole lane — every configuration, the E-NAV
falsification (required to fail), and the extraction gate — fetching
the pinned, checksum-verified TLC jar into the repo root (gitignored)
on first use:

```sh
tools/run_tla.sh
```

The individual steps, if you want them by hand (Java 11+ only):

```sh
curl -fsSL -o tla2tools.jar \
  https://github.com/tlaplus/tlaplus/releases/download/v1.8.0/tla2tools.jar
echo "ab323b79802aedc3203b3f9af37c6aca3ed43f4e0225b36f2aa77b26de46c05f  tla2tools.jar" | shasum -a 256 -c

java -cp tla2tools.jar tlc2.TLC -workers auto -deadlock \
  -config tla/ArchRendezvousSRA.cfg tla/ArchRendezvous.tla
# likewise ArchRendezvousSRNA, Mission, NavApp, NavAppBackpressure;
# NavAppReader must FAIL with "Invariant NoPartialPublish is violated"

python3 tools/tla_extract.py --strict
```

The `-deadlock` flag *disables* TLC's deadlock reporting — required
because the models terminate by stuttering at their horizon. Runs
finish in seconds; state-space growth is a review signal (CI posts the
counts in the job summary). On a property violation TLC writes a
`*_TTrace_*.tla` trace file next to the spec (gitignored).

Trace validation runs separately. The validator copies the specifications
into a temporary directory, runs TLC there, and does not modify the
repository's Eclipse `.project` file:

```sh
# Fast synthetic Mission trace (no optional Python dependencies)
python3 tools/tla_trace_check.py mission --selftest --receipt

# Ten concrete ARCH traces (requires Podium's NumPy dependency)
.venv/bin/python tools/tla_trace_check.py arch --receipt

# Full reference Mission, seed 7 (requires the opt + contact extras)
python3 -m pip install -e ".[opt,contact]"
python3 tools/tla_trace_check.py mission --seed 7 --receipt
```

Use `--jar /path/to/tla2tools.jar` if the pinned TLC JAR is not in the
repository root. The synthetic Mission and ARCH commands must pass and
their tamper receipts must be rejected. The real Mission command also
enforces the BC near-gate estimator bound from `Mission.cfg`; a gate-region
error above that assumption is a validation failure, not an automatically
weakened witness. The global error is diagnostic because the gate cannot act
outside that region.

## Receipts

Falsification receipts — mutations that must fail, in the spirit of
`tools/fault_coverage.py` — are kept as comments in the modules and
documented with their TLC counterexamples in
[docs/add-tla-specs.md](../docs/add-tla-specs.md) Sections 7 and 10:
swapped guard priority (falsifies `AbortDominates`, and *only* that —
urgency masks it from the deadline invariant), removed urgency
(falsifies `AbortByDeadline`), dropped estimator bound (falsifies
`CorridorHoldSound`), the E-NAV reader (falsifies `NoPartialPublish`),
and a spec-variable rename (fails extraction).

## What's deliberately left, and what would trigger it

- **Solver status protocol** ([tla-potential.md](../docs/tla-potential.md)
  candidate 4): deferred until solver statuses drive mode transitions;
  the useful near-term change is a code change (status enums +
  exhaustive branch handling in callers), not a spec.
- **Capture latch** (`src/podium/sim/contact.py`): excluded while it
  has no control authority — it classifies a physics rollout post hoc.
  Revisit if it ever gates a retry or abort.
- **Possibility properties** ("abort always remains reachable"):
  branching-time, outside TLC's linear-time checking; the linear-time
  reformulation (the abort guard is never permanently disabled) is
  available if needed.
- **The main open runway — extending `NavApp.tla` as flight logic
  grows**: mode commands racing telemetry, stale-measurement policy,
  run-time-assurance switching, multiple apps sharing the bus. This is
  where the plan says the biggest value accrues, and it is *supposed*
  to wait — the whole point is to write the spec additions as (or
  before) that logic appears, not after.

Beyond what the documents call for, the honest candidates in rough
value order (trace validation, formerly first on this list, is built —
see the Trace validation section above):

1. **E3 promotion** — an explicit environment-operator allowlist in
   `tools/tla_extract.py`, after which unanchored actions become
   errors, completing the plan's "promoted once coverage is complete."
2. **A `prove()` consumer** — `podium.verify.contracts.prove` still
   has zero call sites; it needs a future invariant whose condition is
   computable in situ (`CorridorHoldSound` isn't, since the controller
   cannot observe truth — that is the point of the property).
3. **Apalache / refinement** — only if state spaces outgrow TLC or a
   finer Mission model appears; unnecessary at 10³–10⁴ states.
