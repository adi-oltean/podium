------------------------------- MODULE Mission -------------------------------
(* Discrete abstraction of the mission controller's decisions in
   src/podium/sim/mission.py fly(): the one-way phase transition A -> BC,
   the corridor hold gate, the replan queue, and the burn cursor.  The
   distinguishing feature is an explicit truth/estimate split: the gate
   reads the ESTIMATE (estOff); the safety claim concerns the TRUTH
   (truthOff).  The failure class is on record — an under-weighted process
   noise left the EKF near-open-loop, a ~0.1 m lateral bias accumulated,
   and the estimate-keyed gate never saw the true offset (mission.py's
   q_accel comment).  Dropping the gate-region estimator bound
   (B = MaxOff, negative control N3) makes TLC reproduce exactly that shape.

   Model structure: one GNC tick is the pc cycle env -> est -> act -> tick,
   mirroring the strictly sequential body of controller() — the
   interleaving freedom is in the VALUES (truthOff, near, estOff), never
   in the order, which is what the code actually guarantees.  Offsets are
   discretized (about one unit per decimeter); truth moves without any
   pace bound — no checked property needs one, the environment is never
   assumed helpful, and unconstrained motion lets concrete fly() traces
   replay literally through this spec (tools/tla_trace_check.py).  The
   range axis is not truth/estimate-split: the
   recorded failure class is lateral, and a split range would add a second
   bound with the same structure.

   Stated assumptions (constraint C4 of docs/add-tla-specs.md):
     - estimator error bound: immediately before a BC near-gate decision,
       after Estimate, |estOff - truthOff| <= B (the EstimateOK invariant).
       The estimate is unconstrained where the gate cannot act; truthOff and
       near remain adversarial everywhere.
     - scheduler fairness: the GNC loop keeps running (WF on the cycle
       actions).  Fairness is on the loop EXECUTING, never on the
       environment being helpful.

   RESULT (checked, not assumed): CorridorHoldSound — never closing while
   misaligned in truth, at the instant the command is issued — holds iff
     Tol >= GateTol + B.
   The shipped config sits exactly on the boundary (GateTol = 1, B = 1,
   Tol = 2); B = MaxOff (N3) or Tol < GateTol + B is falsified by TLC.
   What happens to the truth WITHIN a tick, after the command is issued,
   is continuous drift owned by the reachability/simulation lanes; the
   discrete claim is about the decision.

   Falsification receipt: the unconditional burn-cursor monotonicity
   [][burnIdx' >= burnIdx]_vars is FALSE — Replan resets the cursor for
   the fresh plan (mission.py replan()).  BurnMono therefore permits the
   reset exactly on replan steps; TLC exhibits the reset trace if the
   unconditional form is checked instead. *)
EXTENDS Naturals

CONSTANTS
  HorizonTicks,   \* mission horizon [ticks]
  PhaseATicks,    \* the A -> BC time trigger (t_phase_a)
  NReplans,       \* length of the replan queue (replan_times)
  MaxBurns,       \* burn-cursor bound within one plan
  MaxOff,         \* lateral-offset domain 0..MaxOff [~decimeters]
  GateTol,        \* the gate's threshold on the ESTIMATE (lat > 0.06 m)
  B,              \* estimator-error bound at BC near-gate decisions
  Tol             \* truth-misalignment bound of the safety claim

ASSUME
  /\ HorizonTicks \in Nat \ {0}
  /\ PhaseATicks \in 1..(HorizonTicks - 1)
  /\ NReplans \in Nat \ {0}
  /\ MaxBurns \in Nat
  /\ MaxOff \in Nat \ {0}
  /\ GateTol \in 0..MaxOff /\ B \in 0..MaxOff /\ Tol \in 0..MaxOff

VARIABLES
  phase,     \* "A" | "BC"                      (state["phase"])
  t,         \* mission clock [ticks]
  pc,        \* tick sub-schedule: "env" -> "est" -> "act" -> "tick"
  truthOff,  \* true lateral offset             (environment)
  estOff,    \* estimated lateral offset        (EKF)
  near,      \* inside the 40 m gate region
  closing,   \* commanded closing rate nonzero  (v_close > 0)
  replans,   \* pending replans                 (len(replan_times))
  burnIdx    \* burn cursor into the current plan

vars == <<phase, t, pc, truthOff, estOff, near, closing, replans, burnIdx>>

AbsDiff(a, b) == IF a >= b THEN a - b ELSE b - a
GateRelevant == phase = "BC" /\ near

TypeOK ==
  /\ phase \in {"A", "BC"}
  /\ t \in 0..HorizonTicks
  /\ pc \in {"env", "est", "act", "tick"}
  /\ truthOff \in 0..MaxOff
  /\ estOff \in 0..MaxOff
  /\ near \in BOOLEAN
  /\ closing \in BOOLEAN
  /\ replans \in 0..NReplans
  /\ burnIdx \in 0..MaxBurns

Init ==
  /\ phase = "A" /\ t = 0 /\ pc = "env"
  /\ truthOff \in 0..MaxOff /\ estOff \in 0..MaxOff
  /\ near \in BOOLEAN /\ closing = FALSE
  /\ replans = NReplans /\ burnIdx = 0

\* replan_times[i] = i * t_phase_a / NReplans, consumed front to back
NextReplanTime == ((NReplans - replans) * PhaseATicks) \div NReplans
ReplanDue == phase = "A" /\ replans > 0 /\ t >= NextReplanTime

EnvDrift ==                     \* adversarial: truth and range unconstrained
  /\ pc = "env"
  /\ truthOff' \in 0..MaxOff
  /\ near' \in BOOLEAN
  /\ pc' = "est"
  /\ UNCHANGED <<phase, t, estOff, closing, replans, burnIdx>>

Estimate ==                     \* nav.step: gate-local error assumption
  /\ pc = "est"
  /\ estOff' \in {v \in 0..MaxOff :
                    ~GateRelevant \/ AbsDiff(v, truthOff) <= B}
  /\ pc' = "act"
  /\ UNCHANGED <<phase, t, truthOff, near, closing, replans, burnIdx>>

Replan ==                       \* consume the queue head, reset the cursor
  /\ pc = "act" /\ ReplanDue
  /\ replans' = replans - 1 /\ burnIdx' = 0
  /\ UNCHANGED <<phase, t, pc, truthOff, estOff, near, closing>>

ExecuteBurn ==                  \* cursor advances within the current plan
  /\ pc = "act" /\ phase = "A" /\ burnIdx < MaxBurns
  /\ burnIdx' = burnIdx + 1
  /\ UNCHANGED <<phase, t, pc, truthOff, estOff, near, closing, replans>>

EndActA ==                      \* a due replan cannot be skipped (urgent)
  /\ pc = "act" /\ phase = "A" /\ ~ReplanDue
  /\ pc' = "tick"
  /\ UNCHANGED <<phase, t, truthOff, estOff, near, closing, replans, burnIdx>>

GateEvaluate ==                 \* the corridor hold, keyed on the ESTIMATE
  /\ pc = "act" /\ phase = "BC"
  /\ closing' = ~(near /\ estOff > GateTol)
  /\ pc' = "tick"
  /\ UNCHANGED <<phase, t, truthOff, estOff, near, replans, burnIdx>>

AdvancePhase ==                 \* one-way time trigger (urgent)
  /\ pc = "tick" /\ phase = "A" /\ t >= PhaseATicks
  /\ phase' = "BC"
  /\ UNCHANGED <<t, pc, truthOff, estOff, near, closing, replans, burnIdx>>

Tick ==
  /\ pc = "tick" /\ t < HorizonTicks
  /\ ~(phase = "A" /\ t >= PhaseATicks)
  /\ t' = t + 1 /\ pc' = "env"
  /\ UNCHANGED <<phase, truthOff, estOff, near, closing, replans, burnIdx>>

Done == pc = "tick" /\ t = HorizonTicks /\ UNCHANGED vars

Next == EnvDrift \/ Estimate \/ Replan \/ ExecuteBurn \/ EndActA
          \/ GateEvaluate \/ AdvancePhase \/ Tick \/ Done

Spec == Init /\ [][Next]_vars
             /\ WF_vars(EnvDrift) /\ WF_vars(Estimate) /\ WF_vars(Replan)
             /\ WF_vars(EndActA) /\ WF_vars(GateEvaluate)
             /\ WF_vars(AdvancePhase) /\ WF_vars(Tick)

(* --- checked formulas --------------------------------------------------- *)
\* the stated estimator assumption immediately before GateEvaluate
EstimateOK ==
  (pc = "act" /\ GateRelevant) => AbsDiff(estOff, truthOff) <= B

\* the central conditional-safety claim, at the command instant
CorridorHoldSound ==
  (pc = "tick" /\ phase = "BC" /\ near /\ truthOff > Tol) => ~closing

PhaseMono  == [][phase = "BC" => phase' = "BC"]_vars
ReplanMono == [][replans' <= replans]_vars
BurnMono   == [][burnIdx' >= burnIdx \/ replans' # replans]_vars
\* deliberately NOT in Mission.cfg: false, kept for the receipt above
BurnMonoUnconditional == [][burnIdx' >= burnIdx]_vars

HoldReleases ==                 \* the hold cannot livelock past the horizon
  []((~closing /\ phase = "BC") => <>(closing \/ t = HorizonTicks))
ReplanLive ==                   \* queued replans fire during phase A
  <>(replans = 0 \/ phase = "BC")
===============================================================================
