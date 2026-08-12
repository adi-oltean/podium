---------------------------- MODULE ArchRendezvous ----------------------------
(* Discrete abstraction of the ARCH rendezvous mode logic in
   src/podium/guidance/arch.py (simulate(), urgent guard semantics, and
   the transition structure of export_model()).  Box entry is a
   nondeterministic environment event; the continuous semantics — whether
   the trajectory actually enters the attempt octagon, LOS cone / velocity
   / keep-out safety — are owned by the reachability lane (tools/reach).
   Time unit: minutes, as in the benchmark.

   Urgency is encoded in Tick's guard: time cannot advance while an urgent
   transition is enabled.  This turns the bounded-liveness claim "abort is
   engaged by the deadline" into the state invariant AbortByDeadline — the
   clock-augmentation idiom the code base already uses (arch.py keeps the
   clock as a state variable so aborts are time-triggered guards).

   Checked configurations (see ../.github/workflows/tla.yml):
     ArchRendezvousSRA.cfg   AbortTime = 120  (abort scenario, SR-A01)
     ArchRendezvousSRNA.cfg  AbortTime = 301  (> Horizon: no abort, SR-NA01)

   Falsification receipt (kept as a comment, per the practice of
   tests/test_arch.py): the unconditional EntryFromBox is VIOLATED in the
   SRA configuration — TLC exhibits the trace where EnvEnterBox fires at
   clock = AbortTime, EnterAttempt is already disabled (~AbortEnabled
   fails), and Abort takes over: box entry that races the abort deadline
   in the same instant never yields attempt mode.  The fixed-dt Python
   simulation quantizes exactly this simultaneity away (arch.py:
   "switching times are quantized to dt").  The SRA configuration
   therefore checks the corrected conditional property EntryBeforeAbort;
   EntryFromBox is checked where it is sound, in SRNA. *)
EXTENDS Naturals

CONSTANTS AbortTime, Horizon, BoxLatched
ASSUME AbortTime \in Nat /\ Horizon \in Nat \ {0} /\ BoxLatched \in BOOLEAN
  (* "no abort" (SRNA) is modeled as AbortTime > Horizon.  BoxLatched =
     TRUE is the baseline (entry is one-way, matching the hybrid model
     whose attempt-mode invariant confines the state to the box);
     BoxLatched = FALSE is the box-exit robustness experiment, checked
     by ArchRendezvousExit.cfg — see EnvExitBox below. *)

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

Abort ==                        \* urgent, dominant (arch.py abort branch)
  /\ AbortEnabled /\ mode # "aborting"
  /\ mode' = "aborting" /\ UNCHANGED <<clock, inBox>>

EnterAttempt ==                 \* urgent; dominated by Abort
  /\ mode = "approaching" /\ inBox /\ ~AbortEnabled
  /\ mode' = "attempt" /\ UNCHANGED <<clock, inBox>>

EnvEnterBox ==                  \* environment event
  /\ mode = "approaching" /\ ~inBox
  /\ inBox' = TRUE /\ UNCHANGED <<mode, clock>>

EnvExitBox ==                   \* robustness experiment only (see below)
  /\ ~BoxLatched /\ mode = "approaching" /\ inBox
  /\ inBox' = FALSE /\ UNCHANGED <<mode, clock>>

Tick ==                         \* time advances only when nothing is urgent
  /\ clock < Horizon
  /\ ~(AbortEnabled /\ mode # "aborting")
  /\ ~(mode = "approaching" /\ inBox /\ ~AbortEnabled)
  /\ clock' = clock + 1 /\ UNCHANGED <<mode, inBox>>

Done == clock = Horizon /\ UNCHANGED vars

Next == Abort \/ EnterAttempt \/ EnvEnterBox \/ EnvExitBox \/ Tick \/ Done

(* EnvEnterBox/EnvExitBox carry no fairness: the environment is never
   assumed helpful.  Fairness on Abort/EnterAttempt/Tick is the
   controller's and scheduler's obligation only. *)
Spec == Init /\ [][Next]_vars
             /\ WF_vars(Abort) /\ WF_vars(EnterAttempt) /\ WF_vars(Tick)

(* Box-exit experiment (BoxLatched = FALSE): weak fairness is no longer
   enough for entry — the box can flicker (enter/exit alternating with
   the clock frozen by urgency), so EnterAttempt is enabled infinitely
   often but never continuously, and WF never obliges it to fire.  TLC
   exhibits the flicker trace if EntrySustainedBox is checked under
   Spec.  Under strong fairness (SpecSF) sustained presence converts to
   mode entry: EntrySustainedBox holds.  Plain EntryFromBox fails under
   ANY fairness once exits are allowed — entry can be revoked before the
   controller reacts, and if the box is never re-entered no fairness
   assumption applies.  Both falsification receipts are in
   docs/add-tla-specs.md Section 10. *)
SpecSF == Spec /\ SF_vars(EnterAttempt)

(* --- checked formulas --------------------------------------------------- *)
AbortByDeadline    == clock > AbortTime => mode = "aborting"
AttemptRequiresBox == mode = "attempt" => inBox

(* Abort dominance must be stated at the ACTION level: negative control
   N1 (drop ~AbortEnabled from EnterAttempt) leaves AbortByDeadline and
   Monotone intact — urgency still stops the clock until Abort fires, so
   the swapped priority only inserts a transient attempt step.  This
   formula is the one N1 falsifies. *)
AbortDominates == [][~(AbortEnabled /\ mode' = "attempt")]_vars
Monotone  == [][ModeOrd[mode'] >= ModeOrd[mode]]_vars
Absorbing == [][mode = "aborting" => mode' = "aborting"]_vars
AbortTaken   == <>[](mode = "aborting")            \* SRA config only

(* Conditional liveness: IF the continuous layer delivers box entry
   (reachability/simulation evidence supports it), the discrete layer
   converts it to mode entry.  The unconditional form is sound only when
   no abort deadline exists (SRNA); see the falsification receipt above. *)
EntryFromBox    == (<> inBox) => <>(mode = "attempt")
EntryBeforeAbort == <>(inBox /\ clock < AbortTime) => <>(mode = "attempt")
EntrySustainedBox == ([]<> inBox) => <>(mode = "attempt")  \* Exit config
===============================================================================
