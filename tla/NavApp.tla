------------------------------- MODULE NavApp -------------------------------
(* The message protocol of the generated cFS application
   (examples/cfs_nav_app/podium_nav_app.c, emitted by podium.emit.cfsapp).
   The measurement handler is a two-kernel read-modify-write chain on the
   static module state g: podium_ekf_update_sequential, copy-back, then
   podium_ekf_predict, copy-back.  Between the two copy-backs g is
   INCONSISTENT — updated but not propagated — and must not be observable.
   The kernel-level proofs (EVA zero alarms, CompCert, golden vectors) are
   established under sequential semantics; this module states, for the
   first time, the integration-level contract those proofs assume:

     - OnMeas is atomic with respect to the filter state (NoPartialPublish);
     - bus messages are consumed one at a time, in order (InOrder);
     - every measurement is eventually processed (AllProcessed, under
       fairness on the source and the run loop — stated assumptions,
       constraint C4 of docs/add-tla-specs.md).

   In the baseline (ReaderEnabled = FALSE, NavApp.cfg) all three hold —
   atomicity is enforced by structure alone: one task, one pipe, blocking
   receive (podium_nav_app.c run loop).

   Experiment E-NAV (required; NavAppReader.cfg): ReaderEnabled = TRUE
   adds a second task that reads g without going through the pipe.  TLC
   must then produce the violating interleaving of NoPartialPublish —
   Send, Recv, Update, Reader — the concrete witness that the atomicity
   obligation is real, not hypothetical.  The CI lane asserts this run
   FAILS (a falsification receipt, the practice of tests/test_arch.py);
   it mirrors what embedding the kernels in a multi-task flight build
   would do without a locking or message-ownership discipline.

   pc \in {"idle", "have", "updated"}: "have" is after dispatch, before
   the update chain (g still consistent, one cycle stale); "updated" is
   the inconsistent interior between the copy-backs. *)
EXTENDS Naturals, Sequences

CONSTANTS
  NMsgs,          \* measurements the environment will send (in order)
  QCap,           \* bus pipe capacity (16 in the generated app)
  ReaderEnabled   \* experiment E-NAV: enable the un-piped reader

ASSUME NMsgs \in Nat /\ QCap \in Nat \ {0} /\ ReaderEnabled \in BOOLEAN

VARIABLES
  sent,       \* measurements sent so far (environment)
  queue,      \* the software-bus pipe, FIFO
  pc,         \* the handler's program counter (see above)
  current,    \* the dequeued message being processed (0 = none)
  processed,  \* sequence of published STATE messages
  dirtyRead   \* set iff a reader ever observed g mid-chain (E-NAV)

vars == <<sent, queue, pc, current, processed, dirtyRead>>

TypeOK ==
  /\ sent \in 0..NMsgs
  /\ queue \in Seq(1..NMsgs) /\ Len(queue) <= QCap
  /\ pc \in {"idle", "have", "updated"}
  /\ current \in 0..NMsgs
  /\ processed \in Seq(1..NMsgs)
  /\ dirtyRead \in BOOLEAN

Init == sent = 0 /\ queue = <<>> /\ pc = "idle" /\ current = 0
        /\ processed = <<>> /\ dirtyRead = FALSE

Send ==      \* environment: measurements arrive in order, bus-bounded
  /\ sent < NMsgs /\ Len(queue) < QCap
  /\ sent' = sent + 1 /\ queue' = Append(queue, sent + 1)
  /\ UNCHANGED <<pc, current, processed, dirtyRead>>

Recv ==      \* CFE_SB_ReceiveBuffer + MsgId dispatch: consume the head
  /\ pc = "idle" /\ queue # <<>>
  /\ current' = Head(queue) /\ queue' = Tail(queue) /\ pc' = "have"
  /\ UNCHANGED <<sent, processed, dirtyRead>>

Update ==    \* podium_ekf_update_sequential + copy-back: g now mid-chain
  /\ pc = "have" /\ pc' = "updated"
  /\ UNCHANGED <<sent, queue, current, processed, dirtyRead>>

Predict ==   \* podium_ekf_predict + copy-back + STATE publish
  /\ pc = "updated" /\ pc' = "idle"
  /\ processed' = Append(processed, current)
  /\ UNCHANGED <<sent, queue, current, dirtyRead>>

Reader ==    \* E-NAV only: a second task reads g, bypassing the pipe
  /\ ReaderEnabled
  /\ dirtyRead' = (dirtyRead \/ pc = "updated")
  /\ UNCHANGED <<sent, queue, pc, current, processed>>

Next == Send \/ Recv \/ Update \/ Predict \/ Reader

(* Fairness: the source delivers its NMsgs messages and the run loop
   keeps executing (WF on Send/Recv/Update/Predict).  The reader gets no
   fairness — it needs only one lucky moment, and TLC finds it. *)
Spec == Init /\ [][Next]_vars
             /\ WF_vars(Send) /\ WF_vars(Recv)
             /\ WF_vars(Update) /\ WF_vars(Predict)

(* --- checked formulas --------------------------------------------------- *)
InOrder          == \A i \in 1..Len(processed) : processed[i] = i
NoPartialPublish == ~dirtyRead
AllProcessed     == <>(Len(processed) = NMsgs)
===============================================================================
