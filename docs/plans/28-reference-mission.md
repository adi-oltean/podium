# 28 — End-to-end reference mission + audit bundle

GitHub issue: https://github.com/adi-oltean/podium/issues/28

## Fix (landed) — `podium.sim.mission`

One seeded engine run composing every layer built this cycle:

- START: inside the barrier-certified passively-safe formation
  (podium.verify.barrier SAFE_CASE); dispersed starts are proven to
  stay inside the certified ellipsoid by test, and the certificate is
  re-verified in exact rational arithmetic INTO the audit bundle.
- PHASE A (2 km -> 300 m): PtrDockingPlanner (true keep-out +
  exact-flow CTCS cuts) with mid-course REPLANS at T/3 and 2T/3 from
  the EKF estimate — open-loop pulse plans provably miss under
  actuator error (test_sensors), so the mission replans instead.
- PHASE B/C (300 m -> contact): rate-command terminal feedback with
  the corridor GATE (inside 40 m, closing holds until lateral error is
  in-box — do not advance misaligned).
- PHASE D: the contact-crossing state feeds the MuJoCo probe-drogue
  sim (20 N docking thrust); quaternion-feedback attitude runs in
  parallel and is judged against the rotational IDSS box at contact.
- Nav: Joseph EKF on docking-grade position measurements; actuators:
  2 mm/s MIB + 1% execution error.

## Two measured GNC lessons (both encoded as comments + receipts)

1. **Process noise must budget actuator mismatch.** q_accel = 5e-8
   made the filter near-open-loop (tiny gains); the per-tick
   commanded-vs-applied difference (quantization + execution error)
   integrated into a ~0.1 m lateral estimate bias, and the
   ESTIMATE-keyed corridor gate never saw the true offset — contact at
   0.129 m lateral, outside the box, diagnosed from the endgame trace.
   q_accel = 2e-6 (the known ~1 mm/s/tick disturbance) fixes it.
2. **Docking boxes need docking sensing.** The 0.10 m IDSS lateral box
   is unreachable on 2 m GNSS-class measurements; the mission uses
   fused docking-sensor-grade noise (5 cm), documented.

## Receipts (all green)

Reference run: captured, all IDSS translation+rotation margins
positive, phase STL margins positive, dv 12.4 m/s (measured; terminal
chatter trim noted as polish, bound set at the honest 15), barrier
verdict true. Audit bundle byte-identical across two runs of the same
seed. Dispersed campaign (3 seeds through the full mission): 100%
capture. Initial-state containment in the certified set (20 seeds).

## Deferred

Per-tag release workflow publishing the bundle; dv chatter deadband;
ROE far-range phasing prologue; MC campaign size in CI.

## Push/merge instructions

Single commit on main: `28 — Reference mission + audit bundle (#28)`;
push; close.
