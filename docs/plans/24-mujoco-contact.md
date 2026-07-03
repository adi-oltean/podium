# 24 — MuJoCo contact/capture backend

GitHub issue: https://github.com/adi-oltean/podium/issues/24

## Fix (landed) — `podium.sim.contact`

Probe-drogue capture truth model, MJCF generated programmatically:
funnel of 8 convex box plates (MuJoCo convexifies meshes — a
non-convex funnel must be a union of convex parts) oriented by
slant/tangent frames, an axis-parallel throat sleeve, backstop disc;
chaser = 500 kg free body with probe capsule and tracked tip site;
zero gravity (seconds-scale event; orbital terms are micrometers);
attitude-hold approximated by large rotational inertia (contact
attitude is deferred 6-DOF scope). Capture = tip seated in the throat
for 0.5 s. Optional sustained docking thrust (probe-drogue practice:
keep pushing through capture).

Contact-modeling lessons paid for and encoded in comments:
1. box-EDGE contact normals have axial components that stop the probe
   regardless of friction — the funnel-throat junction needs a sleeve
   so it is a surface, not an edge;
2. cone steepness dominates ballistic capture: the 47-deg first cut
   reflected most of the axial momentum (probe stalled ON the wall at
   0.1 m/s closing); ~29 deg captures the whole IDSS box;
3. default friction (1.0) wedges probes; slick surfaces (mu = 0.05)
   are both realistic and necessary.

## Receipts (all green)

- centered slow approach captures with zero wall contact;
- IDSS translation-box corners (closing 0.05/0.10, offset 0/0.10,
  lateral rate 0/0.04) ALL capture with 20 N docking thrust, nominal
  corner also fully ballistic — the acceptance box from #17 is now
  tied to a physical capture mechanism;
- 0.45 m offset (beyond mouth) misses; 2 m/s bounces out;
- peak contact force strictly monotone in closing rate;
- bitwise deterministic; envelope boundary bracketed (captures through
  0.2 m, misses at 0.45 m).

## Deferred

Contact attitude / angular misalignment (6-DOF PTR + contact, the one
remaining open v0.4 item), latch/retract modeling, compliant drogue,
dispersed capture-envelope MC campaign through sim.monte_carlo.

## Push/merge instructions

Single commit on main: `24 — MuJoCo contact/capture (#24)`; push; close.
