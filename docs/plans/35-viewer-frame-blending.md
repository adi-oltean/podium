# 35 — three.js viewer ECI/LVLH frame blending

GitHub issue: https://github.com/adi-oltean/podium/issues/35

## Fix (landed)

- `Trace.frame_quaternions()` + `to_viewer_json` q_le field: per-sample
  LVLH->ECI rotation quaternions (w,x,y,z) from the target trajectory
  via Shepperd's method, hemisphere-continuous so the viewer can slerp
  sample-to-sample. Schema is now {meta, t, x, q_le, burns}.
- `viewer/3d/index.html`: a LVLH<->inertial blend slider. LVLH-frame
  structures (station, corridor, grid, keep-out) live under a `world`
  group rotated by the CURRENT-time blended frame; the trajectory
  trail is rotated PER-POINT by slerp(I, QREL[i], blend) so at full
  blend the straight glideslope becomes its true curved inertial-space
  arc (not a rigid spin). QREL = q_le[0]^-1 . q_le[k], conjugated
  through the LVLH->scene axis permutation so it stays a proper
  rotation. Chaser and burn glyphs ride their own instants' frames.
  Embedded DATA regenerated (engine replay of the build-1 scenario)
  to carry q_le.

## Why blend ROTATIONS, not positions

The doc sketch was `render_pos = inertial - blend*offset`. Blending
rotations from identity instead keeps the shipped LVLH view
pixel-identical at blend 0 and makes RANGE provably blend-invariant (a
rotation about the target preserves distance) — the physics on screen
cannot change with the view, only the frame does.

## Receipts

- test_sim.test_frame_quaternions_physics: unit norm, hemisphere
  continuity, conversion matches lvlh_rotation^T on samples, and the
  inter-sample angle equals n*dt (circular orbit).
- test_viewer3d (Playwright): q_le exported; at inertial blend the
  world rotates ~n*t (measured 77.8 deg at t=1200 s — exactly the
  frame rotation), and chaser range is blend-invariant to 1e-3 m.

## Deferred

Chaser recentering, log-time map, URL state; re-DEPLOY of the live
site (tools/deploy_viewer.py mints a new immutable build — a user-side
step, not done here).

## Push/merge instructions

Single commit on main: `35 — Viewer frame blending (#35)`; push;
close. The live site updates on the next viewer deploy.
