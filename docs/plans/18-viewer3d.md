# 18 — Interactive three.js viewer

GitHub issue: https://github.com/adi-oltean/podium/issues/18

## Problem

Tier-2 visualization (docs/visualization.md) was design-only; the last
open v0.3 item.

## Fix (landed)

`viewer/3d/index.html` + vendored three.js r172 (MIT, two-file module
build under `viewer/3d/vendor/`, version-pinned, no external requests):
LVLH scene (station body + gold docking ring, 200 m KOZ wireframe,
10-degree approach-corridor cone, V-bar plane grid), trajectory trail
in a preallocated Float32Array with drawRange playback (fermi pattern:
rAF render loop decoupled from a fixed-cadence playback timer), chaser
marker, burn glyphs appearing at burn times, follow camera from the
dock-relative geometry plus free-orbit mode seeded from the current
pose (no jump on switch), play/speed/scrub/HUD sharing the 2-D viewer's
embedded DATA schema. Debug handle `window.__P3D` for the UI tests.

Receipts (tools/ui/test_viewer3d.py, fermi RAM hygiene: one browser,
one server, both torn down in finally): zero console errors; WebGL
draws (draw calls > 0 AND non-blank pixel readback — required adding
preserveDrawingBuffer); playback advances; scrub-to-end matches the
physics (range 20.2 m, dv 3.765 shown); camera toggle. The pixel
check caught a real bug: setT without clamping extrapolated past the
trajectory end (range 2.4e7 m) — fixed with clamping in setT.

Deferred (recorded in docs/visualization.md status note): ECI/LVLH
frame blending (needs target-ECI in to_viewer_json), chaser
recentering for large scenes, log-time map, URL-shareable state,
attitude keyframes.

## Acceptance Criteria

- [x] Playwright receipts ALL PASS locally
- [x] Pages deploy covers viewer/3d/ (site root = viewer/)
- [x] README demos, visualization.md status, roadmap refreshed;
      v0.3 marked COMPLETE

## Push/merge instructions

Single commit on main: `18 — three.js 3-D viewer (#18)`; push (triggers
pages deploy); verify the live URL; close.
