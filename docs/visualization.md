# Visualization design

Two tiers: matplotlib for analysis, a zero-build three.js viewer for
interactive 3-D. The interactive tier adopts patterns proven in the fermi
project (interstellar-mission explorer with a real-time three.js chase cam).

## Tier 1 — analysis plots (matplotlib)

- LVLH trajectory projections (y-x "V-bar view", z-x "cross-track view") with
  approach corridor and keep-out-zone overlays;
- range / range-rate corridor plots against the approach envelope;
- delta-v timelines and cumulative budget;
- Monte Carlo dispersion: trajectory fans, capture-envelope scatter at
  contact, violation histograms.

## Tier 2 — interactive viewer (static HTML + vendored three.js)

A single `viewer/index.html` plus a vendored `three.min.js`; no bundler, no
node toolchain. Simulations export a compact JSON (`sim.to_viewer_json()`)
with decimated trajectory, attitude keyframes, burn events, and geometry
descriptors; the viewer loads it via a file picker or URL parameter.

Patterns adopted from fermi (and why):

1. **Render loop decoupled from the sim clock.** A continuous
   `requestAnimationFrame` loop that only renders, and a separate fixed-cadence
   timer that advances playback time and updates geometry buffers. Camera
   stays smooth during pause and drag.
2. **State → geometry → draw as pure stages with preallocated buffers.**
   Trail polylines live in fixed-size `Float32Array`s mutated in place
   (`needsUpdate` flags), never reallocated per frame; decimation always
   keeps the newest vertex so the trail head never lags.
3. **Continuous reference-frame blending.** RPOD viewing constantly switches
   between inertial (ECI) and target-LVLH perspectives. Instead of hard camera
   cuts, blend: `render_pos = inertial − blend · frame_offset(t)` with a
   smoothstep on `blend`, giving one continuous dolly between "orbit view" and
   "approach view".
4. **Scene recentered on the chaser** so 32-bit GPU coordinates stay small —
   the sim spans ~7000 km orbit radius down to centimeter docking offsets.
5. **Analytic chase camera from the LVLH basis** (position = offset along
   −approach axis tilted up, up = orbit normal, lookAt target), with
   auto-framing from the visible-trail bounding box and per-frame near/far
   recomputation to cover the km→cm dynamic range. Free-orbit mode seeded
   from the current pose so mode switches don't jump.
6. **Single normalized time parameter** `u ∈ [0,1]` with a pluggable
   linear/log time map — RPOD timelines mix multi-hour phasing with
   seconds-scale contact dynamics; log-time keeps both scrubbably visible.
   The timeline doubles as a scrubber.
7. **Shareable state via URL parameters** (scenario, camera, playback time)
   for cheap reproducibility.

RPOD-specific scene elements: target model with docking axis and capture
envelope cone, approach corridor (translucent cone/frustum), keep-out sphere,
V-bar/R-bar grid in the LVLH plane, thruster-firing glyphs at burn events,
day/night terminator lighting from the sun vector.

## Parity discipline

If any dynamics ever get reimplemented viewer-side (e.g., smooth interpolation
between exported states), fermi's rule applies: the Python engine is the
source of truth, the JS side is a fast approximation, and an automated parity
check pins them together. Default posture: the viewer *only interpolates
exported data* and computes no physics.
