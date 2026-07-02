"""Visualization of RPOD trajectories and simulation state.

Two tiers (design borrowed from the fermi project — see docs/visualization.md):

- **Analysis plots** (matplotlib): LVLH trajectory projections, range/range-
  rate corridors, delta-v timelines, Monte Carlo dispersion ellipses.
- **Interactive 3-D viewer** (static HTML + vendored three.js, no build
  system): chase-camera view of the approach with smooth blending between
  inertial and target-LVLH camera frames, scene recentered on the chaser for
  GPU float precision, preallocated trail buffers, and a scrubbable timeline
  decoupled from the render loop. Simulation results export to a compact JSON
  the viewer loads directly.
"""
