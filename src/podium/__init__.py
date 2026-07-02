"""Podium: physics-precise RPOD GNC library and simulation sandbox for LEO/MEO.

Subpackages
-----------
core        Verifiable algorithm core: pure, statically-shaped step functions
            (the subset intended for Python->C translation and abstract
            interpretation).
dynamics    Relative-motion and attitude dynamics models (CW, Tschauner-Hempel,
            nonlinear two-body, J2/drag perturbations).
guidance    Guidance laws and trajectory optimization (glideslope, CW targeting,
            convex/SCP transcription).
control     Feedback controllers and actuator allocation (LQR, MPC hooks,
            quaternion attitude control, thruster mapping).
nav         Relative navigation filters and sensor models.
sim         Deterministic fixed-step simulation engine, events, Monte Carlo.
viz         Visualization of trajectories, approach corridors, and sim state.
verify      Contracts, input-range assumptions, and export harness for the
            external abstract-interpretation validation tool.
"""

__version__ = "0.0.1"
