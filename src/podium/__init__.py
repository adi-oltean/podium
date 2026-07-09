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
emit        C99 code generator for the static subset: ACSL contracts, the EVA
            driver, and the cFS app, checked by golden vectors and CompCert.
verify      Exact-rational certificate checkers (barrier, KKT, optimality-gap,
            Lyapunov, SOS) and input-range contracts.
"""

__version__ = "0.8.0.dev0"
