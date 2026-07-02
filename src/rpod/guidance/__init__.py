"""Guidance laws and trajectory optimization for proximity operations.

Modules
-------
glideslope  Classical glideslope approach guidance (Hablani), core-compliant.
targeting   Multi-impulse CW/TH targeting wrappers around rpod.core.cw.
convex      Convex trajectory optimization: direct SOCP transcription of the
            CW/TH relative dynamics with approach-cone, keep-out-zone (via
            successive convexification), thrust, and plume constraints.
            Prototyping backend: cvxpy; embedded path: generated fixed-
            iteration solver (see docs/trajectory-optimization.md).
safety      Passive-abort / free-drift safety checks for candidate
            trajectories.
"""
