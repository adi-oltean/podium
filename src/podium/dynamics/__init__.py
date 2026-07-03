"""Relative-motion and attitude dynamics models.

Modules
-------
cw          Re-export of the verifiable CW kernel (podium.core.cw).
ya          Re-export of the Yamanaka-Ankersen/Tschauner-Hempel state
            transition for eccentric target orbits (podium.core.ya).
roe         Re-export of the relative-orbital-elements kernel
            (podium.core.roe): quasi-nonsingular ROE, Koenig Keplerian/J2
            STMs, near-circular LVLH maps, impulsive control matrix.
nonlinear   Nonlinear two-body relative dynamics in the target LVLH frame,
            with optional J2 and exponential-atmosphere drag — the LEO/MEO
            truth model for validating linearized guidance.
attitude    Rigid-body rotational dynamics (Euler equations, quaternion
            kinematics, RK4) with body-frame torque input; conservation
            receipts in tests/test_attitude.py.
"""

from podium.core import cw, roe, ya  # noqa: F401
from podium.dynamics import attitude, nonlinear  # noqa: F401
