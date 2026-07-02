"""Relative-motion and attitude dynamics models.

Modules
-------
cw          Re-export of the verifiable CW kernel (rpod.core.cw).
th          Tschauner-Hempel equations for eccentric target orbits (Yamanaka-
            Ankersen state transition).
nonlinear   Nonlinear two-body relative dynamics in the target LVLH frame,
            with optional J2 and exponential-atmosphere drag — the LEO/MEO
            truth model for validating linearized guidance.
attitude    Rigid-body rotational dynamics (Euler equations) with reaction
            wheel and thruster torque inputs.
"""

from rpod.core import cw  # noqa: F401
