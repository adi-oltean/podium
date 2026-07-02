"""Feedback control and actuator allocation.

Modules
-------
lqr        Discrete LQR for CW translational station-keeping/approach.
           Gain synthesis (Riccati iteration) runs offline in full Python;
           the flight-side gain application is core-compliant.
attitude   Quaternion-feedback attitude control (PD with rate limiting).
allocation Thruster mapping: minimum-effort allocation of commanded
           force/torque to a fixed thruster configuration.
"""
