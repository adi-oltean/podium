"""Aerodynamic disturbance-torque receipts (#45): the second dominant
LEO attitude disturbance after gravity gradient, validated against the
analytic tau = r_cp x F_drag and the weathervane-stability behavior.

The drag force is consistent with the truth model
(podium.dynamics.nonlinear): with cd_area = m/bc, F_aero equals
m * a_drag there.
"""

import numpy as np

from podium.dynamics import attitude as att
from podium.dynamics import nonlinear as nl


def test_torque_equals_cp_cross_drag():
    """tau = r_cp x F_drag exactly, F_drag = -1/2 rho Cd A |v| v."""
    v = np.array([7500.0, 20.0, -5.0])
    rho, cd_area = 3.0e-12, 4.4
    r_cp = np.array([-0.5, 0.1, 0.05])
    tau = att.aerodynamic_torque(v, rho, cd_area, r_cp)
    f = -0.5 * rho * cd_area * np.linalg.norm(v) * v
    assert np.allclose(tau, np.cross(r_cp, f), atol=1e-18)


def test_zero_when_cp_aligned_with_flow():
    """No moment arm perpendicular to the drag: cp on the velocity line
    gives zero torque."""
    v = np.array([7500.0, 0.0, 0.0])
    r_cp = np.array([-0.8, 0.0, 0.0])       # along -v
    tau = att.aerodynamic_torque(v, 3e-12, 4.4, r_cp)
    assert np.allclose(tau, 0.0, atol=1e-18)


def test_weathervane_restoring_sign():
    """Center of pressure BEHIND the c.m. (downstream): a small angle
    of attack produces a torque that rotates the body BACK toward the
    flow — passive aero (weathervane) stability."""
    speed = 7500.0
    cd_area, rho = 4.4, 3.0e-12
    r_cp = np.array([-1.0, 0.0, 0.0])       # cp behind c.m. along body -x
    # body flying +x into the wind, tilted by +alpha about body z:
    # relative wind in body frame gains a -y component
    alpha = 0.02
    v_body = speed * np.array([np.cos(alpha), -np.sin(alpha), 0.0])
    tau = att.aerodynamic_torque(v_body, rho, cd_area, r_cp)
    # torque about +z must OPPOSE the +alpha tilt (restore to alignment)
    assert tau[2] < 0.0
    # cp AHEAD of c.m. is the unstable (destabilizing) configuration
    tau_unstable = att.aerodynamic_torque(v_body, rho, cd_area, -r_cp)
    assert tau_unstable[2] > 0.0


def test_force_consistent_with_truth_model_drag():
    """The aero force implied by aerodynamic_torque matches the truth
    model's drag acceleration times mass (cd_area = m/bc)."""
    mass, bc = 500.0, 100.0            # bc = m/(Cd A) -> Cd A = m/bc = 5
    cd_area = mass / bc
    r = np.array([6.8e6, 0.0, 0.0])
    v = np.array([0.0, 7600.0, 10.0])
    cfg = nl.ForceConfig(j2=0.0, drag=nl.DragConfig(rho0=3e-12, h0=400e3,
                                                    scale_height=60e3))
    a_drag = nl.perturb_accel(r, v, cfg, bc)      # ECI drag accel
    rho = cfg.drag.density(float(np.linalg.norm(r)) - cfg.r_body)
    v_rel = v - np.array([-cfg.omega_earth * r[1],
                          cfg.omega_earth * r[0], 0.0])
    f_from_torque = -0.5 * rho * cd_area * np.linalg.norm(v_rel) * v_rel
    assert np.allclose(mass * a_drag, f_from_torque, rtol=1e-12)
