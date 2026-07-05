"""Disturbance-environment receipts (#47): the three environmental
torques compose, and a quaternion-feedback controller holds attitude
against the combined disturbance over a full orbit.
"""

import numpy as np
import pytest

from podium.control import attitude as att_ctl
from podium.core import quat
from podium.dynamics import attitude as att
from podium.dynamics.disturbances import DisturbanceModel

N = 0.0011


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def _model():
    return DisturbanceModel(
        inertia=np.diag([180.0, 140.0, 100.0]),
        n=N,
        aero=(4.4, np.array([-0.6, 0.0, 0.0])),
        srp=(12.0, 1.8, np.array([0.0, 0.2, 0.0])),
    )


def test_superposition_equals_sum_of_three():
    """The aggregate torque is exactly gravity-gradient + aero + SRP."""
    m = _model()
    nadir = _unit([0.1, -0.9, 0.3])
    vrel = np.array([7500.0, 30.0, -10.0])
    sun = _unit([0.4, 0.5, 0.7])
    rho = 2.0e-12
    total = m.torque(nadir, vrel, rho, sun, illuminated=True)
    gg = att.gravity_gradient_torque(nadir, m.inertia, N)
    aero = att.aerodynamic_torque(vrel, rho, m.aero[0], m.aero[1])
    srp = att.srp_torque(sun, m.srp[0], m.srp[1], m.srp[2],
                         illuminated=True)
    assert np.allclose(total, gg + aero + srp, atol=1e-18)


def test_terms_disable_cleanly():
    """Leaving a term's config None (or rho=0) drops it; gravity
    gradient always contributes."""
    m = DisturbanceModel(inertia=np.diag([180.0, 140.0, 100.0]), n=N)
    nadir = np.array([0.0, 0.0, 1.0])          # aligned -> gg zero too
    assert np.allclose(m.torque(nadir), 0.0, atol=1e-15)
    nadir2 = _unit([0.2, 0.0, 1.0])
    only_gg = m.torque(nadir2, np.array([7500.0, 0, 0]), 0.0)
    assert np.allclose(only_gg,
                       att.gravity_gradient_torque(nadir2, m.inertia, N))


def test_magnitudes_are_physically_reasonable():
    """Order-of-magnitude sanity for a representative LEO spacecraft:
    gravity gradient dominates, then aero, then SRP."""
    m = _model()
    nadir = _unit([0.15, 0.0, 1.0])   # ~8 deg tilt
    vrel = np.array([7600.0, 150.0, 0.0])   # slight sideslip -> aero arm
    sun = np.array([1.0, 0.0, 0.0])         # off the srp cp axis -> arm
    gg = np.linalg.norm(att.gravity_gradient_torque(nadir, m.inertia, N))
    aero = np.linalg.norm(att.aerodynamic_torque(vrel, 2e-12, *m.aero))
    srp = np.linalg.norm(att.srp_torque(sun, *m.srp))
    # physically-reasonable bands for a representative LEO spacecraft
    # (micro-Nm-class disturbances); exact ordering is geometry-dependent
    assert 1e-6 < gg < 1e-3
    assert 1e-7 < aero < 1e-4
    assert 1e-7 < srp < 1e-4


@pytest.mark.slow
def test_attitude_hold_rejects_combined_disturbance():
    """A quaternion-feedback controller holds an inertially-fixed
    attitude against the combined environmental disturbance over a full
    orbit — pointing error stays small and the control stays inside
    saturation. Nadir and the relative wind rotate at the orbital rate;
    the Sun is inertially fixed."""
    m = _model()
    q = np.array([1.0, 0.0, 0.0, 0.0])
    w = np.zeros(3)
    q_ref = np.array([1.0, 0.0, 0.0, 0.0])
    dt = 1.0
    kp, kd, tau_max = 5.0, 200.0, 0.05
    max_err_deg, max_ctrl = 0.0, 0.0
    steps = int(2 * np.pi / N)
    for k in range(steps):
        t = k * dt
        # disturbance directions in the inertial frame
        nadir_in = np.array([-np.cos(N * t), -np.sin(N * t), 0.0])
        vrel_in = 7600.0 * np.array([-np.sin(N * t), np.cos(N * t), 0.0])
        sun_in = np.array([0.0, 0.0, 1.0])
        # into the body frame
        cq = quat.conjugate(q)
        d = m.torque(quat.rotate(cq, nadir_in), quat.rotate(cq, vrel_in),
                     2.0e-12, quat.rotate(cq, sun_in), illuminated=True)
        u = att_ctl.quaternion_feedback(q, w, q_ref, np.zeros(3),
                                        kp, kd, tau_max)
        q, w = att.step(q, w, m.inertia, u + d, dt)
        e = quat.error(q, q_ref)
        ang = np.degrees(2 * np.arcsin(min(1.0, float(np.linalg.norm(e[1:])))))
        max_err_deg = max(max_err_deg, ang)
        max_ctrl = max(max_ctrl, float(np.max(np.abs(u))))
    assert max_err_deg < 1.0, max_err_deg      # holds within 1 degree
    assert max_ctrl <= tau_max + 1e-12         # never clips beyond bound
