"""Attitude receipts: torque-free conservation (energy, inertial angular
momentum, unit norm) incl. the intermediate-axis tumble; detumble;
pointing step response pinned to the second-order design prediction."""

import math

import numpy as np

from podium.control import attitude as ac
from podium.core import quat
from podium.dynamics import attitude as att

INERTIA = np.diag([120.0, 90.0, 60.0])  # small-sat class, distinct axes


def propagate(q, w, torque_fn, dt, steps):
    qs = [q.copy()]
    ws = [w.copy()]
    for k in range(steps):
        tau = torque_fn(k * dt, qs[-1], ws[-1])
        q_n, w_n = att.step(qs[-1], ws[-1], INERTIA, tau, dt)
        qs.append(q_n)
        ws.append(w_n)
    return np.array(qs), np.array(ws)


def test_torque_free_conservation_intermediate_axis():
    """Tumble about the intermediate axis (the hard case: unstable,
    trajectory explores the full polhode): energy and inertial angular
    momentum conserved to 1e-9 relative over 10 minutes, |q| = 1."""
    q0 = quat.normalize(np.array([0.9, 0.2, -0.3, 0.1]))
    w0 = np.array([0.001, 0.35, 0.002])  # intermediate axis + perturbation
    qs, ws = propagate(q0, w0, lambda t, q, w: np.zeros(3), dt=0.05,
                       steps=12_000)
    e0 = att.kinetic_energy(ws[0], INERTIA)
    l0 = att.momentum_inertial(qs[0], ws[0], INERTIA)
    for idx in (3000, 6000, 12_000):
        e = att.kinetic_energy(ws[idx], INERTIA)
        l_vec = att.momentum_inertial(qs[idx], ws[idx], INERTIA)
        assert abs(e - e0) / e0 < 1e-9
        assert np.linalg.norm(l_vec - l0) / np.linalg.norm(l0) < 1e-9
        assert abs(np.linalg.norm(qs[idx]) - 1.0) < 1e-12
    # the tumble genuinely departed the initial axis (instability seen)
    assert float(np.max(np.abs(ws[:, 0]))) > 0.05


def test_detumble():
    """10 deg/s tumble driven below 0.01 deg/s with saturated torques."""
    q0 = quat.identity()
    w0 = np.radians([10.0, -7.0, 4.0])
    q_ref = quat.identity()

    def ctl(t, q, w):
        return ac.quaternion_feedback(q, w, q_ref, np.zeros(3),
                                      kp=2.0, kd=60.0, tau_max=1.0)

    qs, ws = propagate(q0, w0, ctl, dt=0.1, steps=6000)
    assert float(np.linalg.norm(ws[-1])) < math.radians(0.01)
    assert float(np.linalg.norm(quat.error(qs[-1], q_ref))) < 1e-3


def test_step_response_matches_second_order_design():
    """20-degree single-axis slew about x (J = 120): gains chosen for
    wn = 0.1 rad/s, zeta = 0.9; the response must show the predicted
    settling and near-zero overshoot of a zeta = 0.9 system."""
    j = float(INERTIA[0, 0])
    wn, zeta = 0.1, 0.9
    kp = wn * wn * j
    kd = 2.0 * zeta * wn * j
    ang = math.radians(20.0)
    q0 = np.array([math.cos(ang / 2), math.sin(ang / 2), 0.0, 0.0])
    q_ref = quat.identity()

    def ctl(t, q, w):
        return ac.quaternion_feedback(q, w, q_ref, np.zeros(3),
                                      kp=kp, kd=kd, tau_max=100.0)

    dt = 0.05
    qs, ws = propagate(q0, np.zeros(3), ctl, dt=dt, steps=3000)
    errs = np.array([quat.error(q, q_ref)[0] for q in qs])  # ~2*half-angle
    # overshoot of a zeta=0.9 second-order step: exp(-pi zeta/sqrt(1-z^2))
    overshoot = float(-np.min(errs)) / float(errs[0])
    pred = math.exp(-math.pi * zeta / math.sqrt(1 - zeta * zeta))
    assert overshoot < pred + 0.02  # quaternion factor-2 nonlinearity helps
    # 2% settling time ~ 4/(zeta wn) = 44 s
    t_settle_pred = 4.0 / (zeta * wn)
    idx = int(1.5 * t_settle_pred / dt)
    assert abs(errs[idx]) < 0.02 * abs(errs[0])
    # symmetric: no torque exceeded the (generous) limit, rates bounded
    assert float(np.max(np.abs(ws))) < 0.1


def test_saturation_respected():
    q0 = quat.normalize(np.array([0.7, 0.7, 0.1, 0.0]))
    tau = ac.quaternion_feedback(q0, np.array([0.5, -0.5, 0.2]),
                                 quat.identity(), np.zeros(3),
                                 kp=100.0, kd=100.0, tau_max=0.25)
    assert float(np.max(np.abs(tau))) <= 0.25 + 1e-15