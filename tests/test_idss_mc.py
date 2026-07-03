"""IDSS acceptance + Monte Carlo receipts: a terminal approach flown
through the engine with noisy sensing and an imperfect actuator must
contact inside the IDSS Rev G box; the attitude loop holds the
rotational box at contact time; seeded campaigns reproduce exactly."""

import math

import numpy as np
import pytest

from podium import constants as const
from podium.control import attitude as ac
from podium.core import quat
from podium.dynamics import attitude as att
from podium.sim import Scenario, circular_target, idss, monte_carlo, run

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)
DOCK = np.zeros(3)
AXIS = np.array([0.0, -1.0, 0.0])  # chaser approaches from -y (V-bar)
BOX = idss.IdssBox()


def terminal_controller(v_contact: float = 0.075):
    """Per-axis rate-command terminal approach: closing profile that
    tapers to v_contact, lateral rate commands toward the axis."""

    def ctl(_t, meas):
        y = float(meas[1])
        rng_to_port = abs(y)
        v_close_des = min(0.5, max(v_contact, rng_to_port / 90.0))
        vy_des = v_close_des if y < 0 else 0.0
        vx_des = -float(meas[0]) / 20.0
        vz_des = -float(meas[2]) / 20.0
        vx_des = max(-0.03, min(0.03, vx_des))
        vz_des = max(-0.03, min(0.03, vz_des))
        dv = np.array([vx_des - meas[3], vy_des - meas[4], vz_des - meas[5]])
        step_cap = 0.05
        mag = float(np.linalg.norm(dv))
        if mag > step_cap:
            dv *= step_cap / mag
        return dv

    return ctl


def docking_scenario(x_rel0, seed):
    return Scenario(
        duration=900.0,
        rv_target0=circular_target(A),
        x_rel0=x_rel0,
        dt_gnc=0.5,
        truth_substeps=2,
        seed=seed,
        meas_pos_std=0.02,
        meas_vel_std=0.003,
        dv_quantum=0.001,
        dv_exec_std_frac=0.01,
    )


def contact_metrics(trace):
    """State at the contact crossing (first tick with y >= -0.05 m)."""
    ys = trace.x_rel[:, 1]
    hits = np.flatnonzero(ys >= -0.05)
    if len(hits) == 0:
        return {"contact": 0.0, "closing_above_min": -1.0,
                "closing_below_max": -1.0, "lateral_offset": -1.0,
                "lateral_rate": -1.0}
    k = int(hits[0])
    m = idss.check_translation(trace.x_rel[k], DOCK, AXIS, BOX)
    m["contact"] = 1.0
    return m


@pytest.mark.slow
def test_idss_translation_box_single_run():
    tr = run(docking_scenario(np.array([0.3, -30.0, -0.2, 0.0, 0.0, 0.0]), 3),
             terminal_controller())
    m = contact_metrics(tr)
    assert m["contact"] == 1.0
    for name in ("closing_above_min", "closing_below_max",
                 "lateral_offset", "lateral_rate"):
        assert m[name] > 0.0, (name, m)


def test_idss_attitude_box_at_contact_time():
    """Docking attitude hold: start 10 deg off with 0.5 deg/s rates; by
    any plausible contact time (>200 s) the rotational box holds with
    margin."""
    inertia = np.diag([120.0, 90.0, 60.0])
    ang = math.radians(10.0)
    axis = np.array([0.6, 0.64, 0.48])
    axis = axis / np.linalg.norm(axis)
    q = np.concatenate([[math.cos(ang / 2)], math.sin(ang / 2) * axis])
    w = np.radians([0.5, -0.3, 0.4])
    q_ref = quat.identity()
    dt = 0.1
    for k in range(int(400.0 / dt)):
        tau = ac.quaternion_feedback(q, w, q_ref, np.zeros(3),
                                     kp=1.2, kd=40.0, tau_max=0.5)
        q, w = att.step(q, w, inertia, tau, dt)
        if k * dt == 200.0:
            m_mid = idss.check_attitude(q, q_ref, w, BOX)
            assert m_mid["misalignment"] > 0.0
            assert m_mid["angular_rate"] > 0.0
    m = idss.check_attitude(q, q_ref, w, BOX)
    assert m["misalignment"] > 0.9 * BOX.misalignment_max  # nearly aligned
    assert m["angular_rate"] > 0.5 * BOX.angular_rate_max


def test_box_defaults_match_idss_rev_g():
    assert BOX.closing_min == 0.05
    assert BOX.closing_max == 0.10
    assert BOX.lateral_rate_max == 0.04
    assert BOX.lateral_offset_max == 0.10
    assert abs(math.degrees(BOX.angular_rate_max) - 0.20) < 1e-12
    assert abs(math.degrees(BOX.misalignment_max) - 4.0) < 1e-12


def _docking_case(i, rng):
    x0 = np.array([
        rng.normal(0.0, 0.5), -30.0 + rng.normal(0.0, 2.0),
        rng.normal(0.0, 0.5), 0.0, rng.normal(0.0, 0.005), 0.0,
    ])
    sc = docking_scenario(x0, seed=0)  # seed overwritten by the campaign
    return sc, terminal_controller(), contact_metrics


@pytest.mark.slow
def test_monte_carlo_docking_campaign():
    """20 dispersed docking runs (position/velocity dispersion, sensor
    noise, MIB + execution error): every contact inside the IDSS
    translation box; the campaign table reproduces bit-identically."""
    table = monte_carlo.run_campaign(20, master_seed=42,
                                     make_case=_docking_case)
    assert float(np.min(table["contact"])) == 1.0
    for name in ("closing_above_min", "closing_below_max",
                 "lateral_offset", "lateral_rate"):
        assert float(np.min(table[name])) > 0.0, name
    s = monte_carlo.summary(table, "lateral_offset")
    assert s["min"] > 0.0
    # exact reproducibility of the whole campaign
    table2 = monte_carlo.run_campaign(20, master_seed=42,
                                      make_case=_docking_case)
    assert np.array_equal(table, table2)
    # different master seed gives a different table
    table3 = monte_carlo.run_campaign(20, master_seed=43,
                                      make_case=_docking_case)
    assert not np.array_equal(table, table3)