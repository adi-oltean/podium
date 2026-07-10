"""Sensor + actuator receipts: budgets match statistics, Jacobian pinned
by FD, camera-EKF closed loop through the engine, MIB/execution-error
actuator effects, determinism."""

import math

import numpy as np
import pytest

from podium import constants as const
from podium.control import lqr
from podium.guidance.glideslope import glideslope_pulses
from podium.nav import ekf, sensors
from podium.sim import Scenario, circular_target, run

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)


def test_camera_jacobian_matches_fd():
    x = np.array([120.0, -800.0, 90.0, 0.1, -0.4, 0.05])
    jac = sensors.camera_jacobian(x)
    h = 1e-5
    fd = np.zeros((3, 6))
    for k in range(6):
        dp = np.zeros(6)
        dp[k] = h
        fd[:, k] = (sensors.camera_h(x + dp) - sensors.camera_h(x - dp)) / (2 * h)
    assert np.allclose(jac, fd, rtol=1e-6, atol=1e-10)


def test_camera_h_rejects_zero_range():
    """Coincident target/chaser (zero range) is non-physical: camera_h must
    raise rather than emit a silent NaN elevation into the EKF."""
    with pytest.raises(ValueError, match="range"):
        sensors.camera_h(np.zeros(6))
    # normal geometry still returns finite [az, el, range]
    z = sensors.camera_h(np.array([10.0, -20.0, 5.0, 0.0, 0.0, 0.0]))
    assert np.all(np.isfinite(z))


def test_camera_jacobian_rejects_degenerate_geometry():
    """Zero range and a line of sight along the cross-track axis both make
    a denominator vanish; the Jacobian must raise, not return NaN."""
    with pytest.raises(ValueError, match="range"):
        sensors.camera_jacobian(np.zeros(6))
    # sx = sy = 0 (LOS on the +/- cross-track axis): azimuth undefined
    with pytest.raises(ValueError, match="cross-track"):
        sensors.camera_jacobian(np.array([0.0, 0.0, 50.0, 0.0, 0.0, 0.0]))
    # normal geometry still yields a finite Jacobian
    jac = sensors.camera_jacobian(np.array([10.0, -20.0, 5.0, 0.0, 0.0, 0.0]))
    assert np.all(np.isfinite(jac))


def test_sensor_budgets_match_statistics():
    rng = np.random.default_rng(0)
    x = np.array([0.0, -500.0, 0.0, 0.0, 0.0, 0.0])
    cam = sensors.DockingCamera()
    zs = np.array([cam.measure(x, rng) for _ in range(4000)])
    truth = sensors.camera_h(x)
    assert abs(float(np.std(zs[:, 0])) - cam.bearing_std) < 0.1 * cam.bearing_std
    assert abs(float(np.mean(zs[:, 2])) - truth[2]) < 0.5
    assert (
        abs(float(np.std(zs[:, 2])) - cam.range_std_frac * truth[2])
        < 0.1 * cam.range_std_frac * truth[2]
    )
    gnss = sensors.RelGnss(pos_std=2.0, vel_std=0.02, bias_pos_std=1.0)
    bias = gnss.start(rng)
    zs_g = np.array([gnss.measure(x, rng, bias) for _ in range(4000)])
    err = zs_g[:, 0:3] - x[0:3]
    assert np.allclose(np.mean(err, axis=0), bias, atol=0.15)
    assert abs(float(np.std(err[:, 0])) - 2.0) < 0.2
    # visibility gating
    assert not cam.visible(np.array([0.0, -3000.0, 0.0, 0, 0, 0]))
    assert sensors.Lidar().range_max == 1_000.0


def test_camera_ekf_through_engine():
    """EKF on bearing+range measurements only (no direct position),
    unforced drift through the engine: converges and tracks."""
    dt = 2.0
    cam = sensors.DockingCamera(range_max=5_000.0)
    meas_rng = np.random.default_rng(77)
    x0 = np.array([50.0, -800.0, 20.0, 0.0, 0.0, 0.0])
    f = ekf.RelNavEkf(N, dt=dt, q_accel=2e-8, r_pos=5.0,
                      x0=x0 + np.array([30.0, -40.0, 10.0, 0, 0, 0]))
    errs = []

    def ctl(t, true_rel):
        z = cam.measure(true_rel, meas_rng)
        rho = float(np.linalg.norm(f.x[0:3])) or 1.0
        f.x, f.p, _, _ = ekf.update_joseph_nonlinear(
            f.x, f.p, z, sensors.camera_h, sensors.camera_jacobian,
            cam.noise_cov(rho), angle_rows=(0, 1),
        )
        errs.append(f.x - true_rel)
        f.x, f.p = ekf.predict(f.x, f.p, f.phi, f.q)
        return np.zeros(3)

    sc = Scenario(duration=1200.0, rv_target0=circular_target(A),
                  x_rel0=x0.copy(), dt_gnc=dt, truth_substeps=4)
    run(sc, ctl)
    err = np.array(errs)
    tail = err[-100:]
    pos_rms = float(np.sqrt(np.mean(tail[:, 0:3] ** 2)))
    # bearing 0.1 deg at ~800 m ~ 1.4 m cross-range; range 1% ~ 8 m
    assert pos_rms < 6.0
    assert float(np.linalg.norm(err[0, 0:3])) > 20.0  # started far off


def test_actuator_mib_and_cap_visible_in_burn_log():
    sc = Scenario(
        duration=40.0, rv_target0=circular_target(A),
        x_rel0=np.zeros(6), dt_gnc=2.0, truth_substeps=2,
        dv_quantum=0.01, dv_max_tick=0.05,
    )

    def ctl(t, _x):
        if t == 0.0:
            return np.array([0.0, 0.123, 0.0])  # gets capped then quantized
        if t == 10.0:
            return np.array([0.004, 0.0, 0.0])  # below half a click: dropped
        if t == 20.0:
            return np.array([0.0, 0.0, 0.017])  # rounds to 0.02
        return np.zeros(3)

    tr = run(sc, ctl)
    assert len(tr.burns) == 2  # the sub-MIB burn vanished
    t0, dv0 = tr.burns[0]
    assert t0 == 0.0
    assert abs(float(dv0[1]) - 0.05) < 1e-12  # capped to 0.05, on-grid
    t1, dv1 = tr.burns[1]
    assert t1 == 20.0
    assert abs(float(dv1[2]) - 0.02) < 1e-12


def test_actuator_execution_error_openloop_vs_feedback():
    """The lesson the actuator model teaches: 2% execution error on a
    1.4 m/s open-loop insertion burn drifts ~3*dv*t ~ hundreds of
    meters (physics, not a bug) — while a feedback law flying through
    the same noisy actuator still converges. Both runs bit-identically
    replayable."""
    x0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])
    dock = np.array([0.0, -10.0, 0.0])
    t_burn, dvs = glideslope_pulses(x0, dock, N, 2400.0, 10)

    def make_openloop():
        fired = [False] * 10

        def ctl(t, _m):
            for i in range(10):
                if not fired[i] and t >= t_burn[i] - 1e-9:
                    fired[i] = True
                    return dvs[i]
            return np.zeros(3)
        return ctl

    def scenario():
        return Scenario(
            duration=2460.0, rv_target0=circular_target(A), x_rel0=x0.copy(),
            dt_gnc=2.0, truth_substeps=4, seed=13,
            dv_quantum=0.001, dv_exec_std_frac=0.02,
        )

    tr1 = run(scenario(), make_openloop())
    tr2 = run(scenario(), make_openloop())
    assert np.array_equal(tr1.x_rel, tr2.x_rel)  # seeded exec error replays
    # actual burns differ from the (quantized) commands and are logged
    act0 = tr1.burns[0][1]
    assert not np.allclose(act0, np.round(dvs[0] / 0.001) * 0.001)
    openloop_miss = abs(float(tr1.channels()["range"][-1]) - 20.0)
    assert openloop_miss > 100.0  # open-loop fragility, quantified

    # feedback through the same imperfect actuator absorbs it
    dt = 2.0
    a_d, b_d = lqr.cw_discrete(N, dt)
    k = lqr.dlqr(a_d, b_d, np.diag([1, 1, 1, 100, 100, 100]).astype(float),
                 np.eye(3) * 1e4)

    def fb(t, meas):
        return lqr.apply_gain(k, meas, u_max=0.05) * dt

    sc = Scenario(
        duration=2 * math.pi / N, rv_target0=circular_target(A),
        x_rel0=np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        dt_gnc=dt, truth_substeps=4, seed=13,
        dv_quantum=0.001, dv_exec_std_frac=0.02,
    )
    tr = run(sc, fb)
    assert float(tr.channels()["range"][-1]) < 5.0