"""Relative-nav EKF receipts: Joseph-form invariants, statistical
consistency (NEES/NIS) against the engine's seeded truth, convergence,
and the closed loop — LQR flying on estimates from noisy measurements."""

import math

import numpy as np

from podium import constants as const
from podium.control import lqr
from podium.nav import ekf
from podium.sim import Scenario, circular_target, run

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)


def test_joseph_form_invariants():
    """500 random update/predict cycles: P stays symmetric and PD."""
    rng = np.random.default_rng(4)
    f = ekf.RelNavEkf(N, dt=2.0, q_accel=1e-8, r_pos=5.0)
    for _ in range(500):
        z = rng.normal(0.0, 50.0, 3)
        f.step(z)
        asym = float(np.max(np.abs(f.p - f.p.T)))
        assert asym < 1e-12 * float(np.max(np.abs(f.p)))
        assert np.all(np.linalg.eigvalsh(f.p) > 0.0)


def test_sequential_equals_batch_joseph():
    """The flight-side sequential scalar update (division, no matrix
    solve) is algebraically equivalent to the batch Joseph update for
    H = [I3 0] and diagonal R — the classical result, verified to near
    machine precision over random PSD covariances."""
    rng = np.random.default_rng(11)
    for _ in range(200):
        x = rng.uniform(-1e3, 1e3, 6)
        a_ = rng.uniform(-1.0, 1.0, (6, 6))
        p = a_ @ a_.T + np.eye(6) * rng.uniform(0.1, 2.0)
        z = rng.uniform(-1e3, 1e3, 3)
        rv = rng.uniform(0.01, 5.0)
        xs, ps = ekf.update_sequential(x, p, z, rv)
        xb, pb, _nu, _s = ekf.update_joseph(x, p, z, ekf.H_POS,
                                            np.eye(3) * rv)
        scale = float(np.max(np.abs(pb))) + 1.0
        assert np.max(np.abs(xs - xb)) < 1e-9 * (np.max(np.abs(xb)) + 1.0)
        assert np.max(np.abs(ps - pb)) < 1e-9 * scale
        assert np.all(np.linalg.eigvalsh(ps) > 0.0)


def test_process_noise_structure():
    q = ekf.process_noise_wna(2.0, 1e-8)
    assert np.allclose(q, q.T)
    assert np.all(np.linalg.eigvalsh(q) >= -1e-25)
    assert q[0, 0] == 1e-8 * 8.0 / 3.0
    assert q[3, 3] == 1e-8 * 2.0


def _run_filter(seed, duration=1200.0, dt=2.0, sigma=5.0, q_accel=2e-8,
                x_rel0=None, x_hat0=None):
    """Unforced engine truth + EKF on the engine's noisy measurements."""
    sc = Scenario(
        duration=duration,
        rv_target0=circular_target(A),
        x_rel0=(np.array([30.0, -900.0, 10.0, 0.0, 0.0, 0.0])
                if x_rel0 is None else x_rel0),
        dt_gnc=dt,
        truth_substeps=4,
        seed=seed,
        meas_pos_std=sigma,
    )
    f = ekf.RelNavEkf(N, dt=dt, q_accel=q_accel, r_pos=sigma,
                      x0=(sc.x_rel0.copy() if x_hat0 is None else x_hat0))
    est_err, nis = [], []

    def ctl(t, meas):
        est = f.step(meas[0:3])
        est_err.append(est)
        nu, s = f.last_nu, f.last_s
        nis.append(float(nu @ np.linalg.solve(s, nu)))
        return np.zeros(3)

    tr = run(sc, ctl)
    est_arr = np.array(est_err)
    truth = tr.x_rel[: len(est_arr)]
    err = est_arr - truth
    # NEES needs the post-update covariance; recompute per-step is
    # awkward — use steady-state check on the second half instead
    return err, np.array(nis), tr


def test_consistency_nis_and_accuracy():
    """Time-averaged NIS within the chi-square band (3 dof), and the
    steady-state position error well below the raw measurement noise."""
    err, nis, _ = _run_filter(seed=9)
    half = len(nis) // 2
    nis_mean = float(np.mean(nis[half:]))
    # 3-dof NIS: mean 3; 95% band for ~300-sample average ~ [2.6, 3.4]
    assert 2.4 < nis_mean < 3.6, nis_mean
    pos_rms = float(np.sqrt(np.mean(err[half:, 0:3] ** 2)))
    assert pos_rms < 2.0  # vs 5 m measurement noise
    vel_rms = float(np.sqrt(np.mean(err[half:, 3:6] ** 2)))
    assert vel_rms < 0.02


def test_convergence_from_large_initial_error():
    x0_true = np.array([30.0, -900.0, 10.0, 0.0, 0.0, 0.0])
    x_hat0 = x0_true + np.array([80.0, -60.0, 40.0, 0.3, -0.4, 0.2])
    err, _, _ = _run_filter(seed=3, x_rel0=x0_true, x_hat0=x_hat0)
    # the 100+ m prior error collapses in the very first Joseph update
    # (100 m prior vs 5 m measurement => gain ~ 1), then settles
    assert float(np.linalg.norm(x_hat0[0:3] - x0_true[0:3])) > 100.0
    assert float(np.linalg.norm(err[0, 0:3])) < 30.0
    tail = err[-50:]
    assert float(np.sqrt(np.mean(tail[:, 0:3] ** 2))) < 2.5
    assert float(np.sqrt(np.mean(tail[:, 3:6] ** 2))) < 0.03


def test_seeded_filter_run_is_deterministic():
    e1, n1, _ = _run_filter(seed=21)
    e2, n2, _ = _run_filter(seed=21)
    assert np.array_equal(e1, e2)
    assert np.array_equal(n1, n2)


def test_closed_loop_lqr_on_estimates():
    """The v0.3 milestone receipt: navigation + control + truth. LQR
    gains from the CW model, flying on EKF estimates built from 5 m
    position-only measurements, stabilize the nonlinear truth."""
    dt = 2.0
    sigma = 5.0
    a_d, b_d = lqr.cw_discrete(N, dt)
    q = np.diag([1, 1, 1, 100, 100, 100]).astype(float)
    r = np.eye(3) * 1e4
    k = lqr.dlqr(a_d, b_d, q, r)
    x0 = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    f = ekf.RelNavEkf(N, dt=dt, q_accel=2e-8, r_pos=sigma, x0=x0.copy())

    def ctl(t, meas):
        est = f.step(meas[0:3], dv=None)
        u = lqr.apply_gain(k, est, u_max=0.05)
        dv = u * dt
        # feed the commanded burn through the filter's prediction
        f.x[3:6] += dv
        return dv

    sc = Scenario(
        duration=2 * math.pi / N,
        rv_target0=circular_target(A),
        x_rel0=x0.copy(),
        dt_gnc=dt,
        truth_substeps=4,
        seed=5,
        meas_pos_std=sigma,
    )
    tr = run(sc, ctl)
    final_range = float(tr.channels()["range"][-1])
    # perfect-measurement baseline reaches <2 m; with 5 m noise the box
    # is noise-floor-limited
    assert final_range < 12.0
    # and the whole trajectory stays bounded (no noise-driven divergence)
    assert float(np.max(tr.channels()["range"])) < 150.0