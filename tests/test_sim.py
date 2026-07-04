"""Sim engine + spec registry receipts: determinism, closed-loop flights
against the nonlinear truth, spec semantics, crossings, viewer export."""

import json
import math

import numpy as np
import pytest

from podium.control import lqr
from podium.core import cw, quat
from podium.dynamics import nonlinear as nl
from podium.guidance.glideslope import glideslope_pulses
from podium.sim import Scenario, circular_target, mean_motion_of, run
from podium.sim import spec as sp

A = 6_778_137.0
TARGET = circular_target(A)
N = mean_motion_of(TARGET)


def no_control(_t, _x):
    return np.zeros(3)


def make_scenario(**kw):
    base = dict(
        duration=600.0,
        rv_target0=TARGET.copy(),
        x_rel0=np.array([50.0, -800.0, 20.0, 0.0, 0.0, 0.0]),
        dt_gnc=2.0,
        truth_substeps=5,
    )
    base.update(kw)
    return Scenario(**base)


def test_bit_identical_replay():
    tr1 = run(make_scenario(seed=7, meas_pos_std=1.0, meas_vel_std=0.01), no_control)
    tr2 = run(make_scenario(seed=7, meas_pos_std=1.0, meas_vel_std=0.01), no_control)
    assert np.array_equal(tr1.x_rel, tr2.x_rel)
    assert np.array_equal(tr1.rv_target, tr2.rv_target)


def test_noise_is_seeded_and_matters():
    captured = {}

    def probe(seed):
        seen = []

        def ctl(t, meas):
            seen.append(meas.copy())
            return np.zeros(3)

        run(make_scenario(duration=20.0, seed=seed, meas_pos_std=2.0), ctl)
        return np.array(seen)

    m1, m2, m3 = probe(1), probe(1), probe(2)
    assert np.array_equal(m1, m2)  # same seed, same noise
    assert not np.array_equal(m1, m3)  # different seed, different noise
    captured  # noqa: B018


def test_unforced_matches_truth_propagation():
    """Engine with no control must reproduce the plain truth propagation."""
    sc = make_scenario(duration=300.0)
    tr = run(sc, no_control)
    _, x_ref, _ = nl.propagate_relative(
        sc.rv_target0, sc.x_rel0, 300.0, sc.dt_gnc / sc.truth_substeps
    )
    assert np.allclose(tr.x_rel[-1], x_ref[-1], rtol=1e-12, atol=1e-9)


@pytest.mark.slow
def test_glideslope_closed_loop_with_specs():
    """The v0.1 canonical scenario flown through the engine, judged by specs."""
    x0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])
    dock = np.array([0.0, -10.0, 0.0])
    duration, pulses = 2400.0, 10
    t_burn, dvs = glideslope_pulses(x0, dock, N, duration, pulses)
    fired = [False] * pulses

    def ctl(t, _meas):
        for i in range(pulses):
            if not fired[i] and t >= t_burn[i] - 1e-9:
                fired[i] = True
                return dvs[i]
        return np.zeros(3)

    specs = (
        sp.always_below("stay_inside_1100m", "range", 1100.0),
        sp.always_below("no_radial_blowup", "x", 60.0),
        sp.eventually_below("arrive_25m", "range", 25.0),
        sp.final_between("terminal_speed", "speed", 0.0, 0.05),
    )
    sc = make_scenario(duration=duration + 60.0, x_rel0=x0, dt_gnc=1.0,
                       truth_substeps=5)
    tr = run(sc, ctl, specs)
    assert all(fired)
    for name, margin in tr.spec_margins.items():
        assert margin > 0.0, f"spec {name} violated: {margin}"
    # arrival: station-keeping ~10 m short of the -10 m hold point
    assert abs(tr.channels()["range"][-1] - 20.0) < 6.0
    assert abs(tr.dv_total() - float(np.sum(np.linalg.norm(dvs, axis=1)))) < 1e-9


@pytest.mark.slow
def test_lqr_stabilizes_truth_through_engine():
    """CW-synthesized gains, impulse-equivalent actuation (dv = u*dt),
    flying the nonlinear truth: 100 m offset driven to <2 m in one orbit."""
    dt = 2.0
    a_d, b_d = lqr.cw_discrete(N, dt)
    q = np.diag([1, 1, 1, 100, 100, 100]).astype(float)
    r = np.eye(3) * 1e4
    k = lqr.dlqr(a_d, b_d, q, r)

    def ctl(_t, meas):
        u = lqr.apply_gain(k, meas, u_max=0.05)
        return u * dt

    sc = make_scenario(
        duration=2 * math.pi / N,
        x_rel0=np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        dt_gnc=dt,
        truth_substeps=4,
    )
    tr = run(sc, ctl)
    assert float(tr.channels()["range"][-1]) < 2.0


def test_spec_semantics_on_synthetic_traces():
    t = np.linspace(0.0, 10.0, 101)
    sig = np.sin(t)  # max ~0.9999, min ~-1.0 over [0,10]
    ch = {"t": t, "s": sig}
    assert sp.evaluate((sp.always_below("m", "s", 2.0),), ch)["m"] == pytest.approx(
        2.0 - np.max(sig)
    )
    assert sp.evaluate((sp.always_above("m", "s", -2.0),), ch)["m"] == pytest.approx(
        np.min(sig) + 2.0
    )
    # windowing: on [0, 1.5] the minimum of sin is at t=0
    assert sp.evaluate(
        (sp.always_above("m", "s", 0.0, t_start=0.5, t_end=1.5),), ch
    )["m"] == pytest.approx(math.sin(0.5))
    ev = sp.evaluate((sp.eventually_below("m", "s", -0.9),), ch)["m"]
    assert ev == pytest.approx(-0.9 - np.min(sig))
    fin = sp.evaluate((sp.final_between("m", "s", -1.0, 0.0),), ch)["m"]
    assert fin == pytest.approx(min(sig[-1] + 1.0, 0.0 - sig[-1]))
    # empty window is loudly violated
    assert sp.evaluate(
        (sp.always_above("m", "s", 0.0, t_start=20.0),), ch
    )["m"] == -math.inf


def test_crossing_times_bracketed():
    sc = make_scenario(duration=400.0,
                       x_rel0=np.array([0.0, -900.0, 0.0, 0.0, 1.0, 0.0]))
    tr = run(sc, no_control)
    crossings = tr.crossing_times("range", 700.0)
    assert len(crossings) >= 1
    t_c = crossings[0]
    ch = tr.channels()
    idx = int(np.searchsorted(tr.times, t_c))
    assert (ch["range"][idx - 1] - 700.0) * (ch["range"][idx] - 700.0) <= 0.0


def test_viewer_json_schema():
    sc = make_scenario(duration=60.0)
    tr = run(sc, lambda t, x: np.array([0.0, 0.01, 0.0]) if t == 0.0 else np.zeros(3))
    js = json.loads(tr.to_viewer_json(name="test", orbit="400 km", n=N))
    assert set(js) == {"meta", "t", "x", "q_le", "burns"}
    assert js["meta"]["dv_total"] == pytest.approx(0.01, abs=1e-6)
    assert len(js["t"]) == len(js["x"]) == len(js["q_le"]) == 31
    assert len(js["x"][0]) == 6
    assert js["burns"][0]["t"] == 0.0


def test_frame_quaternions_physics():
    """The exported LVLH->ECI quaternions: unit norm, hemisphere-
    continuous (slerp-safe), rotation matrix matches lvlh_rotation
    transposed, and consecutive frames rotate by ~n*dt about the orbit
    normal — the physics the viewer's frame blending displays."""
    sc = make_scenario(duration=600.0)
    tr = run(sc, no_control)
    q = tr.frame_quaternions()
    norms = np.linalg.norm(q, axis=1)
    assert np.max(np.abs(norms - 1.0)) < 1e-9
    dots = np.sum(q[1:] * q[:-1], axis=1)
    assert np.all(dots > 0.0)  # continuity
    # conversion correctness on a sample
    for i in (0, 150, 300):
        m = nl.lvlh_rotation(tr.rv_target[i, 0:3], tr.rv_target[i, 3:6]).T
        for col, e in enumerate(np.eye(3)):
            assert np.allclose(quat.rotate(q[i], e), m[:, col], atol=1e-9)
    # frame rate: angle between consecutive quats = n*dt (circular)
    ang = 2.0 * np.arccos(np.clip(dots, -1.0, 1.0))
    dt = tr.times[1] - tr.times[0]
    assert np.allclose(ang, N * dt, rtol=2e-3)


def test_dv_actually_applied():
    """An along-track impulse must change the trajectory consistently with
    the CW response at small times."""
    dv = 0.1
    sc = make_scenario(duration=100.0,
                       x_rel0=np.zeros(6), dt_gnc=1.0, truth_substeps=5)
    tr = run(sc, lambda t, x: np.array([0.0, dv, 0.0]) if t == 0.0 else np.zeros(3))
    x_cw = cw.stm(N, 100.0) @ np.array([0.0, 0.0, 0.0, 0.0, dv, 0.0])
    assert np.allclose(tr.x_rel[-1][:3], x_cw[:3], atol=0.02)
