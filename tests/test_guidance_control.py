"""Glideslope guidance and LQR closed-loop sanity checks."""

import numpy as np

from rpod.control import lqr
from rpod.core import cw
from rpod.guidance.glideslope import glideslope_pulses

MU = 3.986004418e14
N = cw.mean_motion(MU, 6_778_137.0)


def test_glideslope_arrives_at_dock():
    x0 = np.array([0.0, -500.0, 0.0, 0.0, 0.0, 0.0])
    dock = np.array([0.0, -10.0, 0.0])  # hold point 10 m out on V-bar
    duration, pulses = 1800.0, 8
    times, dvs = glideslope_pulses(x0, dock, N, duration, pulses)

    # Replay the impulses through the STM.
    x = x0.copy()
    for i in range(pulses - 1):
        x[3:6] += dvs[i]
        x = cw.stm(N, times[i + 1] - times[i]) @ x
    x[3:6] += dvs[pulses - 1]

    assert np.linalg.norm(x[0:3] - dock) < 5.0  # ~1% of initial range
    assert np.linalg.norm(x[3:6]) < 1e-9  # terminal pulse nulls velocity


def test_glideslope_range_monotone():
    x0 = np.array([50.0, -800.0, 20.0, 0.0, 0.0, 0.0])
    dock = np.zeros(3)
    times, dvs = glideslope_pulses(x0, dock, N, 2400.0, 12)
    ranges = []
    x = x0.copy()
    for i in range(11):
        ranges.append(np.linalg.norm(x[0:3] - dock))
        x[3:6] += dvs[i]
        x = cw.stm(N, times[i + 1] - times[i]) @ x
    assert all(b < a * 1.05 for a, b in zip(ranges, ranges[1:]))


def test_lqr_stabilizes_cw():
    dt = 1.0
    a, b = lqr.cw_discrete(N, dt)
    q = np.diag([1, 1, 1, 100, 100, 100]).astype(float)
    r = np.eye(3) * 1e4
    k = lqr.dlqr(a, b, q, r)

    # Closed-loop spectral radius strictly inside unit circle.
    eig = np.linalg.eigvals(a - b @ k)
    assert np.max(np.abs(eig)) < 1.0

    # Drive a 100 m offset to <1 m within an orbit, honoring saturation.
    x = np.array([100.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    for _ in range(int(2 * np.pi / N / dt)):
        u = lqr.apply_gain(k, x, u_max=0.05)
        x = a @ x + b @ u
    assert np.linalg.norm(x[0:3]) < 1.0


def test_apply_gain_saturates():
    k = np.eye(3, 6) * 1e6
    u = lqr.apply_gain(k, np.ones(6), u_max=0.05)
    assert np.all(np.abs(u) <= 0.05 + 1e-15)
