"""CW kernel: STM vs numerical integration, two-impulse targeting closure."""

import math

import numpy as np
import pytest

from rpod.core import cw, integrators

MU = 3.986004418e14
A_LEO = 6_778_137.0  # ~400 km altitude
N = cw.mean_motion(MU, A_LEO)


def propagate_numeric(x0, n, t, steps=2000):
    u = np.zeros(3)
    f = lambda _t, x: cw.cw_deriv(x, n, u)  # noqa: E731
    dt = t / steps
    x = x0.copy()
    for i in range(steps):
        x = integrators.rk4_step(f, i * dt, x, dt)
    return x


def test_mean_motion_leo():
    # ~92.6 min period at 400 km
    period = 2 * math.pi / N
    assert 5500 < period < 5600


@pytest.mark.parametrize("tof", [60.0, 600.0, 3000.0])
def test_stm_matches_integration(tof):
    x0 = np.array([100.0, -2000.0, 50.0, 0.1, 0.5, -0.05])
    x_stm = cw.stm(N, tof) @ x0
    x_num = propagate_numeric(x0, N, tof)
    assert np.allclose(x_stm, x_num, rtol=1e-9, atol=1e-6)


def test_stm_identity_at_zero():
    assert np.allclose(cw.stm(N, 1e-12), np.eye(6), atol=1e-8)


def test_two_impulse_reaches_target():
    x0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])  # 1 km behind on V-bar
    target = np.zeros(6)  # arrive at origin, at rest
    tof = 1500.0
    dv1, dv2 = cw.two_impulse(x0, target, N, tof)

    xb = x0.copy()
    xb[3:6] += dv1
    x_arr = cw.stm(N, tof) @ xb
    assert np.allclose(x_arr[0:3], target[0:3], atol=1e-6)
    assert np.allclose(x_arr[3:6] + dv2, target[3:6], atol=1e-9)


def test_two_impulse_dv_reasonable():
    # Hop 1 km along V-bar in a quarter period: dv should be cm/s-to-m/s scale.
    x0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])
    tof = math.pi / (2 * N)
    dv1, dv2 = cw.two_impulse(x0, np.zeros(6), N, tof)
    total = np.linalg.norm(dv1) + np.linalg.norm(dv2)
    assert 0.01 < total < 5.0
