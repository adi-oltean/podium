"""Truth-model validation: transforms, J2 secular rates, drag decay,
energy conservation, and cross-check against the YA STM."""

import math

import numpy as np
import pytest

from podium import constants as const
from podium.core import ya
from podium.dynamics import nonlinear as nl

MU = const.MU_EARTH


def make_target(a, e, inc=0.9, raan=0.5, argp=1.2, nu=0.3):
    r, v = nl.elements_to_rv(a, e, inc, raan, argp, nu, MU)
    return np.concatenate([r, v])


def test_elements_roundtrip_radius_speed():
    a, e, nu = 7_000_000.0, 0.1, 0.8
    r, v = nl.elements_to_rv(a, e, 0.9, 0.5, 1.2, nu, MU)
    p = a * (1 - e * e)
    assert math.isclose(np.linalg.norm(r), p / (1 + e * math.cos(nu)), rel_tol=1e-12)
    # vis-viva
    assert math.isclose(
        np.linalg.norm(v), math.sqrt(MU * (2 / np.linalg.norm(r) - 1 / a)), rel_tol=1e-12
    )


def test_lvlh_roundtrip():
    cfg = nl.ForceConfig(j2=const.J2_EARTH, drag=nl.DragConfig())
    rv_t = make_target(6_778_137.0, 0.02)
    x = np.array([120.0, -3500.0, 40.0, 0.4, -1.1, 0.02])
    rv_c = nl.lvlh_to_eci(rv_t, x, cfg, bc_target=120.0)
    back = nl.eci_to_lvlh(rv_t, rv_c, cfg, bc_target=120.0)
    assert np.allclose(back, x, rtol=1e-12, atol=1e-9)


@pytest.mark.slow
def test_twobody_energy_conserved():
    rv_t = make_target(6_778_137.0, 0.1)
    x0 = np.zeros(6)
    tof = 2 * math.pi * math.sqrt(6_778_137.0**3 / MU)  # one period
    _, _, rv = nl.propagate_relative(rv_t, x0, tof, dt=2.0)
    def energy(s):
        return 0.5 * np.dot(s[3:6], s[3:6]) - MU / np.linalg.norm(s[0:3])
    e0, e1 = energy(rv[0]), energy(rv[-1])
    assert abs(e1 - e0) / abs(e0) < 1e-10


@pytest.mark.parametrize("e", [0.0, 0.1, 0.3])
def test_matches_ya_stm_twobody(e):
    """End-to-end cross-check of the ECI pipeline against the independent
    YA implementation: agreement to linearization error."""
    a = 8_000_000.0
    nu0 = 1.4
    n = math.sqrt(MU / a**3)
    rv_t = make_target(a, e, nu=nu0)
    x0 = np.array([30.0, -90.0, 15.0, 0.02, 0.01, -0.02])
    tof = 2500.0
    _, x_rel, _ = nl.propagate_relative(rv_t, x0, tof, dt=1.0)
    xf_ya = ya.stm(n, e, nu0, tof) @ x0
    assert np.linalg.norm(x_rel[-1][:3] - xf_ya[:3]) < 0.05
    assert np.linalg.norm(x_rel[-1][3:] - xf_ya[3:]) < 1e-4


def test_lvlh_velocity_is_position_derivative():
    """Central-difference the relative position history and compare with the
    reported relative velocity — validates the frame angular velocity,
    including the out-of-plane (omega_x) term, with J2 + drag active."""
    cfg = nl.ForceConfig(j2=const.J2_EARTH, drag=nl.DragConfig())
    rv_t = make_target(6_778_137.0, 0.01, inc=0.9)
    x0 = np.array([50.0, -2000.0, 500.0, 0.1, 0.3, -0.4])
    dt = 0.25
    times, x_rel, _ = nl.propagate_relative(
        rv_t, x0, 200.0, dt=dt, cfg=cfg, bc_target=150.0, bc_chaser=60.0
    )
    for k in range(10, 790, 97):
        v_num = (x_rel[k + 1, :3] - x_rel[k - 1, :3]) / (2 * dt)
        assert np.linalg.norm(v_num - x_rel[k, 3:]) < 2e-5  # < 0.02 mm/s


@pytest.mark.slow
def test_j2_raan_drift_matches_analytic():
    a, e, inc = 6_778_137.0, 0.001, math.radians(51.6)
    cfg = nl.ForceConfig(j2=const.J2_EARTH)
    rv0 = make_target(a, e, inc=inc, raan=1.0)
    period = 2 * math.pi * math.sqrt(a**3 / MU)
    tof = 10 * period
    _, _, rv = nl.propagate_relative(rv0, np.zeros(6), tof, dt=5.0, cfg=cfg)

    def raan_of(s):
        h = np.cross(s[0:3], s[3:6])
        node = np.cross([0.0, 0.0, 1.0], h)
        return math.atan2(node[1], node[0])

    d_raan = raan_of(rv[-1]) - raan_of(rv[0])
    d_raan = math.atan2(math.sin(d_raan), math.cos(d_raan))
    n = math.sqrt(MU / a**3)
    p = a * (1 - e * e)
    analytic = -1.5 * n * const.J2_EARTH * (const.R_EARTH / p) ** 2 * math.cos(inc) * tof
    # ~ -0.056 rad over 10 orbits; short-period osculating wiggle ~1e-3 rad.
    assert abs(d_raan - analytic) < 0.05 * abs(analytic)


@pytest.mark.slow
def test_drag_decays_semimajor_axis():
    a = 6_778_137.0
    bc = 50.0
    drag = nl.DragConfig(rho0=1e-11, h0=400e3, scale_height=60e3)
    cfg = nl.ForceConfig(drag=drag)
    rv0 = make_target(a, 0.0005, inc=0.9)
    period = 2 * math.pi * math.sqrt(a**3 / MU)
    _, _, rv = nl.propagate_relative(
        rv0, np.zeros(6), 5 * period, dt=5.0, cfg=cfg, bc_target=bc, bc_chaser=bc
    )

    def sma(s):
        r = np.linalg.norm(s[0:3])
        v2 = np.dot(s[3:6], s[3:6])
        return 1.0 / (2.0 / r - v2 / MU)

    da = sma(rv[-1]) - sma(rv[0])
    # Analytic per-orbit decay for circular orbit: 2*pi*rho*a^2/BC (inertial
    # velocity approximation; co-rotation reduces it ~10% at LEO).
    rho = drag.density(a - const.R_EARTH)
    da_analytic = -2 * math.pi * rho * a * a / bc * 5
    assert da < 0
    assert abs(da - da_analytic) / abs(da_analytic) < 0.25


@pytest.mark.slow
def test_differential_drag_along_track():
    """Lower-BC (draggier) chaser falls to a lower, faster orbit and pulls
    ahead along-track (+y drift) — the classic differential-drag signature."""
    cfg = nl.ForceConfig(drag=nl.DragConfig(rho0=1e-11))
    rv0 = make_target(6_778_137.0, 0.0005)
    period = 2 * math.pi * math.sqrt(6_778_137.0**3 / MU)
    _, x_rel, _ = nl.propagate_relative(
        rv0, np.zeros(6), 5 * period, dt=5.0, cfg=cfg, bc_target=200.0, bc_chaser=40.0
    )
    assert x_rel[-1, 1] > 100.0  # chaser drifts ahead by hundreds of meters
