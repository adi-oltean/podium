"""Stochastic atmosphere receipts: OU statistics, determinism, exact
baseline reduction, and trajectory-level effect within the envelope."""

import math

import numpy as np
import pytest

from podium import constants as const
from podium.dynamics import nonlinear as nl

MU = const.MU_EARTH


def test_ou_is_seeded_and_deterministic():
    p1 = nl.DensityPerturbation(seed=42, duration=3600.0)
    p2 = nl.DensityPerturbation(seed=42, duration=3600.0)
    p3 = nl.DensityPerturbation(seed=43, duration=3600.0)
    ts = np.linspace(0.0, 3600.0, 61)
    f1 = [p1.factor(t) for t in ts]
    assert f1 == [p2.factor(t) for t in ts]
    assert f1 != [p3.factor(t) for t in ts]


def test_ou_statistics_match_parameters():
    """Stationary std and lag-tau autocorrelation of the exact OU
    discretization (long series, fixed seed => deterministic test)."""
    tau, sigma = 21_600.0, 0.35
    dt = tau / 20.0
    pert = nl.DensityPerturbation(seed=7, duration=400 * tau, dt=dt,
                                  sigma_log=sigma, tau=tau)
    p = pert._p
    assert abs(float(np.std(p)) - sigma) < 0.2 * sigma
    lag = int(round(tau / dt))
    c = float(np.corrcoef(p[:-lag], p[lag:])[0, 1])
    assert abs(c - math.exp(-1.0)) < 0.1


def test_calibration_is_in_storm_band():
    """+2-sigma excursion of the default calibration sits in the
    documented +50-125% storm band."""
    pert = nl.DensityPerturbation(seed=0, duration=1.0)
    two_sigma_factor = math.exp(2.0 * pert.sigma_log)
    assert 1.5 <= two_sigma_factor <= 2.25


def test_zero_sigma_is_exactly_baseline():
    drag_base = nl.DragConfig(rho0=1e-11)
    pert = nl.DensityPerturbation(seed=5, duration=10_000.0, sigma_log=0.0)
    drag_zero = nl.DragConfig(rho0=1e-11, perturbation=pert)
    for h in (350e3, 400e3, 500e3):
        for t in (0.0, 1234.5, 9999.0):
            assert drag_zero.density(h, t) == drag_base.density(h)


def test_clamps_beyond_grid():
    pert = nl.DensityPerturbation(seed=3, duration=600.0, dt=60.0)
    assert pert.factor(1e9) == pert.factor(1e8)  # holds last value
    assert pert.factor(-5.0) == math.exp(float(pert._p[0]))


@pytest.mark.slow
def test_perturbed_decay_within_envelope_and_replays():
    """Trajectory-level receipt: perturbed drag decay is bit-identically
    reproducible under the same seed, differs from baseline, and stays
    within the factor bounds implied by the perturbation's own extremes."""
    a = 6_778_137.0
    bc = 50.0
    period = 2 * math.pi * math.sqrt(a**3 / MU)
    tof = 3 * period
    rv0 = np.concatenate(
        nl.elements_to_rv(a, 0.0005, 0.9, 0.5, 1.2, 0.3, MU)
    )

    def decay(drag):
        cfg = nl.ForceConfig(drag=drag)
        _, _, rv = nl.propagate_relative(
            rv0, np.zeros(6), tof, dt=5.0, cfg=cfg, bc_target=bc, bc_chaser=bc
        )
        def sma(s):
            r = np.linalg.norm(s[0:3])
            return 1.0 / (2.0 / r - np.dot(s[3:6], s[3:6]) / MU)
        return sma(rv[-1]) - sma(rv[0]), rv

    base = nl.DragConfig(rho0=1e-11)
    pert = nl.DensityPerturbation(seed=11, duration=tof)
    stoch = nl.DragConfig(rho0=1e-11, perturbation=pert)

    da_base, _ = decay(base)
    da_p1, rv_p1 = decay(stoch)
    # replay: fresh perturbation object, same seed -> identical trajectory
    stoch2 = nl.DragConfig(
        rho0=1e-11, perturbation=nl.DensityPerturbation(seed=11, duration=tof)
    )
    da_p2, rv_p2 = decay(stoch2)
    assert np.array_equal(rv_p1, rv_p2)
    assert da_p1 != da_base
    # decay ratio bounded by the realized factor extremes
    ts = np.arange(0.0, tof, 60.0)
    factors = np.array([pert.factor(t) for t in ts])
    ratio = da_p1 / da_base
    assert float(factors.min()) - 1e-9 <= ratio <= float(factors.max()) + 1e-9