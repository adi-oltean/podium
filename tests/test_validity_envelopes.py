"""CW / Yamanaka-Ankersen validity envelopes against the nonlinear truth.

These tests are the quantified statement of where the linearized guidance
models may be used (roadmap v0.1). Bounds are physics-derived, not tuned:
YA error is pure linearization error, O(sep^2 / a); CW adds an O(e * sep)
eccentricity error on top.
"""

import math

import numpy as np
import pytest

from podium import constants as const
from podium.core import cw, ya
from podium.dynamics import nonlinear as nl

MU = const.MU_EARTH
A = 7_000_000.0
N = math.sqrt(MU / A**3)
PERIOD = 2 * math.pi / N
NU0 = 0.9

# Direction chosen so all separations excite radial, along-track, and
# cross-track dynamics; velocity scaled to natural-motion magnitude (n*sep).
UNIT_STATE = np.array([0.3, -0.9, 0.3, 0.0, 0.0, 0.0])
UNIT_VEL = np.array([0.0, 0.0, 0.0, 0.5, -0.5, 0.7])


def propagate_all(sep, e, tof):
    x0 = sep * UNIT_STATE + sep * N * UNIT_VEL
    rv_t = np.concatenate(nl.elements_to_rv(A, e, 0.9, 0.5, 1.2, NU0, MU))
    _, x_rel, _ = nl.propagate_relative(rv_t, x0, tof, dt=2.0)
    truth = x_rel[-1]
    x_ya = ya.stm(N, e, NU0, tof) @ x0
    x_cw = cw.stm(N, tof) @ x0
    return truth, x_ya, x_cw


# Envelope constants C(e): YA position error after one orbit is below
# C * sep^2 / a. Measured on this grid: ~5-25 for e <= 0.05; ~130 at e = 0.2
# (the nonlinear term is amplified near perigee, where r shrinks and the
# relative dynamics speed up). Bounds carry ~50% margin over measurement.
YA_ENVELOPE_C = {0.0: 40.0, 0.05: 40.0, 0.2: 200.0}


@pytest.mark.parametrize("e", [0.0, 0.05, 0.2])
@pytest.mark.slow
def test_ya_error_quadratic_envelope(e):
    """Two claims, both physics-derived: (1) YA error scales quadratically
    with separation (pure linearization error — a model bug would scale
    linearly); (2) the absolute error stays below C(e) * sep^2 / a."""
    seps = [100.0, 1_000.0, 10_000.0]
    errs = []
    for sep in seps:
        truth, x_ya, _ = propagate_all(sep, e, PERIOD)
        err = np.linalg.norm(x_ya[:3] - truth[:3])
        errs.append(err)
        assert err < YA_ENVELOPE_C[e] * sep * sep / A
    # Quadratic scaling: each 10x in separation gives ~100x in error.
    assert 50.0 < errs[1] / errs[0] < 200.0
    assert 50.0 < errs[2] / errs[1] < 200.0


@pytest.mark.parametrize("sep", [100.0, 1_000.0, 10_000.0])
@pytest.mark.slow
def test_cw_equals_ya_at_zero_ecc(sep):
    truth, x_ya, x_cw = propagate_all(sep, 0.0, PERIOD)
    assert np.allclose(x_cw, x_ya, rtol=1e-9, atol=1e-9)
    assert np.linalg.norm(x_cw[:3] - truth[:3]) < 40.0 * sep * sep / A


@pytest.mark.parametrize("e", [0.05, 0.2])
@pytest.mark.slow
def test_cw_degrades_with_eccentricity(e):
    """CW error at eccentricity e is dominated by the O(e*sep) model error
    and dwarfs YA's linearization error — the reason YA exists. Documented
    ratio: >30x at e=0.05, >100x at e=0.2 for 1 km separation."""
    sep = 1_000.0
    truth, x_ya, x_cw = propagate_all(sep, e, PERIOD)
    err_ya = np.linalg.norm(x_ya[:3] - truth[:3])
    err_cw = np.linalg.norm(x_cw[:3] - truth[:3])
    assert err_cw > 30.0 * err_ya
    # CW eccentricity error is O(e * sep) per orbit, with an O(10) dynamics
    # amplification factor; it must at least exceed the bare e*sep scale.
    assert err_cw > e * sep
