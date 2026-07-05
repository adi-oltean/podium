"""Solar-radiation-pressure torque receipts (#46): the third classic
environmental attitude disturbance (dominant at GEO), validated against
tau = r_cp x F_srp with F_srp = -P Cr A s_hat, plus the eclipse cutoff.
"""

import numpy as np

from podium import constants as const
from podium.dynamics import attitude as att


def test_torque_equals_cp_cross_srp_force():
    """tau = r_cp x (-P Cr A s_hat) exactly."""
    s = np.array([0.6, -0.8, 0.0])
    area, cr = 12.0, 1.8
    r_cp = np.array([0.2, 0.1, -0.4])
    tau = att.srp_torque(s, area, cr, r_cp)
    f = -const.SOLAR_PRESSURE * cr * area * (s / np.linalg.norm(s))
    assert np.allclose(tau, np.cross(r_cp, f), atol=1e-20)


def test_force_pushes_away_from_sun():
    """The SRP force is anti-sunward; its magnitude is P Cr A."""
    s = np.array([1.0, 0.0, 0.0])           # sun along +x
    area, cr = 10.0, 2.0
    # torque from a cp offset perpendicular to the sun line reveals the
    # force direction: r_cp = +z, F along -x -> tau = z x (-x) = -y ... check
    tau = att.srp_torque(s, area, cr, np.array([0.0, 0.0, 1.0]))
    f_expected = np.array([-const.SOLAR_PRESSURE * cr * area, 0.0, 0.0])
    assert np.allclose(tau, np.cross([0, 0, 1], f_expected), atol=1e-20)
    assert np.linalg.norm(f_expected) == const.SOLAR_PRESSURE * cr * area


def test_zero_when_cp_on_sun_line():
    s = np.array([0.0, 0.0, 1.0])
    tau = att.srp_torque(s, 10.0, 1.5, np.array([0.0, 0.0, -0.7]))
    assert np.allclose(tau, 0.0, atol=1e-20)


def test_eclipse_zeroes_the_torque():
    """In Earth's shadow SRP vanishes."""
    s = np.array([0.3, 0.4, 0.5])
    r_cp = np.array([0.5, -0.2, 0.1])
    lit = att.srp_torque(s, 20.0, 1.9, r_cp, illuminated=True)
    dark = att.srp_torque(s, 20.0, 1.9, r_cp, illuminated=False)
    assert np.linalg.norm(lit) > 0.0
    assert np.array_equal(dark, np.zeros(3))


def test_reflectivity_scales_force_linearly():
    """A perfect reflector (Cr=2) pushes exactly twice as hard as a
    perfect absorber (Cr=1)."""
    s = np.array([1.0, 0.2, 0.0])
    r_cp = np.array([0.0, 0.0, 0.5])
    t1 = att.srp_torque(s, 8.0, 1.0, r_cp)
    t2 = att.srp_torque(s, 8.0, 2.0, r_cp)
    assert np.allclose(t2, 2.0 * t1, atol=1e-20)
