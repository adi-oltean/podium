"""Magnetic disturbance-torque receipts (#48): the fourth classic
environmental attitude disturbance, tau = m x B, plus a centered-dipole
geomagnetic field model.
"""

import numpy as np

from podium import constants as const
from podium.dynamics import attitude as att
from podium.dynamics.disturbances import DisturbanceModel


def test_torque_equals_dipole_cross_field():
    m = np.array([0.5, -0.2, 0.1])          # residual dipole [A m^2]
    b = np.array([2.0e-5, -1.0e-5, 3.0e-5])  # field [T]
    tau = att.magnetic_torque(m, b)
    assert np.allclose(tau, np.cross(m, b), atol=1e-30)


def test_zero_when_dipole_aligned_with_field():
    b = np.array([1.0e-5, 2.0e-5, -1.0e-5])
    tau = att.magnetic_torque(2.5 * b, b)   # m parallel B
    assert np.allclose(tau, 0.0, atol=1e-30)


def test_dipole_field_magnitude_and_scaling():
    """Equatorial surface field ~3.1e-5 T, and B falls as 1/r^3."""
    # equator: r perpendicular to the -z dipole axis -> |B| = coeff*|m_hat|
    r_eq = np.array([const.R_EARTH, 0.0, 0.0])
    b_eq = att.dipole_field(r_eq)
    mag = float(np.linalg.norm(b_eq))
    assert 2.8e-5 < mag < 3.4e-5
    # double the radius -> field down by 2^3 = 8
    b_2r = att.dipole_field(2.0 * r_eq)
    assert abs(np.linalg.norm(b_2r) - mag / 8.0) < 1e-12


def test_dipole_field_stronger_at_pole():
    """On the dipole axis the field is twice the equatorial value at
    the same radius (the classic 3(m.r)r - m -> 2 factor)."""
    r_eq = np.array([const.R_EARTH, 0.0, 0.0])
    r_pole = np.array([0.0, 0.0, const.R_EARTH])
    ratio = np.linalg.norm(att.dipole_field(r_pole)) / \
        np.linalg.norm(att.dipole_field(r_eq))
    assert abs(ratio - 2.0) < 1e-9


def test_torque_magnitude_reasonable_in_leo():
    """A ~1 A m^2 residual dipole in a ~3e-5 T field gives a
    micro-Nm-class torque, comparable to the other disturbances."""
    r = np.array([const.R_EARTH + 400e3, 0.0, 0.0])
    b = att.dipole_field(r)
    m = np.array([1.0, 0.0, 0.0])           # 1 A m^2, worst-case perp
    tau = att.magnetic_torque(m, b)
    assert 1e-6 < np.linalg.norm(tau) < 1e-4


def test_aggregator_includes_magnetic_term():
    m = DisturbanceModel(inertia=np.diag([180.0, 140.0, 100.0]),
                         n=0.0011,
                         residual_dipole=np.array([0.3, -0.1, 0.2]))
    nadir = np.array([0.0, 0.0, 1.0])       # gg zero at this alignment
    b = np.array([1.0e-5, 2.0e-5, -1.0e-5])
    total = m.torque(nadir, b_field_body=b)
    assert np.allclose(total, att.magnetic_torque(m.residual_dipole, b))
