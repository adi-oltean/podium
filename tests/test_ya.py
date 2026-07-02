"""Yamanaka-Ankersen STM: validated against nonlinear elliptic relative
dynamics (truth), the CW limit at e=0, and structural STM properties."""

import math

import numpy as np
import pytest

from rpod.core import cw, ya

MU = 3.986004418e14


def truth_deriv(state, mu, h):
    """Nonlinear LVLH relative dynamics about an elliptic reference.

    state = [x, y, z, vx, vy, vz, R, Rdot, theta]; the reference orbit
    (R, Rdot, theta) is integrated self-consistently alongside.
    """
    x, y, z, vx, vy, vz, r_ref, rd_ref, _ = state
    thd = h / (r_ref * r_ref)
    thdd = -2.0 * rd_ref * thd / r_ref
    rc = math.sqrt((r_ref + x) ** 2 + y * y + z * z)
    g = mu / rc**3
    ax = 2.0 * thd * vy + thdd * y + thd * thd * x - g * (r_ref + x) + mu / r_ref**2
    ay = -2.0 * thd * vx - thdd * x + thd * thd * y - g * y
    az = -g * z
    rdd = -mu / r_ref**2 + h * h / r_ref**3
    return np.array([vx, vy, vz, ax, ay, az, rd_ref, rdd, thd])


def integrate_truth(x0, a, e, theta0, tof, steps):
    p = a * (1.0 - e * e)
    h = math.sqrt(MU * p)
    rho0 = 1.0 + e * math.cos(theta0)
    r0 = p / rho0
    rd0 = math.sqrt(MU / p) * e * math.sin(theta0)
    s = np.concatenate([x0, [r0, rd0, theta0]])
    dt = tof / steps
    for _ in range(steps):
        k1 = truth_deriv(s, MU, h)
        k2 = truth_deriv(s + 0.5 * dt * k1, MU, h)
        k3 = truth_deriv(s + 0.5 * dt * k2, MU, h)
        k4 = truth_deriv(s + dt * k3, MU, h)
        s = s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return s[:6]


@pytest.mark.parametrize(
    ("a", "e", "theta0"),
    [
        (6_778_137.0, 0.0, 0.0),
        (7_200_000.0, 0.05, 0.7),
        (8_000_000.0, 0.15, 2.1),
        (10_500_000.0, 0.35, 4.5),
    ],
)
def test_matches_nonlinear_truth(a, e, theta0):
    n = math.sqrt(MU / a**3)
    x0 = np.array([20.0, -80.0, 10.0, 0.01, 0.02, -0.01])
    tof = 2500.0
    xf_truth = integrate_truth(x0, a, e, theta0, tof, steps=5000)
    xf = ya.stm(n, e, theta0, tof) @ x0
    # Residual vs truth is pure linearization error, ~ |x|^2 / R: at ~600 m
    # final separation that is centimeters. An implementation bug would show
    # up meters-to-kilometers scale.
    assert np.linalg.norm(xf[:3] - xf_truth[:3]) < 0.05
    assert np.linalg.norm(xf[3:] - xf_truth[3:]) < 1e-4


def test_truth_residual_scales_quadratically():
    """Halving separation must quarter the truth residual — distinguishes
    linearization error (quadratic) from an STM bug (linear)."""
    a, e, theta0, tof = 8_000_000.0, 0.15, 2.1, 2500.0
    n = math.sqrt(MU / a**3)
    x1 = np.array([200.0, -800.0, 100.0, 0.1, 0.2, -0.1])
    errs = []
    for scale in (1.0, 0.5, 0.25):
        x0 = scale * x1
        xf_truth = integrate_truth(x0, a, e, theta0, tof, steps=5000)
        xf = ya.stm(n, e, theta0, tof) @ x0
        errs.append(np.linalg.norm(xf[:3] - xf_truth[:3]))
    assert errs[0] / errs[1] > 3.5  # ~4 expected for quadratic
    assert errs[1] / errs[2] > 3.5


def test_reduces_to_cw_at_zero_ecc():
    n = cw.mean_motion(MU, 6_778_137.0)
    for theta0 in (0.0, 1.3, 4.0):
        for tof in (100.0, 1700.0, 6000.0):
            assert np.allclose(
                ya.stm(n, 0.0, theta0, tof), cw.stm(n, tof), rtol=1e-9, atol=1e-12
            )


def test_identity_at_zero_dt():
    assert np.allclose(ya.stm(1e-3, 0.2, 1.0, 0.0), np.eye(6), atol=1e-12)


def test_composition():
    n, e, theta0 = 9e-4, 0.25, 0.8
    t1, t2 = 900.0, 1400.0
    theta_mid = ya.propagate_true_anomaly(n, e, theta0, t1)
    lhs = ya.stm(n, e, theta0, t1 + t2)
    rhs = ya.stm(n, e, theta_mid, t2) @ ya.stm(n, e, theta0, t1)
    assert np.allclose(lhs, rhs, rtol=1e-8, atol=1e-10)


def test_backward_inverts_forward():
    n, e, theta0, tof = 1.1e-3, 0.4, 2.0, 1800.0
    theta1 = ya.propagate_true_anomaly(n, e, theta0, tof)
    fwd = ya.stm(n, e, theta0, tof)
    back = ya.stm(n, e, theta1, -tof)
    assert np.allclose(back @ fwd, np.eye(6), rtol=1e-8, atol=1e-9)


def test_det_is_one():
    for e in (0.0, 0.1, 0.5, 0.9):
        phi = ya.stm(1e-3, e, 0.6, 3000.0)
        assert math.isclose(np.linalg.det(phi), 1.0, rel_tol=1e-9)


def test_kepler_roundtrip():
    for e in (0.0, 0.3, 0.9):
        for theta in np.linspace(-3.0, 3.0, 13):
            ecc = ya.eccentric_from_true(theta, e)
            mean = ecc - e * math.sin(ecc)
            ecc2 = ya.kepler_eccentric(mean, e)
            theta2 = ya.true_from_eccentric(ecc2, e)
            assert math.isclose(
                math.atan2(math.sin(theta2 - theta), math.cos(theta2 - theta)),
                0.0,
                abs_tol=1e-10,
            )


def test_anomaly_propagation_full_period():
    n, e, theta0 = 1e-3, 0.3, 1.1
    period = 2.0 * math.pi / n
    theta1 = ya.propagate_true_anomaly(n, e, theta0, period)
    assert math.isclose(theta1, theta0, abs_tol=1e-8)
