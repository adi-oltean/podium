"""Yamanaka-Ankersen state transition for Tschauner-Hempel relative dynamics.

Closed-form STM for relative motion about an **elliptic** target orbit
(0 <= e < 1), exact within the linearization. Reference: Yamanaka &
Ankersen, "New State Transition Matrix for Relative Motion on an Arbitrary
Elliptical Orbit," JGCD 25(1):60-66, 2002.

Frame and state as in :mod:`podium.core.cw`: target-centered LVLH, x radial
(zenith), y along-track, z cross-track; state [rx, ry, rz, vx, vy, vz], SI.

Method: physical states map to Tschauner-Hempel "tilde" variables
(positions scaled by rho = 1 + e cos(theta), derivatives taken with respect
to true anomaly). In tilde space the in-plane equations

    x~'' = 2 y~' + (3/rho) x~,     y~'' = -2 x~',     z~'' = -z~

admit the closed-form fundamental solutions below (verified by direct
substitution); the out-of-plane motion is a pure harmonic in theta. The
STM over a time step is assembled as Lambda^-1(theta1) Phi~ Lambda(theta0).

Static-subset compliant: fixed shapes, bounded loops (Kepler solve runs a
fixed 20 Newton iterations; the denominator 1 - e cos E >= 1 - e stays
positive under the eccentricity contract).
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

from podium.verify import Interval, contract

F64 = NDArray[np.float64]

_TWO_PI = 2.0 * math.pi
_KEPLER_ITERS = 20
# In-plane tilde components (x~, y~, x~', y~') within the 6-state ordering.
_IP = (0, 1, 3, 4)


def eccentric_from_true(theta: float, e: float) -> float:
    """Eccentric anomaly E from true anomaly theta (same half-plane)."""
    return 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(0.5 * theta),
        math.sqrt(1.0 + e) * math.cos(0.5 * theta),
    )


def true_from_eccentric(ecc_anom: float, e: float) -> float:
    """True anomaly theta from eccentric anomaly E (same half-plane)."""
    return 2.0 * math.atan2(
        math.sqrt(1.0 + e) * math.sin(0.5 * ecc_anom),
        math.sqrt(1.0 - e) * math.cos(0.5 * ecc_anom),
    )


def kepler_eccentric(mean_anom: float, e: float) -> float:
    """Solve Kepler's equation M = E - e sin E for E.

    Fixed 20 Newton iterations from E0 = M + e sin M — far past quadratic
    convergence for the contracted range e <= 0.9.
    """
    ecc = mean_anom + e * math.sin(mean_anom)
    for _ in range(_KEPLER_ITERS):
        f = ecc - e * math.sin(ecc) - mean_anom
        fp = 1.0 - e * math.cos(ecc)
        ecc = ecc - f / fp
    return ecc


def propagate_true_anomaly(n: float, e: float, theta0: float, dt: float) -> float:
    """True anomaly (in [0, 2*pi)) after time dt from true anomaly theta0."""
    ecc0 = eccentric_from_true(theta0, e)
    mean0 = ecc0 - e * math.sin(ecc0)
    mean1 = math.fmod(mean0 + n * dt, _TWO_PI)
    if mean1 < 0.0:
        mean1 += _TWO_PI
    ecc1 = kepler_eccentric(mean1, e)
    theta1 = true_from_eccentric(ecc1, e)
    return theta1


def _phi_inplane(e: float, theta: float, j: float) -> F64:
    """In-plane fundamental matrix on (x~, y~, x~', y~') at anomaly theta.

    j is the Yamanaka-Ankersen integral J = k^2 (t - t0); the columns are
    the four fundamental solutions (two periodic, the secular J-solution,
    and the constant along-track offset).
    """
    sin_t = math.sin(theta)
    cos_t = math.cos(theta)
    sin_2t = math.sin(2.0 * theta)
    cos_2t = math.cos(2.0 * theta)
    rho = 1.0 + e * cos_t
    s = rho * sin_t
    c = rho * cos_t
    sp = cos_t + e * cos_2t
    cp = -(sin_t + e * sin_2t)

    m = np.zeros((4, 4))
    m[0, 0] = s
    m[0, 1] = c
    m[0, 2] = 2.0 - 3.0 * e * s * j
    m[1, 0] = c * (1.0 + 1.0 / rho)
    m[1, 1] = -s * (1.0 + 1.0 / rho)
    m[1, 2] = -3.0 * rho * rho * j
    m[1, 3] = 1.0
    m[2, 0] = sp
    m[2, 1] = cp
    m[2, 2] = -3.0 * e * (sp * j + sin_t / rho)
    m[3, 0] = -2.0 * s
    m[3, 1] = e - 2.0 * c
    m[3, 2] = -3.0 * (1.0 - 2.0 * e * s * j)
    return m


def _lambda_to_tilde(e: float, k2: float, theta: float) -> F64:
    """Physical LVLH state -> tilde state at anomaly theta."""
    rho = 1.0 + e * math.cos(theta)
    es = e * math.sin(theta)
    lam = np.zeros((6, 6))
    for i in range(3):
        lam[i, i] = rho
        lam[i + 3, i] = -es
        lam[i + 3, i + 3] = 1.0 / (k2 * rho)
    return lam


def _lambda_from_tilde(e: float, k2: float, theta: float) -> F64:
    """Tilde state -> physical LVLH state at anomaly theta."""
    rho = 1.0 + e * math.cos(theta)
    es = e * math.sin(theta)
    lam = np.zeros((6, 6))
    for i in range(3):
        lam[i, i] = 1.0 / rho
        lam[i + 3, i] = k2 * es
        lam[i + 3, i + 3] = k2 * rho
    return lam


@contract(n=Interval(1e-5, 1e-2), e=Interval(0.0, 0.9))
def stm(n: float, e: float, theta0: float, dt: float) -> F64:
    """Yamanaka-Ankersen state transition matrix Phi(dt), 6x6.

    x(t0 + dt) = Phi @ x(t0) for unforced Tschauner-Hempel motion about a
    target orbit with mean motion n and eccentricity e, starting at target
    true anomaly theta0. Reduces exactly to the CW STM for e = 0.

    dt may be negative (backward transition).
    """
    k2 = n / math.pow(1.0 - e * e, 1.5)
    theta1 = propagate_true_anomaly(n, e, theta0, dt)
    j = k2 * dt

    # In-plane tilde transition: Phi_f(theta1, J) @ Phi_f(theta0, 0)^-1.
    p0 = _phi_inplane(e, theta0, 0.0)
    p1 = _phi_inplane(e, theta1, j)
    pin = np.linalg.solve(p0.T, p1.T).T

    # Out-of-plane: harmonic in theta (periodic, so mod-2*pi anomalies work).
    dth = theta1 - theta0
    cz = math.cos(dth)
    sz = math.sin(dth)

    tilde = np.zeros((6, 6))
    for a in range(4):
        for b in range(4):
            tilde[_IP[a], _IP[b]] = pin[a, b]
    tilde[2, 2] = cz
    tilde[2, 5] = sz
    tilde[5, 2] = -sz
    tilde[5, 5] = cz

    return _lambda_from_tilde(e, k2, theta1) @ tilde @ _lambda_to_tilde(e, k2, theta0)
