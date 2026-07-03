"""Passive-safety metrics on relative orbital elements.

The flight-proven passive-safety idea (D'Amico & Montenbruck 2006; TAFF,
PRISMA heritage): the radial/cross-track (RN-plane) separation of a
near-circular relative orbit is independent of the along-track offset,
so if the RN-plane trajectory never enters the keep-out radius, the
formation is safe under arbitrary along-track drift — the dominant
uncertainty direction. Alignment of the relative e- and i-vectors keeps
the RN-plane ellipse away from the origin.

Functions here are pure and bounded (fixed 360-point scan), usable as
guidance constraints, sim monitors, and test oracles.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]

_N_SCAN = 360  # fixed scan resolution over one orbit of argument of latitude
_TWO_PI_OVER_N = 2.0 * math.pi / _N_SCAN


def ei_separation_angle(roe: F64) -> float:
    """Angle [rad, 0..pi] between the relative e- and i-vectors.

    Passive safety wants (anti-)parallel vectors: angle near 0 or pi.
    Returns pi/2 (worst case) if either vector is negligibly small, since
    a vanishing vector gives no phasing protection by itself.
    """
    de = math.hypot(roe[2], roe[3])
    di = math.hypot(roe[4], roe[5])
    if de < 1e-12 or di < 1e-12:
        return 0.5 * math.pi
    dot = (roe[2] * roe[4] + roe[3] * roe[5]) / (de * di)
    return math.acos(min(1.0, max(-1.0, dot)))


def min_rn_separation(roe: F64, a: float) -> float:
    """Minimum radial/cross-track-plane separation [m] over one orbit.

    Uses the near-circular first-order map: x/a = da - dex cos u - dey
    sin u, z/a = dix sin u - diy cos u; the along-track coordinate is
    excluded by construction, so this is a lower bound on 3-D separation
    for ANY along-track offset — the e/i-vector separation concept
    operationalized. Fixed 360-point scan (0.02% worst-case bound error
    for pure harmonics at 1-degree resolution).
    """
    min_sep2 = math.inf
    for k in range(_N_SCAN):
        u = _TWO_PI_OVER_N * k
        x = roe[0] - roe[2] * math.cos(u) - roe[3] * math.sin(u)
        z = roe[4] * math.sin(u) - roe[5] * math.cos(u)
        sep2 = x * x + z * z
        min_sep2 = min(min_sep2, sep2)
    return a * math.sqrt(min_sep2)


def min_rn_separation_analytic(roe: F64, a: float) -> float:
    """Exact minimum RN-plane separation [m] (no scan).

    The RN-plane trajectory is r(u) = c + P [cos u, sin u]^T with
    c = (da, 0), P = [[-dex, -dey], [-diy, dix]] (units of a). Stationary
    points of ||r(u)||^2 satisfy Re(alpha z + beta z^2) = 0 on |z| = 1
    (z = e^{iu}), i.e. the quartic beta z^4 + alpha z^3 + conj(alpha) z
    + conj(beta) = 0 with alpha = b1 + i b0, beta = M01 - 0.5 (M11 -
    M00) i, where M = P^T P and b = P^T c. We take the minimum of
    ||r(u)|| over the unit-circle roots — exact up to root-finding
    precision, validated against dense scans in the tests. The 360-point
    scan remains the bounded static-subset-friendly variant.
    """
    c = np.array([float(roe[0]), 0.0])
    p = np.array([[-float(roe[2]), -float(roe[3])],
                  [-float(roe[5]), float(roe[4])]])
    m = p.T @ p
    b = p.T @ c
    alpha = complex(b[1], b[0])
    beta = complex(m[0, 1], -0.5 * (m[1, 1] - m[0, 0]))
    coeffs = np.array([beta, alpha, 0.0, np.conj(alpha), np.conj(beta)])
    if np.max(np.abs(coeffs)) < 1e-300:
        # f is constant: circle centered on the origin (or zero geometry)
        return a * float(np.linalg.norm(c))
    # strip leading (numerical) zeros so np.roots sees the true degree
    nz = int(np.argmax(np.abs(coeffs) > 1e-18 * float(np.max(np.abs(coeffs)))))
    roots = np.roots(coeffs[nz:])
    best = math.inf
    for z in roots:
        if abs(abs(z) - 1.0) < 1e-6:
            u = float(np.angle(z))
            w = np.array([math.cos(u), math.sin(u)])
            best = min(best, float(np.linalg.norm(c + p @ w)))
    if math.isinf(best):  # numerically degenerate: fall back to the scan
        return min_rn_separation(roe, a)
    return a * best


def rn_margin(roe: F64, a: float, keep_out_radius: float) -> float:
    """Passive-safety margin [m]: min RN-plane separation minus KOZ radius.

    Positive => the free-drift relative orbit cannot enter the keep-out
    sphere regardless of along-track drift (to first order, unperturbed
    e/i geometry). Use as a guidance constraint or a sim monitor.
    """
    return min_rn_separation(roe, a) - keep_out_radius
