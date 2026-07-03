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


def rn_margin(roe: F64, a: float, keep_out_radius: float) -> float:
    """Passive-safety margin [m]: min RN-plane separation minus KOZ radius.

    Positive => the free-drift relative orbit cannot enter the keep-out
    sphere regardless of along-track drift (to first order, unperturbed
    e/i geometry). Use as a guidance constraint or a sim monitor.
    """
    return min_rn_separation(roe, a) - keep_out_radius
