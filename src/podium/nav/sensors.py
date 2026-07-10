"""Sensor truth models with error budgets, plus flight-side measurement
functions/Jacobians for the EKF.

Truth models are classes with `measure(t, x_rel, rng)` — the engine calls
them with ITS seeded Generator, so determinism is inherited, never
re-seeded here. Budgets are 1-sigma white noise per axis unless stated.
The flight-side pieces (camera_h / camera_jacobian) are pure fixed-shape
functions in the static-subset style.

Camera/lidar geometry: measurement = [azimuth, elevation, range] of the
TARGET as seen from the chaser: line of sight s = -r (r = chaser
position in target LVLH), az = atan2(s_y, s_x), el = asin(s_z/|s|),
rho = |s| = |r|.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]


def camera_h(x: F64) -> F64:
    """Flight-side measurement prediction [az, el, range] from the
    relative LVLH state (position part used)."""
    s = -x[0:3]
    rho = float(np.linalg.norm(s))
    if rho <= 0.0:
        raise ValueError(
            "camera_h: degenerate geometry, relative range is zero "
            "(target coincident with chaser)"
        )
    az = math.atan2(float(s[1]), float(s[0]))
    el = math.asin(float(s[2]) / rho)
    return np.array([az, el, rho])


def camera_jacobian(x: F64) -> F64:
    """d[az, el, range]/d(state), (3, 6); velocity columns are zero.
    Pinned by finite differences in the tests."""
    s = -x[0:3]
    sx, sy, sz = float(s[0]), float(s[1]), float(s[2])
    rho2 = sx * sx + sy * sy + sz * sz
    rho = math.sqrt(rho2)
    rxy2 = sx * sx + sy * sy
    rxy = math.sqrt(rxy2)
    if rho <= 0.0 or rxy <= 0.0:
        raise ValueError(
            "camera_jacobian: degenerate geometry, relative range is zero or "
            "the line of sight lies on the cross-track axis (bearing undefined)"
        )
    h = np.zeros((3, 6))
    # d(az)/ds then chain ds/dr = -I
    h[0, 0] = -(-sy / rxy2)
    h[0, 1] = -(sx / rxy2)
    # d(el)/ds = d asin(sz/rho): (rho^2 e_z - sz s) / (rho^2 * rxy)
    h[1, 0] = -(-sz * sx / (rho2 * rxy))
    h[1, 1] = -(-sz * sy / (rho2 * rxy))
    h[1, 2] = -(rxy / rho2)
    # d(range)/ds = s/rho
    h[2, 0] = -(sx / rho)
    h[2, 1] = -(sy / rho)
    h[2, 2] = -(sz / rho)
    return h


@dataclass(frozen=True)
class RelGnss:
    """Relative GNSS: position (+ optional velocity) with white noise and
    an optional constant per-run bias drawn once from the engine rng."""

    pos_std: float = 1.0  # [m]
    vel_std: float = 0.01  # [m/s]
    bias_pos_std: float = 0.0  # constant-bias budget [m]

    def start(self, rng: np.random.Generator) -> F64:
        """Draw the per-run constant bias (call once at scenario start)."""
        if self.bias_pos_std > 0.0:
            return rng.normal(0.0, self.bias_pos_std, 3)
        return np.zeros(3)

    def measure(self, x_rel: F64, rng: np.random.Generator, bias: F64) -> F64:
        z = np.empty(6)
        z[0:3] = x_rel[0:3] + bias + rng.normal(0.0, self.pos_std, 3)
        z[3:6] = x_rel[3:6] + rng.normal(0.0, self.vel_std, 3)
        return z


@dataclass(frozen=True)
class DockingCamera:
    """Bearing + range of the target; valid inside range_max only."""

    bearing_std: float = math.radians(0.1)  # [rad] az and el
    range_std_frac: float = 0.01  # range noise as a fraction of range
    range_max: float = 2_000.0  # [m]

    def visible(self, x_rel: F64) -> bool:
        return float(np.linalg.norm(x_rel[0:3])) <= self.range_max

    def measure(self, x_rel: F64, rng: np.random.Generator) -> F64:
        z = camera_h(x_rel)
        z[0] += rng.normal(0.0, self.bearing_std)
        z[1] += rng.normal(0.0, self.bearing_std)
        z[2] *= 1.0 + rng.normal(0.0, self.range_std_frac)
        return z

    def noise_cov(self, rho: float) -> F64:
        return np.diag([
            self.bearing_std**2,
            self.bearing_std**2,
            (self.range_std_frac * rho) ** 2,
        ])


@dataclass(frozen=True)
class Lidar(DockingCamera):
    """Same geometry as the camera, tighter budget, shorter reach."""

    bearing_std: float = math.radians(0.02)
    range_std_frac: float = 0.0005
    range_max: float = 1_000.0
