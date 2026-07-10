"""IDSS contact-condition box (IDD Rev G soft-capture class) + checkers.

The acceptance geometry: the chaser approaches the dock point along the
approach axis (unit vector pointing from the port toward the incoming
chaser). Closing rate is the velocity component INTO the port; lateral
quantities are perpendicular to the axis. Margins > 0 mean inside the
box. Attitude conditions are checked separately (translation and
rotation are decoupled until the 6-DOF engine lands in v0.4 — the
docking test couples them by evaluating both at the contact time).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from podium.core import quat

F64 = NDArray[np.float64]


@dataclass(frozen=True)
class IdssBox:
    """IDSS IDD Rev G contact-condition defaults (SI)."""

    closing_min: float = 0.05  # [m/s]
    closing_max: float = 0.10  # [m/s]
    lateral_rate_max: float = 0.04  # [m/s]
    lateral_offset_max: float = 0.10  # [m]
    angular_rate_max: float = math.radians(0.20)  # [rad/s]
    misalignment_max: float = math.radians(4.0)  # [rad]


def check_translation(
    x_rel: F64,
    dock_point: F64,
    approach_axis: F64,
    box: IdssBox = IdssBox(),
) -> dict[str, float]:
    """Margins of the translational contact conditions (>0 = inside)."""
    axis = np.asarray(approach_axis, dtype=np.float64)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm == 0.0:
        raise ValueError("approach_axis must be nonzero")
    axis = axis / axis_norm
    rel_pos = x_rel[0:3] - np.asarray(dock_point, dtype=np.float64)
    closing = -float(x_rel[3:6] @ axis)  # positive = moving into the port
    lat_pos = rel_pos - (rel_pos @ axis) * axis
    lat_vel = x_rel[3:6] - (x_rel[3:6] @ axis) * axis
    return {
        "closing_above_min": closing - box.closing_min,
        "closing_below_max": box.closing_max - closing,
        "lateral_offset": box.lateral_offset_max - float(np.linalg.norm(lat_pos)),
        "lateral_rate": box.lateral_rate_max - float(np.linalg.norm(lat_vel)),
    }


def check_attitude(
    q: F64, q_ref: F64, w: F64, box: IdssBox = IdssBox()
) -> dict[str, float]:
    """Margins of the rotational contact conditions (>0 = inside)."""
    dq = quat.multiply(quat.conjugate(q_ref), q)
    misalign = 2.0 * math.acos(min(1.0, abs(float(dq[0]))))
    return {
        "misalignment": box.misalignment_max - misalign,
        "angular_rate": box.angular_rate_max - float(np.linalg.norm(w)),
    }
