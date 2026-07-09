"""Quaternion-feedback attitude control (flight-side control law).

The classical quaternion-feedback regulator: torque = -Kp * e - Kd * w
with e the shortest-way small-angle error vector from the core kernel
(2 * sign(dq_w) * vec(dq)) and per-axis saturation. For small angles
about a principal axis with inertia J this is a second-order system with
wn = sqrt(kp/J), zeta = kd / (2 sqrt(kp J)) — the tests hold the step
response to that prediction, so the gains mean what the textbook says.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from podium.core import quat
from podium.verify import Interval, contract

F64 = NDArray[np.float64]


@contract(kp=Interval(1e-6, 1e3), kd=Interval(1e-6, 1e4), tau_max=Interval(1e-6, 1e3))
def quaternion_feedback(
    q: F64, w: F64, q_ref: F64, w_ref: F64, kp: float, kd: float, tau_max: float
) -> F64:
    """Body-frame torque command [N m], per-axis saturated."""
    e = quat.error(q, q_ref)
    out = np.empty(3)
    for i in range(3):
        u = -kp * e[i] - kd * (w[i] - w_ref[i])
        out[i] = min(tau_max, max(-tau_max, u))
    return out
