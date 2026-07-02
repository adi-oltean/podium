"""Glideslope approach guidance (Hablani et al., JGCD 2002).

The chaser closes range along a straight line to the docking port with range
rate proportional to range: rdot = -k * r + rdot_T, giving an exponential
decay to a soft arrival rate. Impulses are computed against the CW state
transition over each guidance interval.

Static-subset compliant: the number of pulses is a compile-time parameter.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

from rpod.core import cw
from rpod.verify import Interval, contract

F64 = NDArray[np.float64]


@contract(
    n=Interval(1e-5, 1e-2),          # LEO/MEO mean motions [rad/s]
    duration=Interval(10.0, 86_400.0),
    num_pulses=Interval(2, 64),
)
def glideslope_pulses(
    x0: F64,
    dock: F64,
    n: float,
    duration: float,
    num_pulses: int,
) -> tuple[F64, F64]:
    """Compute impulsive glideslope guidance toward position `dock` (3,).

    Parameters
    ----------
    x0 : (6,) chaser LVLH state at start.
    dock : (3,) target-fixed docking point in LVLH [m].
    n : target mean motion [rad/s].
    duration : total approach time [s].
    num_pulses : number of impulses (uniform spacing).

    Returns
    -------
    times : (num_pulses,) burn times from start [s].
    dvs : (num_pulses, 3) LVLH velocity increments [m/s].

    The final pulse nulls relative velocity at arrival (station-keep at dock).
    """
    dt = duration / float(num_pulses - 1)
    rho0 = x0[0:3] - dock
    r0 = math.sqrt(rho0[0] ** 2 + rho0[1] ** 2 + rho0[2] ** 2)
    # Decay rate chosen so range contracts to ~1% over the approach.
    k = -math.log(0.01) / duration

    times = np.empty(num_pulses)
    dvs = np.zeros((num_pulses, 3))
    x = x0.copy()
    for i in range(num_pulses - 1):
        t_i = dt * float(i)
        # Waypoint on the glideslope line at the next pulse time.
        r_next = r0 * math.exp(-k * (t_i + dt))
        u_hat = rho0 / r0 if r0 > 0.0 else np.zeros(3)
        wp = np.empty(6)
        wp[0:3] = dock + u_hat * r_next
        wp[3:6] = 0.0
        dv1, _ = cw.two_impulse(x, wp, n, dt)
        times[i] = t_i
        dvs[i, :] = dv1
        # Propagate the post-burn state to the next pulse.
        xb = x.copy()
        xb[3:6] = xb[3:6] + dv1
        x = cw.stm(n, dt) @ xb
    # Terminal pulse: null residual velocity.
    times[num_pulses - 1] = duration
    dvs[num_pulses - 1, :] = -x[3:6]
    return times, dvs
