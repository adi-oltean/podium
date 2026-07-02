"""Clohessy-Wiltshire (Hill) relative-motion dynamics and state transition.

Frame: target-centered LVLH with x radial (away from Earth), y along-track,
z cross-track (orbit normal completes the right-handed triad). All quantities
SI (m, m/s, rad/s). Valid for near-circular target orbits and separations
small relative to the orbit radius; use :mod:`podium.dynamics.th` for eccentric
targets and :mod:`podium.dynamics.nonlinear` as truth.

Static-subset compliant: fixed shapes, no allocation beyond small constant
temporaries, no data-dependent branching.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]


def mean_motion(mu: float, a: float) -> float:
    """Mean motion n = sqrt(mu/a^3) of the target orbit [rad/s]."""
    return math.sqrt(mu / (a * a * a))


def cw_deriv(x: F64, n: float, u: F64) -> F64:
    """CW state derivative. x = [rx, ry, rz, vx, vy, vz], u = accel [m/s^2] in LVLH.

    xdotdot = 3n^2 x + 2n ydot + ux
    ydotdot = -2n xdot + uy
    zdotdot = -n^2 z + uz
    """
    out = np.empty(6)
    out[0] = x[3]
    out[1] = x[4]
    out[2] = x[5]
    out[3] = 3.0 * n * n * x[0] + 2.0 * n * x[4] + u[0]
    out[4] = -2.0 * n * x[3] + u[1]
    out[5] = -n * n * x[2] + u[2]
    return out


def stm(n: float, t: float) -> F64:
    """Closed-form CW state-transition matrix Phi(t), 6x6.

    x(t) = Phi(t) @ x(0) for unforced motion.
    """
    s = math.sin(n * t)
    c = math.cos(n * t)
    phi = np.zeros((6, 6))
    # position rows
    phi[0, 0] = 4.0 - 3.0 * c
    phi[0, 3] = s / n
    phi[0, 4] = 2.0 * (1.0 - c) / n
    phi[1, 0] = 6.0 * (s - n * t)
    phi[1, 1] = 1.0
    phi[1, 3] = 2.0 * (c - 1.0) / n
    phi[1, 4] = (4.0 * s - 3.0 * n * t) / n
    phi[2, 2] = c
    phi[2, 5] = s / n
    # velocity rows
    phi[3, 0] = 3.0 * n * s
    phi[3, 3] = c
    phi[3, 4] = 2.0 * s
    phi[4, 0] = 6.0 * n * (c - 1.0)
    phi[4, 3] = -2.0 * s
    phi[4, 4] = 4.0 * c - 3.0
    phi[5, 2] = -n * s
    phi[5, 5] = c
    return phi


def two_impulse(x0: F64, target: F64, n: float, tof: float) -> tuple[F64, F64]:
    """Two-impulse CW targeting: velocity increments to reach `target` position
    (with `target` velocity) in time-of-flight `tof`.

    Returns (dv1, dv2), each shape (3,), applied at t=0 and t=tof in LVLH.
    Raises no exceptions for singular transfer times in the static subset;
    callers must respect the contract n*tof not a multiple of pi (declared in
    podium.verify contracts) — the matrix inverse is then well-conditioned.
    """
    phi = stm(n, tof)
    prr = phi[0:3, 0:3]
    prv = phi[0:3, 3:6]
    pvr = phi[3:6, 0:3]
    pvv = phi[3:6, 3:6]

    # Required initial velocity so that r(tof) = r_target.
    v_req = np.linalg.solve(prv, target[0:3] - prr @ x0[0:3])
    dv1 = v_req - x0[3:6]
    # Velocity arriving at tof, and the trim burn to match target velocity.
    v_arr = pvr @ x0[0:3] + pvv @ v_req
    dv2 = target[3:6] - v_arr
    return dv1, dv2
