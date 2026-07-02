"""Fixed-step integrators for the verifiable core.

Only fixed-step schemes are provided here: adaptive step control introduces
data-dependent loop bounds, which the static subset forbids. The simulation
truth models in :mod:`podium.sim` may use SciPy adaptive integrators for
cross-validation, but flight-representative propagation uses these.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]
Deriv = Callable[[float, F64], F64]


def rk4_step(f: Deriv, t: float, x: F64, dt: float) -> F64:
    """Single classical Runge-Kutta 4 step."""
    k1 = f(t, x)
    k2 = f(t + 0.5 * dt, x + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, x + 0.5 * dt * k2)
    k4 = f(t + dt, x + dt * k3)
    return x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def euler_step(f: Deriv, t: float, x: F64, dt: float) -> F64:
    """Single explicit Euler step (for cheap onboard predictors)."""
    return x + dt * f(t, x)
