"""Environmental disturbance-torque aggregator.

Composes the three dominant environmental attitude disturbances —
gravity gradient (#44), aerodynamic (#45), and solar radiation pressure
(#46) — into a single body-frame torque, so a control loop can be
exercised against realistic proximity-ops disturbances. Torques
superpose (they are independent physical effects), so the aggregate is
simply their sum; the value of the model is one place that holds the
spacecraft's fixed geometry (inertia, ballistic/area coefficients,
centers of pressure) and turns per-instant direction vectors into a
total torque.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from podium import constants as const
from podium.dynamics import attitude as att

F64 = NDArray[np.float64]


@dataclass
class DisturbanceModel:
    """Fixed spacecraft geometry for the environmental torques.

    inertia : (3,3) body inertia tensor.
    n       : orbital mean motion [rad/s] (gravity gradient).
    aero    : (cd_area, r_cp_body) or None to disable.
    srp     : (area, cr, r_cp_body) or None to disable.
    """

    inertia: F64
    n: float
    aero: tuple[float, F64] | None = None
    srp: tuple[float, float, F64] | None = None
    residual_dipole: F64 | None = None   # body-frame residual dipole [A m^2]

    def torque(
        self,
        nadir_body: F64,
        v_rel_body: F64 | None = None,
        rho: float = 0.0,
        sun_body: F64 | None = None,
        illuminated: bool = True,
        b_field_body: F64 | None = None,
    ) -> F64:
        """Total environmental body torque at one instant. Directions
        are in the body frame; disable any term by leaving its config
        None (or, for aero, rho=0)."""
        tau: F64 = att.gravity_gradient_torque(nadir_body, self.inertia,
                                               self.n)
        if self.aero is not None and v_rel_body is not None and rho > 0.0:
            cd_area, r_cp = self.aero
            tau = tau + att.aerodynamic_torque(v_rel_body, rho, cd_area,
                                               r_cp)
        if self.srp is not None and sun_body is not None:
            area, cr, r_cp = self.srp
            tau = tau + att.srp_torque(sun_body, area, cr, r_cp,
                                       const.SOLAR_PRESSURE, illuminated)
        if self.residual_dipole is not None and b_field_body is not None:
            tau = tau + att.magnetic_torque(self.residual_dipole,
                                            b_field_body)
        return tau
