"""The canonical flight-kernel list: everything the emitter translates.

Single source of truth for the golden-vector tests, the EVA gate, and
the release audit-bundle builder. Order is stable (it determines
emission order in kernels.c).
"""

from __future__ import annotations

from podium.core import cw, quat, roe, ya
from podium.nav import ekf

FLIGHT_KERNELS = [
    quat.normalize, quat.multiply, quat.conjugate, quat.rotate,
    quat.deriv, quat.error,
    cw.mean_motion, cw.cw_deriv, cw.stm,
    ya.kepler_eccentric, ya.true_from_eccentric,
    ya.eccentric_from_true, ya.propagate_true_anomaly,
    roe.stm_keplerian, roe.map_roe_to_lvlh, roe.map_lvlh_to_roe,
    roe.control_matrix,
    ekf.predict, ekf.update_sequential, ekf.process_noise_wna,
]
