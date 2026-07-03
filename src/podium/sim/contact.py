"""Probe-drogue docking contact via MuJoCo (truth-model layer).

Import explicitly (`from podium.sim import contact`) — mujoco is an
optional extra, deliberately not pulled in by `podium.sim`.

Geometry (MuJoCo frame; approach axis = +x, mapping to the LVLH
along-track axis is the caller's concern): the drogue funnel opens
toward -x with mouth radius R_MOUTH at x=0, narrowing to the throat at
x=DEPTH; built from N_PLATES convex box plates (MuJoCo convexifies
meshes, so a non-convex funnel must be a union of convex parts), each
oriented by its slant/tangent frame via `xyaxes`. A backstop disc ends
the throat. The chaser is a 500 kg free body with a probe capsule; the
probe TIP is a tracked site. Gravity is zero: the contact event lasts
seconds, over which orbital dynamics contribute micrometers
(0.5*n^2*r*t^2 class terms) — documented approximation.

Capture = probe tip seated in the throat (x_tip > SEAT_X, radial
distance < SEAT_R) continuously for SEAT_DWELL seconds. Outcomes also
report the peak contact force (grows with closing rate — a receipt)
and whether the probe bounced back out.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]

R_MOUTH = 0.30
R_THROAT = 0.08
DEPTH = 0.40  # half-angle atan(0.22/0.40) ~ 29 deg: steep cones (~47 deg
#               in the first cut) reflect too much axial momentum and
#               ballistic probes stall on the wall — measured, fixed
N_PLATES = 8
PROBE_R = 0.035
PROBE_LEN = 0.6
SEAT_X = 0.36
SEAT_R = 0.07
SEAT_DWELL = 0.5
TIP_START_X = -0.15  # tip clearance from the mouth plane at t=0


def _funnel_plates() -> str:
    """The N convex plates forming the funnel, as MJCF geom strings."""
    rows = []
    slant = np.array([DEPTH, R_THROAT - R_MOUTH])
    slant_len = float(np.linalg.norm(slant))
    r_mid = 0.5 * (R_MOUTH + R_THROAT)
    x_mid = 0.5 * DEPTH
    half_w = math.pi * r_mid / N_PLATES  # slight overlap closes gaps
    for i in range(N_PLATES):
        th = 2.0 * math.pi * i / N_PLATES
        u = np.array([0.0, math.cos(th), math.sin(th)])  # radial
        w = np.array([0.0, -math.sin(th), math.cos(th)])  # tangential
        t = (np.array([DEPTH, 0.0, 0.0])
             + (R_THROAT - R_MOUTH) * u) / slant_len  # slant direction
        pos = np.array([x_mid, 0.0, 0.0]) + r_mid * u
        xy = " ".join(f"{v:.6f}" for v in np.concatenate([t, w]))
        p = " ".join(f"{v:.6f}" for v in pos)
        rows.append(
            f'<geom type="box" size="{slant_len / 2:.4f} {half_w:.4f} 0.01"'
            f' pos="{p}" xyaxes="{xy}" friction="0.05 0.001 0.0001"'
            f' rgba="0.6 0.6 0.7 1"/>'
        )
    # throat sleeve: axis-parallel plates so the funnel-to-throat
    # transition is a surface, not a box EDGE (edge contact normals have
    # axial components that stop the probe regardless of friction —
    # observed: the IDSS-corner case stalled right at the junction)
    r_sleeve = R_THROAT + 0.01
    half_ws = math.pi * r_sleeve / N_PLATES
    for i in range(N_PLATES):
        th = 2.0 * math.pi * i / N_PLATES
        u = np.array([0.0, math.cos(th), math.sin(th)])
        w = np.array([0.0, -math.sin(th), math.cos(th)])
        pos = np.array([DEPTH + 0.035, 0.0, 0.0]) + r_sleeve * u
        xy = f"1 0 0 0 {w[1]:.6f} {w[2]:.6f}"  # x along axis, y tangent
        p = " ".join(f"{v:.6f}" for v in pos)
        rows.append(
            f'<geom type="box" size="0.045 {half_ws:.4f} 0.01"'
            f' pos="{p}" xyaxes="{xy}" friction="0.05 0.001 0.0001"'
            f' rgba="0.5 0.5 0.6 1"/>'
        )
    return "\n        ".join(rows)


def _mjcf() -> str:
    return f"""
<mujoco model="podium-dock">
  <option timestep="0.002" gravity="0 0 0"/>
  <worldbody>
    <body name="target">
        {_funnel_plates()}
        <geom name="backstop" type="cylinder" size="0.12 0.01"
              pos="{DEPTH + 0.045:.3f} 0 0" zaxis="1 0 0"
              rgba="0.8 0.5 0.3 1"/>
    </body>
    <body name="chaser" pos="0 0 0">
      <freejoint/>
      <!-- attitude-hold approximation: very large rotational inertia
           stands in for the attitude controller during the seconds-long
           contact event (6-DOF contact attitude is deferred scope) -->
      <inertial pos="0 0 0" mass="500" diaginertia="1e5 1e5 1e5"/>
      <geom name="bus" type="box" size="0.5 0.5 0.5" pos="-0.5 0 0"
            rgba="0.4 0.7 0.5 1"/>
      <geom name="probe" type="capsule" size="{PROBE_R}"
            fromto="0 0 0 {PROBE_LEN} 0 0"
            friction="0.05 0.001 0.0001" rgba="0.9 0.9 0.4 1"/>
      <site name="tip" pos="{PROBE_LEN} 0 0" size="0.005"/>
    </body>
  </worldbody>
</mujoco>
"""


@dataclass
class ContactOutcome:
    captured: bool
    bounced: bool
    max_tip_x: float
    seat_time: float  # continuous seated dwell achieved [s]
    peak_force: float  # peak contact normal force [N]


def simulate_contact(
    closing_rate: float,
    lateral_offset: float = 0.0,
    lateral_rate: float = 0.0,
    offset_dir: tuple[float, float] = (1.0, 0.0),  # (y, z) unit-ish
    duration: float = 10.0,
    thrust: float = 0.0,
) -> ContactOutcome:
    """One deterministic contact case (positive closing_rate moves the
    chaser toward the drogue along +x). `thrust` [N] is a sustained
    axial docking thrust during the event — probe-drogue practice
    (Soyuz-class 'hot dock'): the chaser keeps pushing through capture
    so wall contacts can't stall the probe."""
    import mujoco

    model = mujoco.MjModel.from_xml_string(_mjcf())
    data = mujoco.MjData(model)
    d = np.array([offset_dir[0], offset_dir[1]])
    d = d / (np.linalg.norm(d) or 1.0)
    # free joint qpos: [x y z qw qx qy qz]
    data.qpos[0] = TIP_START_X - PROBE_LEN
    data.qpos[1] = lateral_offset * d[0]
    data.qpos[2] = lateral_offset * d[1]
    data.qpos[3] = 1.0
    data.qvel[0] = closing_rate
    data.qvel[1] = lateral_rate * d[0]
    data.qvel[2] = lateral_rate * d[1]
    tip_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, "tip")

    steps = int(duration / model.opt.timestep)
    seat_time = 0.0
    best_dwell = 0.0
    max_tip_x = -math.inf
    peak_force = 0.0
    entered = False
    force_buf = np.zeros(6)
    chaser_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "chaser")
    for _ in range(steps):
        if thrust != 0.0:
            data.xfrc_applied[chaser_id, 0] = thrust
        mujoco.mj_step(model, data)
        tip = data.site_xpos[tip_id]
        max_tip_x = max(max_tip_x, float(tip[0]))
        if float(tip[0]) > 0.0:
            entered = True
        r_tip = math.hypot(float(tip[1]), float(tip[2]))
        if float(tip[0]) > SEAT_X and r_tip < SEAT_R:
            seat_time += model.opt.timestep
            best_dwell = max(best_dwell, seat_time)
        else:
            seat_time = 0.0
        for c in range(data.ncon):
            mujoco.mj_contactForce(model, data, c, force_buf)
            peak_force = max(peak_force, float(force_buf[0]))
        if best_dwell >= SEAT_DWELL:
            break
    tip = data.site_xpos[tip_id]
    bounced = entered and float(tip[0]) < -0.05 and best_dwell < SEAT_DWELL
    return ContactOutcome(
        captured=best_dwell >= SEAT_DWELL,
        bounced=bounced,
        max_tip_x=max_tip_x,
        seat_time=best_dwell,
        peak_force=peak_force,
    )


def capture_envelope(
    offsets: F64,
    closing_rates: F64,
    lateral_rate: float = 0.0,
) -> list[dict]:
    """Grid sweep: capture success over (lateral offset, closing rate)."""
    out = []
    for off in offsets:
        for cr in closing_rates:
            o = simulate_contact(float(cr), float(off), lateral_rate)
            out.append({
                "offset": float(off), "closing": float(cr),
                "captured": o.captured, "peak_force": o.peak_force,
            })
    return out
