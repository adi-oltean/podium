"""Terminal docking control with pulsed (minimum-impulse-bit) actuation.

Designed for vehicles whose translation/rotation authority is a fixed
velocity increment per thruster pulse — which is also exactly the control
model of SpaceX's public ISS docking simulator (each button click is one
pulse). The same three functions drive the browser autopilot demo
(``viewer/iss-sim/``); the JS port mirrors these formulas and constants
line for line, with Python as the source of truth.

Axes follow the docking-port frame: axial = distance to port along the
docking axis (positive away from port, decreasing during approach);
lateral = the two perpendicular offsets to null. Angles in degrees here —
the vehicle HUD convention — since these laws act on displayed telemetry.

Flight-side control laws in the restricted style (pure functions, fixed bounds,
contracted inputs); the emitted, verified kernel set is core + nav, not these.
"""

from __future__ import annotations

from podium.verify import Interval, contract

# IDSS-like soft-capture numbers used as defaults; the SpaceX sim's success
# gates are 0.2 m lateral, 0.2 deg per axis, 0.2 m/s closing rate.
CONTACT_RATE = 0.16  # commanded closing rate at contact [m/s] (margin vs 0.2)
FAR_RATE = 3.0  # closing-rate ceiling far from the port [m/s]
APPROACH_TAU = 60.0  # range/rate time constant on final [s]
LATERAL_V_MAX = 0.5  # lateral correction speed ceiling [m/s]
LATERAL_TAU = 25.0  # lateral time constant, far zone [s]
LATERAL_TAU_NEAR = 4.0  # lateral time constant inside the near zone [s]
LATERAL_NEAR_ZONE = 2.0  # near-zone boundary [m]
ANGLE_RATE_MAX = 0.5  # attitude correction rate ceiling [deg/s]
ANGLE_TAU = 10.0  # attitude time constant, far zone [s]
ANGLE_TAU_NEAR = 2.0  # attitude time constant inside the near zone [s]
ANGLE_NEAR_ZONE = 2.0  # near-zone boundary [deg]

# With pulsed actuation a proportional law stalls once |cmd| < dv_pulse/2,
# leaving a standing offset of tau * dv_pulse / 2. The near-zone taus are
# sized so that stall offset sits well inside the 0.2 m / 0.2 deg docking
# gates (e.g. 4 s * 0.05 / 2 = 0.1 m laterally at the sim's fine increment).


@contract(range_m=Interval(0.0, 1e6))
def approach_rate(range_m: float) -> float:
    """Commanded closing rate [m/s] as a function of range to port [m].

    Profile: v = clamp(range/tau, CONTACT_RATE, FAR_RATE) — an exponential-
    decay glideslope in range with a floor so contact happens at a firm but
    capture-safe rate, and a ceiling for the far field.
    """
    v = range_m / APPROACH_TAU
    if v < CONTACT_RATE:
        v = CONTACT_RATE
    elif v > FAR_RATE:
        v = FAR_RATE
    return v


@contract(pos=Interval(-1e6, 1e6), v_max=Interval(1e-3, 100.0), tau=Interval(1e-2, 1e4))
def axis_rate_cmd(pos: float, v_max: float, tau: float) -> float:
    """Saturated proportional command: drive `pos` to zero at rate pos/tau.

    Used for lateral offsets (m -> m/s) and attitude errors (deg -> deg/s)
    alike; the caller picks (v_max, tau) per axis class.
    """
    v = -pos / tau
    if v > v_max:
        v = v_max
    elif v < -v_max:
        v = -v_max
    return v


@contract(pos=Interval(-1e6, 1e6))
def lateral_rate_cmd(pos: float) -> float:
    """Lateral-offset velocity command [m/s], gain-scheduled by zone."""
    tau = LATERAL_TAU if abs(pos) > LATERAL_NEAR_ZONE else LATERAL_TAU_NEAR
    return axis_rate_cmd(pos, LATERAL_V_MAX, tau)


@contract(err_deg=Interval(-360.0, 360.0))
def angle_rate_cmd(err_deg: float) -> float:
    """Attitude-error rate command [deg/s], gain-scheduled by zone."""
    tau = ANGLE_TAU if abs(err_deg) > ANGLE_NEAR_ZONE else ANGLE_TAU_NEAR
    return axis_rate_cmd(err_deg, ANGLE_RATE_MAX, tau)


@contract(dv_pulse=Interval(1e-6, 10.0))
def pulses_needed(v_cmd: float, v_est: float, dv_pulse: float) -> int:
    """Minimum-impulse-bit allocator: signed pulse count to reach v_cmd.

    Rounds toward zero with a half-pulse deadband, so the vehicle never
    chatters around the commanded rate; the residual is always below
    dv_pulse/2. Count is clamped to +/-25 per control cycle (the caller
    re-plans every cycle, so large errors are worked off over a few cycles).
    """
    n = int((v_cmd - v_est) / dv_pulse + (0.5 if v_cmd > v_est else -0.5))
    if n > 25:
        n = 25
    elif n < -25:
        n = -25
    return n
