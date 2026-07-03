"""Deterministic fixed-step simulation engine.

One master GNC clock; truth dynamics (dual-ECI nonlinear model,
podium.dynamics.nonlinear) integrate between ticks with a fixed number of
RK4 substeps. Flight blocks are pure step functions called through the
same interface they will have after C translation; the v0 actuation
interface is impulsive: the controller returns a Δv (LVLH, m/s) applied
instantaneously at the tick — which covers impulsive guidance plans,
discrete LQR (Δv = u·dt), and pulsed docking control.

Determinism is non-negotiable: the single `numpy` Generator seeded from
the scenario is the only randomness (measurement noise); identical
scenario + seed give bit-identical traces, which the test suite enforces.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from podium import constants as const
from podium.core import integrators
from podium.dynamics import nonlinear as nl
from podium.sim import spec as spec_mod

F64 = NDArray[np.float64]

# (t, measured relative LVLH state) -> impulsive dv in LVLH [m/s], shape (3,)
Controller = Callable[[float, F64], F64]


@dataclass
class Scenario:
    """Everything that defines a run; two equal scenarios replay identically."""

    duration: float
    rv_target0: F64  # target ECI state at t=0, shape (6,)
    x_rel0: F64  # chaser relative LVLH state at t=0, shape (6,)
    dt_gnc: float = 1.0  # master clock period [s]
    truth_substeps: int = 10  # RK4 steps per tick
    cfg: nl.ForceConfig = field(default_factory=nl.ForceConfig)
    bc_target: float = 100.0
    bc_chaser: float = 100.0
    seed: int = 0
    meas_pos_std: float = 0.0  # per-axis measurement noise [m]
    meas_vel_std: float = 0.0  # per-axis measurement noise [m/s]
    # actuator imperfections (hardware truth, applied to commanded burns;
    # the burn log records what was ACTUALLY applied):
    dv_quantum: float = 0.0  # minimum impulse bit [m/s]; 0 = ideal
    dv_max_tick: float = math.inf  # per-tick impulse magnitude cap [m/s]
    dv_exec_std_frac: float = 0.0  # per-axis proportional execution error


@dataclass
class Trace:
    """Recorded run: states on the master clock plus burn log and margins."""

    times: F64  # (N+1,)
    x_rel: F64  # (N+1, 6) relative LVLH
    rv_target: F64  # (N+1, 6) target ECI
    burns: list[tuple[float, F64]]
    spec_margins: dict[str, float]

    def channels(self) -> dict[str, F64]:
        """Named scalar channels for specs, monitors, and plotting."""
        x = self.x_rel
        rng = np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2 + x[:, 2] ** 2)
        speed = np.sqrt(x[:, 3] ** 2 + x[:, 4] ** 2 + x[:, 5] ** 2)
        safe = np.where(rng > 1e-9, rng, 1.0)
        rate = (x[:, 0] * x[:, 3] + x[:, 1] * x[:, 4] + x[:, 2] * x[:, 5]) / safe
        return {
            "t": self.times,
            "x": x[:, 0], "y": x[:, 1], "z": x[:, 2],
            "vx": x[:, 3], "vy": x[:, 4], "vz": x[:, 5],
            "range": rng, "range_rate": rate, "speed": speed,
        }

    def dv_total(self) -> float:
        return float(sum(float(np.linalg.norm(dv)) for _, dv in self.burns))

    def crossing_times(self, channel: str, threshold: float) -> list[float]:
        """Times where the channel crosses the threshold (linear interp)."""
        ch = self.channels()
        s = ch[channel] - threshold
        t = self.times
        out: list[float] = []
        for i in range(len(s) - 1):
            if s[i] == 0.0:
                out.append(float(t[i]))
            elif s[i] * s[i + 1] < 0.0:
                frac = s[i] / (s[i] - s[i + 1])
                out.append(float(t[i] + frac * (t[i + 1] - t[i])))
        return out

    def to_viewer_json(self, name: str = "podium scenario", orbit: str = "",
                       dock: tuple[float, float, float] = (0.0, 0.0, 0.0),
                       n: float = 0.0) -> str:
        """Serialize to the schema the live viewer loads."""
        data = {
            "meta": {
                "name": name,
                "orbit": orbit,
                "n": n,
                "dt": float(self.times[1] - self.times[0]) if len(self.times) > 1 else 0.0,
                "dv_total": round(self.dv_total(), 4),
                "dock": list(dock),
                "truth": "nonlinear ECI two-craft, LVLH differencing",
            },
            "t": [round(float(t), 3) for t in self.times],
            "x": [[round(float(v), 4) for v in row] for row in self.x_rel],
            "burns": [
                {"t": round(float(t), 3), "dv": [round(float(v), 5) for v in dv]}
                for t, dv in self.burns
            ],
        }
        return json.dumps(data, separators=(",", ":"))


def run(
    scenario: Scenario,
    controller: Controller,
    specs: tuple[spec_mod.Spec, ...] = (),
) -> Trace:
    """Run the closed loop; returns the recorded Trace with spec margins."""
    sc = scenario
    n_ticks = int(round(sc.duration / sc.dt_gnc))
    rng = np.random.default_rng(sc.seed)
    f = nl._deriv(sc.cfg, sc.bc_target, sc.bc_chaser)
    h = sc.dt_gnc / sc.truth_substeps

    rv_chaser0 = nl.lvlh_to_eci(sc.rv_target0, sc.x_rel0, sc.cfg, sc.bc_target)
    y = np.concatenate([sc.rv_target0, rv_chaser0])

    times = np.zeros(n_ticks + 1)
    x_rel = np.zeros((n_ticks + 1, 6))
    rv_t = np.zeros((n_ticks + 1, 6))
    burns: list[tuple[float, F64]] = []

    for k in range(n_ticks + 1):
        t = k * sc.dt_gnc
        times[k] = t
        rel = nl.eci_to_lvlh(y[0:6], y[6:12], sc.cfg, sc.bc_target)
        x_rel[k] = rel
        rv_t[k] = y[0:6]
        if k == n_ticks:
            break

        meas = rel.copy()
        if sc.meas_pos_std > 0.0:
            meas[0:3] += rng.normal(0.0, sc.meas_pos_std, 3)
        if sc.meas_vel_std > 0.0:
            meas[3:6] += rng.normal(0.0, sc.meas_vel_std, 3)

        dv = np.asarray(controller(t, meas), dtype=np.float64)
        if float(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]) > 0.0:
            # actuator hardware: magnitude cap, MIB quantization, then
            # proportional per-axis execution error (in that order)
            mag = float(np.linalg.norm(dv))
            if mag > sc.dv_max_tick:
                dv = dv * (sc.dv_max_tick / mag)
            if sc.dv_quantum > 0.0:
                dv = np.round(dv / sc.dv_quantum) * sc.dv_quantum
            if sc.dv_exec_std_frac > 0.0:
                dv = dv * (1.0 + rng.normal(0.0, sc.dv_exec_std_frac, 3))
        if float(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]) > 0.0:
            # Impulse commanded in the rotating LVLH frame: position is
            # unchanged, so the frame terms drop and the ECI velocity jump
            # is the rotated dv.
            rot = nl.lvlh_rotation(y[0:3], y[3:6])
            y[9:12] = y[9:12] + rot.T @ dv
            burns.append((t, dv.copy()))

        for i in range(sc.truth_substeps):
            y = integrators.rk4_step(f, t + i * h, y, h)

    trace = Trace(times, x_rel, rv_t, burns, {})
    if specs:
        trace.spec_margins = spec_mod.evaluate(specs, trace.channels())
    return trace


def circular_target(a: float, inc: float = 0.9, raan: float = 0.5,
                    argp: float = 1.2, nu: float = 0.0) -> F64:
    """Convenience: near-circular target ECI state for scenarios."""
    r, v = nl.elements_to_rv(a, 0.0, inc, raan, argp, nu, const.MU_EARTH)
    return np.concatenate([r, v])


def mean_motion_of(rv_target: F64) -> float:
    """Osculating mean motion from an ECI state (two-body)."""
    r = float(np.linalg.norm(rv_target[0:3]))
    v2 = float(np.dot(rv_target[3:6], rv_target[3:6]))
    a = 1.0 / (2.0 / r - v2 / const.MU_EARTH)
    return math.sqrt(const.MU_EARTH / (a * a * a))
