"""Terminal guidance to a tumbling target (scoped study, sandbox).

Structural point this module demonstrates: for a KNOWN tumble state,
the dock port's position and velocity at any time are deterministic
kinematics — so terminal port capture is an ordinary boundary condition
and the rotating approach corridor is a per-node second-order cone with
a known axis. The whole planning problem stays CONVEX on the exact CW
STM; successive convexification only becomes necessary once the tumble
itself is uncertain or torque-coupled (recorded as follow-on scope in
docs/plans/23-tumbling-study.md).

Tumble model: planar rotation about the cross-track (z) axis at rate
w_spin, port at radius rho_p from the target's center of mass:
    p(t) = rho_p [cos(w t + phi), sin(w t + phi), 0]
    pdot(t) = rho_p w [-sin(w t + phi), cos(w t + phi), 0]
The chaser must arrive at p(t_f) WITH velocity pdot(t_f) (grapple/berth
condition) having stayed inside the rotating corridor cone around the
port direction for the final approach nodes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import cvxpy as cp
import numpy as np
from numpy.typing import NDArray

from podium.core import cw

F64 = NDArray[np.float64]

_B = np.zeros((6, 3))
_B[3:6, :] = np.eye(3)


@dataclass(frozen=True)
class Tumble:
    """Planar tumble about +z at w_spin [rad/s]; port radius/phase."""

    w_spin: float
    rho_port: float
    phase0: float = 0.0

    def port_state(self, t: float) -> F64:
        """[position(3), velocity(3)] of the port in LVLH at time t."""
        a = self.w_spin * t + self.phase0
        c, s = math.cos(a), math.sin(a)
        return np.array([
            self.rho_port * c, self.rho_port * s, 0.0,
            -self.rho_port * self.w_spin * s,
            self.rho_port * self.w_spin * c, 0.0,
        ])


@dataclass
class TumblingPlan:
    times: F64
    dvs: F64
    states: F64
    status: str
    terminal_pos_err: float
    terminal_vel_err: float

    def total_dv(self) -> float:
        return float(np.sum(np.linalg.norm(self.dvs, axis=1)))


def plan_tumbling_dock(
    times: F64,
    x0: F64,
    tumble: Tumble,
    n: float,
    corridor_half_angle: float = math.radians(20.0),
    corridor_from_frac: float = 0.5,
    dv_max: float | None = None,
) -> TumblingPlan:
    """Convex rotating-corridor planner: exact CW dynamics, terminal
    match of the rotating port state, per-node cones about the KNOWN
    port direction for the final approach segment."""
    times = np.asarray(times, dtype=np.float64)
    k = len(times) - 1
    xf = tumble.port_state(float(times[-1]))

    x = cp.Variable((6, k + 1))
    v = cp.Variable((3, k + 1))
    cons = [x[:, 0] == x0]
    for i in range(k):
        phi = cw.stm(n, float(times[i + 1] - times[i]))
        cons.append(x[:, i + 1] == phi @ x[:, i] + (phi @ _B) @ v[:, i])
    cons.append(x[:, k] + _B @ v[:, k] == xf)
    if dv_max is not None:
        for i in range(k + 1):
            cons.append(cp.norm(v[:, i], 2) <= dv_max)

    tan_a = math.tan(corridor_half_angle)
    from_idx = max(1, int(corridor_from_frac * k))
    for i in range(from_idx, k + 1):
        p = tumble.port_state(float(times[i]))[0:3]
        axis = p / np.linalg.norm(p)  # outward port direction
        proj = np.eye(3) - np.outer(axis, axis)
        # chaser must sit inside the cone opening OUTWARD from the port
        r_rel = x[0:3, i] - p
        cons.append(cp.norm(proj @ r_rel, 2) <= tan_a * (axis @ r_rel)
                    + 1e-6)

    fuel = cp.sum([cp.norm(v[:, i], 2) for i in range(k + 1)])
    prob = cp.Problem(cp.Minimize(fuel), cons)
    prob.solve(solver=cp.CLARABEL)
    if x.value is None:
        return TumblingPlan(times.copy(), np.zeros((k + 1, 3)),
                            np.zeros((k + 1, 6)), str(prob.status),
                            math.inf, math.inf)
    states = x.value.T.copy()
    dvs = v.value.T.copy()
    xf_reached = states[k] + np.concatenate([np.zeros(3), dvs[k]])
    return TumblingPlan(
        times.copy(), dvs, states, str(prob.status),
        float(np.linalg.norm(xf_reached[0:3] - xf[0:3])),
        float(np.linalg.norm(xf_reached[3:6] - xf[3:6])),
    )


def envelope_sweep(
    rates: F64,
    times: F64,
    x0: F64,
    n: float,
    rho_port: float = 10.0,
    dv_max: float | None = 0.5,
    fix_arrival_phase: float | None = None,
) -> list[dict]:
    """Fuel/feasibility vs tumble rate — the study's headline table.

    STUDY FINDING: with a free phase, fuel vs rate is NON-monotone —
    the port's orientation at arrival (phase0 + w*t_f) changes with the
    rate and arrival geometry dominates the co-rotation cost. For a
    clean envelope, fix the ARRIVAL phase per rate
    (phase0 = fix_arrival_phase - w*t_f); then dv isolates the true
    rate cost: co-rotation velocity rho*w plus corridor chasing."""
    t_f = float(np.asarray(times)[-1])
    out = []
    for w in rates:
        phase0 = 0.0 if fix_arrival_phase is None \
            else fix_arrival_phase - float(w) * t_f
        tumble = Tumble(w_spin=float(w), rho_port=rho_port, phase0=phase0)
        plan = plan_tumbling_dock(times, x0, tumble, n, dv_max=dv_max)
        feasible = plan.status == "optimal" and plan.terminal_pos_err < 1e-3
        out.append({
            "w_spin": float(w),
            "feasible": feasible,
            "dv": plan.total_dv() if feasible else math.inf,
        })
    return out
