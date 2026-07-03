"""ARCH-COMP spacecraft rendezvous benchmark: model, simulation, export.

The community-standard closed-loop verification case (Chan & Mitra,
ARCH17): planar CW relative motion about a GEO target with a switched
controller — mode 1 "approaching" (x <= -100 m), mode 2 "rendezvous
attempt" (inside the 100 m box), mode 3 "aborting" (passive, after the
abort time). Safety properties: inside mode 2 the chaser stays in the
30-degree line-of-sight cone with speed below 0.055 m/s (an octagonal
overapproximation of the speed disc); after an abort it never enters the
0.2 m target box.

Units follow the benchmark: positions in meters, velocities in meters
per MINUTE, time in minutes. State vector: [x, y, vx, vy, t] — the clock
is part of the state so abort transitions are time-triggered guards.

Matrix provenance: the closed-loop mode matrices are the published
benchmark values (ARCH-COMP repeatability packages). The abort mode is
exactly planar CW at GEO mean motion: its matrix entries are
[3n^2, 2n; -2n, 0] with n = 0.00438138 rad/min — asserted against
podium.core.cw in the tests, tying the benchmark to our kernel.

This module is sandbox-side. `export_model()` writes the hybrid
automaton as JSON; `tools/reach/arch_rendezvous.jl` consumes it and
re-proves the properties with ReachabilityAnalysis.jl — Podium exports
the model, the reachability tool verifies it, CI gates on the result.
"""

from __future__ import annotations

import json
import math

import numpy as np
from numpy.typing import NDArray

from podium.control import lqr as _lqr
from podium.core import cw as _cw

F64 = NDArray[np.float64]

# mean motion of the GEO target orbit [rad/min], from the benchmark's
# abort-mode matrix (2n = 0.00876276)
N_RAD_MIN = 0.00438138

# --- closed-loop mode matrices (5-state: x, y, vx, vy, t) --------------
# mode 1: "approaching" (LQR gains K1 folded in, published values)
A_APPROACH = np.array([
    [0.0, 0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 0.0],
    [-0.057599765881773, 0.000200959896519766, -2.89995083970656,
     0.00877200894463775, 0.0],
    [-0.000174031357370456, -0.0665123984901026, -0.00875351105536225,
     -2.90300269286856, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0],
])
# mode 2: "rendezvous attempt"
A_ATTEMPT = np.array([
    [0.0, 0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 0.0],
    [-0.575999943070835, 0.000262486079431672, -19.2299795908647,
     0.00876275931760007, 0.0],
    [-0.000262486080737868, -0.575999940191886, -0.00876276068239993,
     -19.2299765959399, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0],
])
# mode 3: "aborting" — passive planar CW at GEO
A_ABORT = np.array([
    [0.0, 0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 0.0],
    [3.0 * N_RAD_MIN**2, 0.0, 0.0, 2.0 * N_RAD_MIN, 0.0],
    [0.0, 0.0, -2.0 * N_RAD_MIN, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0],
])
_B_CLOCK = np.array([0.0, 0.0, 0.0, 0.0, 1.0])  # dt/dt = 1

MODE_NAMES = ("approaching", "attempt", "aborting")
_MATS = (A_APPROACH, A_ATTEMPT, A_ABORT)

# initial set (mode 1): center +/- radius, velocities and clock exact
X0_CENTER = np.array([-900.0, -400.0, 0.0, 0.0, 0.0])
X0_RADIUS = np.array([25.0, 25.0, 0.0, 0.0, 0.0])

TAN30 = math.tan(math.pi / 6.0)
V_MAX = 0.055 * 60.0  # velocity ceiling [m/min] = 3.3
_CX = V_MAX * math.cos(math.pi / 8.0)
_CY = V_MAX * math.sin(math.pi / 8.0)
TARGET_HALF_WIDTH = 0.2  # abort keep-out box [m]
HORIZON = 300.0  # minutes


def _in_attempt_box(s: F64) -> bool:
    """Guard/invariant region of the rendezvous-attempt mode (octagonal)."""
    x, y = float(s[0]), float(s[1])
    return (
        -100.0 <= x <= 100.0
        and -100.0 <= y <= 100.0
        and abs(x + y) <= 141.1
        and abs(x - y) <= 141.1
    )


# --- Podium-synthesized controller variant ----------------------------
# Plant: planar CW at the benchmark's GEO mean motion, built from
# podium.core.cw (never transcribed); inputs are accelerations along
# x, y. Gains from podium.control.lqr.clqr with the Q/R below — chosen
# so the closed-loop bandwidth/damping land in the same regime the
# published (verified) controller occupies, but the numbers are OURS.
_B_ACCEL = np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
Q_APPROACH = np.diag([0.0033, 0.0044, 8.3, 8.3])
Q_ATTEMPT = np.diag([0.33, 0.33, 369.0, 369.0])
R_CTRL = np.eye(2)


def cw_planar(n: float) -> F64:
    """Planar CW system matrix (columns via podium.core.cw.cw_deriv)."""
    a = np.zeros((4, 4))
    for j, e in enumerate(np.eye(4)):
        full = np.array([e[0], e[1], 0.0, e[2], e[3], 0.0])
        d = _cw.cw_deriv(full, n, np.zeros(3))
        a[:, j] = [d[0], d[1], d[3], d[4]]
    return a


def podium_gains() -> tuple[F64, F64]:
    """(K_approach, K_attempt) acceleration gains from Podium's CARE."""
    a = cw_planar(N_RAD_MIN)
    k1 = _lqr.clqr(a, _B_ACCEL, Q_APPROACH, R_CTRL)
    k2 = _lqr.clqr(a, _B_ACCEL, Q_ATTEMPT, R_CTRL)
    return k1, k2


def podium_mode_matrices() -> tuple[F64, F64, F64]:
    """Closed-loop 5-state matrices with Podium-synthesized gains.

    The abort mode is passive CW in both variants by construction.
    """
    a = cw_planar(N_RAD_MIN)
    k1, k2 = podium_gains()
    out = []
    for k in (k1, k2):
        acl = np.zeros((5, 5))
        acl[:4, :4] = a - _B_ACCEL @ k
        out.append(acl)
    return out[0], out[1], A_ABORT.copy()


def implied_reference_gains() -> tuple[F64, F64]:
    """Acceleration gains implied by the published closed loops
    (A_published - A_cw), for comparison/documentation."""
    a = cw_planar(N_RAD_MIN)
    k1 = -(A_APPROACH[2:4, :4] - a[2:4, :4])
    k2 = -(A_ATTEMPT[2:4, :4] - a[2:4, :4])
    return k1, k2


def _matrices(gains: str) -> tuple[F64, F64, F64]:
    if gains == "reference":
        return A_APPROACH, A_ATTEMPT, A_ABORT
    if gains == "podium":
        return podium_mode_matrices()
    raise ValueError("gains must be 'reference' or 'podium'")


def simulate(
    x0: F64,
    abort_time: float = -1.0,
    dt: float = 0.01,
    horizon: float = HORIZON,
    gains: str = "reference",
) -> tuple[F64, F64, NDArray[np.int64]]:
    """Deterministic hybrid simulation (RK4, urgent guard semantics).

    abort_time < 0 disables the abort (SRNA scenario); guards are checked
    on the fixed grid, so switching times are quantized to dt — fine for
    spec-margin checks, and the reachability tool owns the exact story.
    Returns (times, states (M,5), modes (M,)) with modes in {1,2,3}.
    """
    mats = _matrices(gains)
    steps = int(round(horizon / dt))
    times = np.zeros(steps + 1)
    states = np.zeros((steps + 1, 5))
    modes = np.zeros(steps + 1, dtype=np.int64)
    s = np.asarray(x0, dtype=np.float64).copy()
    mode = 1
    for k in range(steps + 1):
        t = k * dt
        # transitions (urgent, abort dominates)
        if abort_time >= 0.0 and mode != 3 and t >= abort_time:
            mode = 3
        elif mode == 1 and _in_attempt_box(s):
            mode = 2
        times[k] = t
        states[k] = s
        modes[k] = mode
        if k == steps:
            break
        a = mats[mode - 1]

        def f(y: F64) -> F64:
            out: F64 = a @ y + _B_CLOCK
            return out

        k1 = f(s)
        k2 = f(s + 0.5 * dt * k1)
        k3 = f(s + 0.5 * dt * k2)
        k4 = f(s + dt * k3)
        s = s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return times, states, modes


def spec_margins(states: F64, modes: NDArray[np.int64]) -> dict[str, float]:
    """Worst-case margins of the three benchmark properties (>0 = satisfied).

    - los_cone: in mode 2, x >= -100 and |y| <= -x*tan(30 deg);
    - velocity: in mode 2, (vx, vy) inside the 3.3 m/min octagon;
    - abort_avoidance: in mode 3, outside the 0.2 m target box.
    Modes that never occur report +inf (vacuously satisfied, visibly so).
    """
    inf = math.inf
    los, vel, avoid = inf, inf, inf
    for s, m in zip(states, modes):
        x, y, vx, vy = s[0], s[1], s[2], s[3]
        if m == 2:
            los = min(los, x + 100.0, -TAN30 * x - y, -TAN30 * x + y)
            vel = min(
                vel,
                _CX - abs(vx),
                _CX - abs(vy),
                (_CX + _CY) - abs(vx + vy),
                (_CX + _CY) - abs(vx - vy),
            )
        elif m == 3:
            avoid = min(avoid, max(abs(x), abs(y)) - TARGET_HALF_WIDTH)
    return {"los_cone": los, "velocity": vel, "abort_avoidance": avoid}


def initial_corners() -> list[F64]:
    """Center + the four position corners of the initial hyperrectangle."""
    pts = [X0_CENTER.copy()]
    for sx in (-1.0, 1.0):
        for sy in (-1.0, 1.0):
            p = X0_CENTER.copy()
            p[0] += sx * X0_RADIUS[0]
            p[1] += sy * X0_RADIUS[1]
            pts.append(p)
    return pts


def export_model(
    abort_time: float = -1.0, gains: str = "reference"
) -> dict[str, object]:
    """Hybrid-automaton export for the external reachability tool.

    Halfspaces are (a, b) meaning a . s <= b over the 5-state vector.
    gains: "reference" (published benchmark controller) or "podium"
    (gains synthesized by podium.control.lqr.clqr — same FSM, guards,
    initial set, and properties).
    """
    def hs(idx: list[int], coef: list[float], b: float) -> dict[str, object]:
        a = [0.0] * 5
        for i, c in zip(idx, coef):
            a[i] = c
        return {"a": a, "b": b}

    attempt_box = [
        hs([0], [-1.0], 100.0), hs([0], [1.0], 100.0),
        hs([1], [-1.0], 100.0), hs([1], [1.0], 100.0),
        hs([0, 1], [-1.0, -1.0], 141.1), hs([0, 1], [1.0, 1.0], 141.1),
        hs([0, 1], [1.0, -1.0], 141.1), hs([0, 1], [-1.0, 1.0], 141.1),
    ]
    aborting = abort_time >= 0.0
    inv1 = [hs([0], [1.0], -100.0)]
    inv2 = list(attempt_box)
    if aborting:
        inv1.append(hs([4], [1.0], abort_time))
        inv2.append(hs([4], [1.0], abort_time))

    a1, a2, a3 = _matrices(gains)
    modes = [
        {"name": "approaching", "A": a1.tolist(),
         "b": _B_CLOCK.tolist(), "invariant": inv1},
        {"name": "attempt", "A": a2.tolist(),
         "b": _B_CLOCK.tolist(), "invariant": inv2},
    ]
    transitions = [{"from": 1, "to": 2, "guard": attempt_box}]
    if aborting:
        modes.append({"name": "aborting", "A": a3.tolist(),
                      "b": _B_CLOCK.tolist(), "invariant": []})
        abort_guard = [hs([4], [-1.0], -abort_time)]
        transitions.append({"from": 1, "to": 3, "guard": abort_guard})
        transitions.append({"from": 2, "to": 3, "guard": abort_guard})

    return {
        "name": f"ARCH spacecraft rendezvous ({gains} gains)",
        "units": {"position": "m", "velocity": "m/min", "time": "min"},
        "state": ["x", "y", "vx", "vy", "t"],
        "horizon": HORIZON,
        "abort_time": abort_time,
        "modes": modes,
        "transitions": transitions,
        "initial": {"mode": 1, "center": X0_CENTER.tolist(),
                    "radius": X0_RADIUS.tolist()},
        "properties": {
            "attempt_mode": 2,
            "abort_mode": 3 if aborting else None,
            "tan30": TAN30,
            "v_octagon_cx": _CX,
            "v_octagon_cy": _CY,
            "target_half_width": TARGET_HALF_WIDTH,
        },
    }


def write_model(
    path: str, abort_time: float = -1.0, gains: str = "reference"
) -> None:
    with open(path, "w") as fh:
        json.dump(export_model(abort_time, gains), fh, indent=1)
