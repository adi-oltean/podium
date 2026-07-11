"""Thruster allocation: realize a commanded body wrench with a cluster
of discrete, push-only thrusters.

The 6-DOF guidance (podium.guidance.sixdof) commands a continuous body
wrench — thrust along +x_body plus a 3-axis torque. Hardware realizes
that wrench through N fixed thrusters, each producing force along a
fixed unit direction d_i applied at a fixed body position r_i, and each
able only to PUSH (magnitude u_i >= 0). The per-thruster wrench is

    [ F ]   [        d_i        ]
    [ T ] = [    r_i x d_i      ] u_i          (u_i >= 0)

so the cluster's effectiveness matrix B (6 x N) stacks those columns
and the realized wrench is B u. Allocation inverts this: find the
minimum-propellant non-negative u producing the demanded wrench.

Non-negativity is what makes this non-trivial — a plain pseudoinverse
happily returns negative "thrust", which no thruster can produce. The
LP (min sum u, B u = w, 0 <= u <= u_max) respects the physics and
minimizes propellant; the least-squares fallback reports the closest
achievable wrench when a demand exceeds the cluster's authority.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]


@dataclass(frozen=True)
class ThrusterConfig:
    """A cluster of push-only thrusters in the body frame."""

    positions: F64   # (N, 3) application points [m]
    directions: F64  # (N, 3) unit thrust directions (force is +d)

    def __post_init__(self) -> None:
        p = np.asarray(self.positions, dtype=np.float64)
        d = np.asarray(self.directions, dtype=np.float64)
        if p.shape != d.shape or p.ndim != 2 or p.shape[1] != 3:
            raise ValueError("positions/directions must both be (N, 3)")
        norms = np.linalg.norm(d, axis=1, keepdims=True)
        if np.any(norms == 0.0):
            raise ValueError("thruster directions must be nonzero")
        object.__setattr__(self, "positions", p)
        object.__setattr__(self, "directions", d / norms)

    @property
    def n(self) -> int:
        return int(self.positions.shape[0])

    def effectiveness(self) -> F64:
        """B (6 x N): column i = [d_i ; r_i x d_i]."""
        d = self.directions
        tau = np.cross(self.positions, d)
        return np.vstack([d.T, tau.T])


@dataclass
class Allocation:
    u: F64            # (N,) non-negative thrust magnitudes
    realized: F64     # (6,) wrench actually produced, B u
    residual: float   # ||realized - demand||
    feasible: bool    # exact (within tol) and within bounds
    propellant: float  # sum u


def allocate(cfg: ThrusterConfig, wrench: F64, u_max: float = np.inf,
             tol: float = 1e-7) -> Allocation:
    """Minimum-propellant non-negative allocation of a body wrench.

    Solves  min sum(u)  s.t.  B u = wrench,  0 <= u <= u_max  (LP).
    If the LP is infeasible (the demand exceeds the cluster's authority
    or bounds), falls back to a bounded non-negative least-squares fit
    and returns feasible=False with the achievable wrench + residual.
    """
    from scipy.optimize import linprog, lsq_linear  # type: ignore[import-untyped]

    b = cfg.effectiveness()
    w = np.asarray(wrench, dtype=np.float64)
    n = cfg.n
    bounds = [(0.0, None if not np.isfinite(u_max) else u_max)] * n
    res = linprog(c=np.ones(n), A_eq=b, b_eq=w, bounds=bounds,
                  method="highs")
    if res.success:
        u = np.clip(res.x, 0.0, None)
        realized = b @ u
        r = float(np.linalg.norm(realized - w))
        return Allocation(u=u, realized=realized, residual=r,
                          feasible=r <= tol, propellant=float(u.sum()))

    # infeasible: closest achievable wrench subject to 0 <= u <= u_max, a
    # TRUE box-bounded least-squares fit -- not an unbounded NNLS clipped to
    # u_max, which need not be the closest feasible point. When u_max is inf
    # this reduces to non-negative least squares.
    sol = lsq_linear(b, w, bounds=(0.0, u_max))
    u = np.clip(sol.x, 0.0, u_max if np.isfinite(u_max) else None)
    realized = b @ u
    r = float(np.linalg.norm(realized - w))
    return Allocation(u=u, realized=realized, residual=r, feasible=False,
                      propellant=float(u.sum()))


def standard_cluster(half: float = 1.0) -> ThrusterConfig:
    """A 24-thruster cluster with full 6-DOF wrench authority (verified
    rank-6 B). At each of the 8 corners of a cube (edge 2*half), three
    thrusters fire INWARD along the body axes (a corner at (+,+,+) has
    thrusters along -x, -y, -z). Opposing corners supply both signs of
    every force component, and the corner offsets supply both signs of
    every torque — so any small wrench, including a pure couple, is
    realizable with push-only (u >= 0) commands."""
    corners = np.array([[sx, sy, sz]
                        for sx in (-1.0, 1.0)
                        for sy in (-1.0, 1.0)
                        for sz in (-1.0, 1.0)]) * half
    pos, dirs = [], []
    for c in corners:
        for axis in range(3):
            d = [0.0, 0.0, 0.0]
            d[axis] = -float(np.sign(c[axis]))  # inward along this axis
            pos.append(c)
            dirs.append(d)
    return ThrusterConfig(positions=np.array(pos),
                          directions=np.array(dirs))
