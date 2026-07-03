"""Layer-0 convex guidance: impulsive planning on exact STM discretizations.

Sandbox layer: this module may use the full Python stack (cvxpy). It is
the prototyping half of the trajectory-optimization design; the embedded
half (generated fixed-iteration solvers) arrives with the C emitter.
Problems are DPP-compiled once per (grid, constraint topology): STMs,
boundary conditions, and keep-out-zone normals are cvxpy Parameters, so
repeated solves (MPC-style replanning, KOZ reference refinement, Monte
Carlo) never rebuild or re-canonicalize the problem.

Dynamics are exact linear equality constraints — the CW and
Yamanaka-Ankersen STMs discretize the relative motion with zero
integration error, so the only model error is linearization about the
target orbit (quantified in tests/test_validity_envelopes.py).

Keep-out-zone handling follows Mueller & Larsson: linear hyperplane
constraints n_k^T r_k >= R with unit normals n_k taken from a reference
trajectory. Since ||r|| >= n^T r for unit n, satisfaction implies true
distance >= R (conservative). Normals default to a straight-line
reference and are refined by one fixed re-solve pass using the previous
solution — a bounded two-pass scheme, not an open iteration.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import cvxpy as cp
import numpy as np
from numpy.typing import NDArray

from podium.core import cw, ya
from podium.core import roe as roe_mod

F64 = NDArray[np.float64]

_B = np.zeros((6, 3))
_B[3:6, :] = np.eye(3)


@dataclass(frozen=True)
class ConeSpec:
    """Approach cone: for active nodes, positions must satisfy
    ||P_perp (r - apex)|| <= tan(half_angle) * axis.(r - apex)."""

    apex: tuple[float, float, float]
    axis: tuple[float, float, float]  # unit vector, opening direction
    half_angle: float  # rad
    from_time: float = 0.0  # active for node times >= from_time


@dataclass(frozen=True)
class KozSpec:
    """Spherical keep-out zone enforced via rotating hyperplanes."""

    radius: float
    until_time: float = math.inf  # active for node times <= until_time


@dataclass(frozen=True)
class PlumeSpec:
    """No thrusting toward the target when closer than range_m: for active
    nodes, d_k . dv_k <= 0 with d_k the reference direction to target."""

    range_m: float


@dataclass
class Plan:
    """Impulsive plan on the node grid (burns at every node, incl. arrival)."""

    times: F64  # (K+1,)
    dvs: F64  # (K+1, 3)
    states: F64  # (K+1, 6) pre-burn states
    objective: float
    status: str
    # (node index, unit direction toward target) for each plume constraint
    # that was active in the final solve — recorded for transparency/tests.
    plume_dirs: list[tuple[int, F64]] = field(default_factory=list)

    def total_dv(self) -> float:
        return float(np.sum(np.linalg.norm(self.dvs, axis=1)))


def plan_to_controller(
    plan: Plan, tick_tol: float = 0.5
) -> "Callable[[float, F64], F64]":
    """Adapter: engine controller firing each planned burn once at the
    first tick within tick_tol of its scheduled time."""
    fired = [False] * len(plan.times)

    def controller(t: float, _meas: F64) -> F64:
        for i, tb in enumerate(plan.times):
            if not fired[i] and abs(t - tb) <= tick_tol:
                fired[i] = True
                return plan.dvs[i]
        return np.zeros(3)

    return controller


class RendezvousPlanner:
    """DPP-compiled impulsive rendezvous planner on CW/YA STM dynamics.

    Build once for a node grid and constraint topology; solve many times
    with different boundary conditions / models via Parameters.
    """

    def __init__(
        self,
        times: F64,
        objective: str = "l2",
        cone: ConeSpec | None = None,
        koz: KozSpec | None = None,
        plume: PlumeSpec | None = None,
        dv_max: float | None = None,
    ) -> None:
        self.times = np.asarray(times, dtype=np.float64)
        k = len(self.times) - 1
        if k < 1:
            raise ValueError("need at least two nodes")
        self.k = k
        self.cone = cone
        self.koz = koz
        self.plume = plume

        x = cp.Variable((6, k + 1), name="x")
        v = cp.Variable((3, k + 1), name="dv")
        self._x, self._v = x, v
        self._x0 = cp.Parameter(6, name="x0")
        self._xf = cp.Parameter(6, name="xf")
        # Phi_k and Phi_k @ B as parameters keeps the problem DPP while
        # letting the dynamics model (CW/YA, any n/e) vary between solves.
        self._phis = [cp.Parameter((6, 6), name=f"phi{i}") for i in range(k)]
        self._phibs = [cp.Parameter((6, 3), name=f"phib{i}") for i in range(k)]

        cons = [x[:, 0] == self._x0]
        for i in range(k):
            cons.append(
                x[:, i + 1] == self._phis[i] @ x[:, i] + self._phibs[i] @ v[:, i]
            )
        cons.append(x[:, k] + _B @ v[:, k] == self._xf)

        if dv_max is not None:
            for i in range(k + 1):
                cons.append(cp.norm(v[:, i], 2) <= dv_max)

        if cone is not None:
            axis = np.asarray(cone.axis) / np.linalg.norm(cone.axis)
            proj = np.eye(3) - np.outer(axis, axis)
            apex = np.asarray(cone.apex)
            tan_a = math.tan(cone.half_angle)
            self.cone_nodes = [
                i for i, t in enumerate(self.times) if t >= cone.from_time
            ]
            for i in self.cone_nodes:
                r = x[0:3, i] - apex
                cons.append(cp.norm(proj @ r, 2) <= tan_a * (axis @ r))

        self.koz_nodes: list[int] = []
        self._koz_normals: list[cp.Parameter] = []
        if koz is not None:
            self.koz_nodes = [
                i for i, t in enumerate(self.times)
                if t <= koz.until_time and 0 < i < k + 1
            ]
            for i in self.koz_nodes:
                nrm = cp.Parameter(3, name=f"koz_n{i}")
                self._koz_normals.append(nrm)
                cons.append(nrm @ x[0:3, i] >= koz.radius)

        self.plume_nodes: list[int] = []
        self._plume_dirs: list[cp.Parameter] = []
        self._plume_active: list[tuple[int, F64]] = []
        if plume is not None:
            # Activity decided per solve from the reference; compile a
            # direction parameter per node and neutralize inactive ones by
            # setting the direction to zero. The arrival node is exempt:
            # the braking burn at contact cannot honor a plume half-space
            # (real vehicles brake with lateral/canted thrusters there).
            self.plume_nodes = list(range(k))
            for i in self.plume_nodes:
                d = cp.Parameter(3, name=f"plume_d{i}")
                self._plume_dirs.append(d)
                cons.append(d @ v[:, i] <= 0.0)

        if objective == "l2":
            obj = cp.sum([cp.norm(v[:, i], 2) for i in range(k + 1)])
        elif objective == "l1":
            obj = cp.sum([cp.norm(v[:, i], 1) for i in range(k + 1)])
        else:
            raise ValueError("objective must be 'l1' or 'l2'")
        self._problem = cp.Problem(cp.Minimize(obj), cons)
        assert self._problem.is_dcp(dpp=True), "planner must stay DPP"

    # -- model loading -------------------------------------------------
    def _set_stms(self, n: float, e: float, theta0: float) -> None:
        for i in range(self.k):
            dt = float(self.times[i + 1] - self.times[i])
            if e == 0.0:
                phi = cw.stm(n, dt)
            else:
                th_i = ya.propagate_true_anomaly(n, e, theta0, float(self.times[i]))
                phi = ya.stm(n, e, th_i, dt)
            self._phis[i].value = phi
            self._phibs[i].value = phi @ _B

    def _reference_positions(self, x0: F64, xf: F64) -> F64:
        """Straight-line position reference for constraint normals."""
        ref = np.zeros((self.k + 1, 3))
        for i in range(self.k + 1):
            frac = i / self.k
            ref[i] = (1 - frac) * x0[0:3] + frac * xf[0:3]
        return ref

    def _set_reference(self, ref_positions: F64) -> None:
        for nrm_param, i in zip(self._koz_normals, self.koz_nodes):
            r = ref_positions[i]
            mag = float(np.linalg.norm(r))
            nrm_param.value = r / mag if mag > 1e-9 else np.array([1.0, 0.0, 0.0])
        if self.plume is not None:
            self._plume_active = []
            for d_param, i in zip(self._plume_dirs, self.plume_nodes):
                r = ref_positions[i]
                mag = float(np.linalg.norm(r))
                if 1e-9 < mag < self.plume.range_m:
                    d = -r / mag  # toward the target
                    d_param.value = d
                    self._plume_active.append((i, d))
                else:
                    d_param.value = np.zeros(3)  # inactive node

    def solve(
        self,
        x0: F64,
        xf: F64,
        n: float,
        e: float = 0.0,
        theta0: float = 0.0,
        refine_passes: int = 1,
    ) -> Plan:
        """Solve; with KOZ/plume, run `refine_passes` extra solves whose
        constraint normals come from the previous solution (bounded, not
        an open iteration)."""
        self._set_stms(n, e, theta0)
        self._x0.value = np.asarray(x0, dtype=np.float64)
        self._xf.value = np.asarray(xf, dtype=np.float64)
        self._set_reference(self._reference_positions(x0, xf))
        self._problem.solve(solver=cp.CLARABEL)
        passes = refine_passes if (self.koz or self.plume) else 0
        for _ in range(passes):
            assert self._x.value is not None
            self._set_reference(self._x.value[0:3, :].T.copy())
            self._problem.solve(solver=cp.CLARABEL)
        assert self._v.value is not None and self._x.value is not None
        return Plan(
            times=self.times.copy(),
            dvs=self._v.value.T.copy(),
            states=self._x.value.T.copy(),
            objective=float(self._problem.value),
            status=str(self._problem.status),
            plume_dirs=list(self._plume_active),
        )


@dataclass
class RoePlan:
    times: F64
    dvs: F64  # (K+1, 3) RTN burns
    roes: F64  # (K+1, 6) pre-burn ROE states
    objective: float
    status: str

    def total_dv(self) -> float:
        return float(np.sum(np.linalg.norm(self.dvs, axis=1)))


class RoePlanner:
    """Impulsive ROE reconfiguration on Keplerian/J2 STM dynamics.

    States live in ROE space; burns act through the control-input matrix
    at each node's argument of latitude. The natural problem class for
    formation acquisition/reconfiguration and safety-geometry shaping.
    """

    def __init__(self, times: F64, objective: str = "l2") -> None:
        self.times = np.asarray(times, dtype=np.float64)
        self.k = len(self.times) - 1
        r = cp.Variable((6, self.k + 1), name="roe")
        v = cp.Variable((3, self.k + 1), name="dv")
        self._r, self._v = r, v
        self._r0 = cp.Parameter(6, name="roe0")
        self._rf = cp.Parameter(6, name="roef")
        self._phis = [cp.Parameter((6, 6), name=f"phi{i}") for i in range(self.k)]
        self._phigs = [cp.Parameter((6, 3), name=f"phig{i}") for i in range(self.k)]
        self._gf = cp.Parameter((6, 3), name="gamma_final")

        cons = [r[:, 0] == self._r0]
        for i in range(self.k):
            cons.append(
                r[:, i + 1] == self._phis[i] @ r[:, i] + self._phigs[i] @ v[:, i]
            )
        cons.append(r[:, self.k] + self._gf @ v[:, self.k] == self._rf)
        if objective == "l2":
            obj = cp.sum([cp.norm(v[:, i], 2) for i in range(self.k + 1)])
        else:
            obj = cp.sum([cp.norm(v[:, i], 1) for i in range(self.k + 1)])
        self._problem = cp.Problem(cp.Minimize(obj), cons)
        assert self._problem.is_dcp(dpp=True)

    def solve(
        self,
        roe0: F64,
        roef: F64,
        a: float,
        n: float,
        u0: float,
        mu: float = 0.0,
        j2: float = 0.0,
        r_body: float = 0.0,
        e: float = 0.0,
        inc: float = 1.0,
        argp: float = 0.0,
    ) -> RoePlan:
        for i in range(self.k):
            dt = float(self.times[i + 1] - self.times[i])
            if j2 > 0.0:
                phi = roe_mod.stm_j2(mu, j2, r_body, a, e, inc, argp, dt)
            else:
                phi = roe_mod.stm_keplerian(n, dt)
            u_i = u0 + n * float(self.times[i])
            gamma = roe_mod.control_matrix(a, n, u_i)
            self._phis[i].value = phi
            self._phigs[i].value = phi @ gamma
        self._gf.value = roe_mod.control_matrix(
            a, n, u0 + n * float(self.times[self.k])
        )
        self._r0.value = np.asarray(roe0, dtype=np.float64)
        self._rf.value = np.asarray(roef, dtype=np.float64)
        self._problem.solve(solver=cp.CLARABEL)
        assert self._v.value is not None and self._r.value is not None
        return RoePlan(
            times=self.times.copy(),
            dvs=self._v.value.T.copy(),
            roes=self._r.value.T.copy(),
            objective=float(self._problem.value),
            status=str(self._problem.status),
        )
