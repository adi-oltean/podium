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

from podium.control import lqr as _lqr
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


@dataclass(frozen=True)
class PassiveSafetySpec:
    """Breger-How failure scenarios: if the thruster dies at a failure
    node, the free-drift trajectory from that node's (pre-burn) state must
    stay outside the keep-out sphere for the safety horizon. Drift states
    are linear in the decision variables, so each (failure node, drift
    sample) pair is one linear hyperplane constraint — the normal folded
    with the drift STM into a single parameter row (stays DPP).

    failure_nodes=None protects every interior node. Currently requires
    the CW model (e=0): drift propagation uses the CW STM.
    """

    radius: float
    horizon: float  # drift time protected after a failure [s]
    n_samples: int = 6
    failure_nodes: tuple[int, ...] | None = None


@dataclass(frozen=True)
class AnnulusSpec:
    """Thrust annulus for finite-burn planning: the engine is always on
    with acceleration magnitude in [rho_min, rho_max] (min-throttle
    engines cannot run below rho_min)."""

    rho_min: float  # m/s^2
    rho_max: float  # m/s^2


@dataclass(frozen=True)
class SafeSetSpec:
    """Convex terminal safe set for the ROE planner: relative e- and
    i-vectors aligned with `direction` within `cone_angle` and at least
    e_min/i_min long — a convex sufficient condition for e/i-vector
    separation. The exact RN-plane margin is verified post-solve with
    guidance.safety.rn_margin (the scan is the receipt, the cone is the
    constraint)."""

    direction: tuple[float, float]  # unit vector in the e/i plane
    e_min: float
    i_min: float
    cone_angle: float  # max misalignment [rad]
    da_tol: float = 0.0  # |relative semi-major axis| bound (drift control)


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
        passive_safety: PassiveSafetySpec | None = None,
        track_state_weight: float = 1.0,
        track_control_weight: float = 1.0,
    ) -> None:
        self.times = np.asarray(times, dtype=np.float64)
        k = len(self.times) - 1
        if k < 1:
            raise ValueError("need at least two nodes")
        self.k = k
        self.cone = cone
        self.koz = koz
        self.plume = plume
        self.passive_safety = passive_safety

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

        # Breger-How passive-safety scenarios: one folded parameter row
        # per (failure node, drift sample); c @ x_j >= radius.
        self._ps_rows: list[tuple[int, int, cp.Parameter]] = []
        self._drift_stms: list[F64] = []
        if passive_safety is not None:
            nodes = passive_safety.failure_nodes
            self.ps_nodes: list[int] = (
                list(nodes) if nodes is not None else list(range(1, k))
            )
            taus = np.linspace(0.0, passive_safety.horizon,
                               passive_safety.n_samples + 1)[1:]
            self.ps_taus: F64 = taus
            for j in self.ps_nodes:
                for i in range(len(taus)):
                    c = cp.Parameter(6, name=f"ps_{j}_{i}")
                    self._ps_rows.append((j, i, c))
                    cons.append(c @ x[:, j] >= passive_safety.radius)

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

        self._x_ref: cp.Parameter | None = None
        if objective == "l2":
            obj = cp.sum([cp.norm(v[:, i], 2) for i in range(k + 1)])
        elif objective == "l1":
            obj = cp.sum([cp.norm(v[:, i], 1) for i in range(k + 1)])
        elif objective == "qp_tracking":
            # MPC-style QP: quadratic tracking of a reference trajectory
            # (a Parameter — re-solves stay on the compiled problem) plus
            # control effort. Terminal equality still applies.
            self._x_ref = cp.Parameter((6, k + 1), name="x_ref")
            obj = track_control_weight * cp.sum_squares(v) \
                + track_state_weight * cp.sum_squares(x - self._x_ref)
        else:
            raise ValueError("objective must be 'l1', 'l2', or 'qp_tracking'")
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

    def _reference_states(self, x0: F64, xf: F64) -> F64:
        """Straight-line full-state reference for constraint normals."""
        ref = np.zeros((self.k + 1, 6))
        for i in range(self.k + 1):
            frac = i / self.k
            ref[i, 0:3] = (1 - frac) * x0[0:3] + frac * xf[0:3]
        return ref

    def _set_reference(self, ref_states: F64) -> None:
        for nrm_param, i in zip(self._koz_normals, self.koz_nodes):
            r = ref_states[i, 0:3]
            mag = float(np.linalg.norm(r))
            nrm_param.value = r / mag if mag > 1e-9 else np.array([1.0, 0.0, 0.0])
        if self.plume is not None:
            self._plume_active = []
            for d_param, i in zip(self._plume_dirs, self.plume_nodes):
                r = ref_states[i, 0:3]
                mag = float(np.linalg.norm(r))
                if 1e-9 < mag < self.plume.range_m:
                    d = -r / mag  # toward the target
                    d_param.value = d
                    self._plume_active.append((i, d))
                else:
                    d_param.value = np.zeros(3)  # inactive node
        for j, i, c_param in self._ps_rows:
            phi = self._drift_stms[i]
            r_ref = (phi @ ref_states[j])[0:3]
            mag = float(np.linalg.norm(r_ref))
            nrm = r_ref / mag if mag > 1e-9 else np.array([1.0, 0.0, 0.0])
            c_param.value = nrm @ phi[0:3, :]

    def solve(
        self,
        x0: F64,
        xf: F64,
        n: float,
        e: float = 0.0,
        theta0: float = 0.0,
        refine_passes: int = 1,
        x_ref: F64 | None = None,
    ) -> Plan:
        """Solve; with KOZ/plume/passive-safety, run `refine_passes` extra
        solves whose constraint normals come from the previous solution
        (bounded, not an open iteration)."""
        self._set_stms(n, e, theta0)
        if self.passive_safety is not None:
            if e != 0.0:
                raise ValueError("passive-safety scenarios require e=0 (CW drift)")
            self._drift_stms = [cw.stm(n, float(t)) for t in self.ps_taus]
        if self._x_ref is not None:
            self._x_ref.value = (
                np.zeros((6, self.k + 1)) if x_ref is None else np.asarray(x_ref).T
            )
        self._x0.value = np.asarray(x0, dtype=np.float64)
        self._xf.value = np.asarray(xf, dtype=np.float64)
        self._set_reference(self._reference_states(x0, xf))
        self._problem.solve(solver=cp.CLARABEL)
        passes = (
            refine_passes
            if (self.koz or self.plume or self.passive_safety)
            else 0
        )
        for _ in range(passes):
            assert self._x.value is not None
            self._set_reference(self._x.value.T.copy())
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
class FiniteBurnPlan:
    """ZOH finite-burn plan with the LCvx losslessness audit attached.

    lcvx_gaps[k] = Gamma_k - ||u_k||: zero (to solver tolerance) wherever
    the relaxation is lossless. Discrete-time LCvx is not unconditionally
    lossless; the theory bounds the number of non-tight nodes by the
    state dimension, so `lcvx_inactive` larger than 6 (or a large
    `lcvx_max_gap`) means the relaxed solution is NOT a valid thrust
    profile and must not be flown."""

    times: F64  # (K+1,)
    u: F64  # (K, 3) ZOH accelerations
    gammas: F64  # (K,) slack magnitudes
    states: F64  # (K+1, 6)
    objective: float
    status: str
    lcvx_gaps: F64
    lcvx_inactive: list[int]
    lcvx_max_gap: float
    controllable: bool

    def total_dv(self, dt: float) -> float:
        return float(np.sum(np.linalg.norm(self.u, axis=1)) * dt)


class FiniteBurnPlanner:
    """Min-fuel finite-burn rendezvous with a thrust annulus via lossless
    convexification (Acikmese-Blackmore relaxation).

    The nonconvex lower bound rho_min <= ||u_k|| is relaxed with a slack:
    ||u_k|| <= Gamma_k, rho_min <= Gamma_k <= rho_max, min sum Gamma_k dt.
    In continuous time the relaxation is exact under controllability/
    normality; in discrete time that guarantee weakens, so every solve
    ships a validity audit (see FiniteBurnPlan) and a controllability
    check instead of assuming the theorem. Dynamics: exact ZOH CW
    discretization on a uniform grid (CW-only v0)."""

    def __init__(self, times: F64, annulus: AnnulusSpec) -> None:
        self.times = np.asarray(times, dtype=np.float64)
        self.k = len(self.times) - 1
        dts = np.diff(self.times)
        if not np.allclose(dts, dts[0]):
            raise ValueError("finite-burn grid must be uniform (ZOH)")
        self.dt = float(dts[0])
        self.annulus = annulus

        x = cp.Variable((6, self.k + 1), name="x")
        u = cp.Variable((3, self.k), name="u")
        g = cp.Variable(self.k, name="gamma")
        self._x, self._u, self._g = x, u, g
        self._x0 = cp.Parameter(6, name="x0")
        self._xf = cp.Parameter(6, name="xf")
        self._ad = cp.Parameter((6, 6), name="ad")
        self._bd = cp.Parameter((6, 3), name="bd")

        cons = [x[:, 0] == self._x0, x[:, self.k] == self._xf]
        for i in range(self.k):
            cons.append(
                x[:, i + 1] == self._ad @ x[:, i] + self._bd @ u[:, i]
            )
        for i in range(self.k):
            cons.append(cp.norm(u[:, i], 2) <= g[i])
        cons.append(g >= annulus.rho_min)
        cons.append(g <= annulus.rho_max)
        self._problem = cp.Problem(cp.Minimize(cp.sum(g) * self.dt), cons)
        assert self._problem.is_dcp(dpp=True)

    def solve(self, x0: F64, xf: F64, n: float) -> FiniteBurnPlan:
        ad, bd = _lqr.cw_discrete(n, self.dt)
        # controllability precondition (Kalman rank)
        ctrb = np.hstack([np.linalg.matrix_power(ad, i) @ bd for i in range(6)])
        controllable = int(np.linalg.matrix_rank(ctrb)) == 6
        self._ad.value = ad
        self._bd.value = bd
        self._x0.value = np.asarray(x0, dtype=np.float64)
        self._xf.value = np.asarray(xf, dtype=np.float64)
        self._problem.solve(solver=cp.CLARABEL)
        assert self._u.value is not None and self._x.value is not None
        assert self._g.value is not None
        u = self._u.value.T.copy()
        g = self._g.value.copy()
        gaps = g - np.linalg.norm(u, axis=1)
        tol = 1e-6 * self.annulus.rho_max
        inactive = [int(i) for i in np.flatnonzero(gaps > tol)]
        return FiniteBurnPlan(
            times=self.times.copy(),
            u=u,
            gammas=g,
            states=self._x.value.T.copy(),
            objective=float(self._problem.value),
            status=str(self._problem.status),
            lcvx_gaps=gaps,
            lcvx_inactive=inactive,
            lcvx_max_gap=float(np.max(gaps)),
            controllable=controllable,
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

    Internally the ROE state is scaled by _SCALE: raw ROE magnitudes
    (~1e-5) sit at interior-point absolute tolerances, and the safe-set
    variant returns 'optimal_inaccurate' without it. The scaling is exact
    (linear dynamics); inputs and outputs are unscaled.
    """

    _SCALE = 1.0e5

    def __init__(
        self,
        times: F64,
        objective: str = "l2",
        safe_set: SafeSetSpec | None = None,
    ) -> None:
        self.times = np.asarray(times, dtype=np.float64)
        self.k = len(self.times) - 1
        self.safe_set = safe_set
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
        r_term = r[:, self.k] + self._gf @ v[:, self.k]
        if safe_set is None:
            cons.append(r_term == self._rf)
        else:
            uhat = np.asarray(safe_set.direction, dtype=np.float64)
            uhat = uhat / np.linalg.norm(uhat)
            perp = np.eye(2) - np.outer(uhat, uhat)
            tan_c = math.tan(safe_set.cone_angle)
            de = r_term[2:4]
            di = r_term[4:6]
            cons.append(uhat @ de >= safe_set.e_min * self._SCALE)
            cons.append(uhat @ di >= safe_set.i_min * self._SCALE)
            cons.append(cp.norm(perp @ de, 2) <= tan_c * (uhat @ de))
            cons.append(cp.norm(perp @ di, 2) <= tan_c * (uhat @ di))
            cons.append(cp.abs(r_term[0]) <= safe_set.da_tol * self._SCALE)
        if objective == "l2":
            obj = cp.sum([cp.norm(v[:, i], 2) for i in range(self.k + 1)])
        else:
            obj = cp.sum([cp.norm(v[:, i], 1) for i in range(self.k + 1)])
        self._problem = cp.Problem(cp.Minimize(obj), cons)
        assert self._problem.is_dcp(dpp=True)

    def solve(
        self,
        roe0: F64,
        roef: F64 | None,
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
        """roef is the terminal ROE target (required unless the planner
        was built with a safe_set, in which case it is ignored)."""
        if roef is None and self.safe_set is None:
            raise ValueError("roef required without a safe_set")
        for i in range(self.k):
            dt = float(self.times[i + 1] - self.times[i])
            if j2 > 0.0:
                phi = roe_mod.stm_j2(mu, j2, r_body, a, e, inc, argp, dt)
            else:
                phi = roe_mod.stm_keplerian(n, dt)
            u_i = u0 + n * float(self.times[i])
            gamma = roe_mod.control_matrix(a, n, u_i)
            self._phis[i].value = phi
            self._phigs[i].value = self._SCALE * (phi @ gamma)
        self._gf.value = self._SCALE * roe_mod.control_matrix(
            a, n, u0 + n * float(self.times[self.k])
        )
        self._r0.value = self._SCALE * np.asarray(roe0, dtype=np.float64)
        if self.safe_set is None:
            self._rf.value = self._SCALE * np.asarray(roef, dtype=np.float64)
        self._problem.solve(solver=cp.CLARABEL)
        assert self._v.value is not None and self._r.value is not None
        return RoePlan(
            times=self.times.copy(),
            dvs=self._v.value.T.copy(),
            roes=self._r.value.T.copy() / self._SCALE,
            objective=float(self._problem.value),
            status=str(self._problem.status),
        )
