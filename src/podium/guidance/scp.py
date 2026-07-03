"""PTR successive convexification for impulsive rendezvous (Layer 1).

Solves the NONCONVEX problem Layer-0 approximates: min-fuel impulsive
transfer with a true spherical keep-out constraint ||r(t)|| >= R —
enforced in CONTINUOUS time, not only at nodes.

Method (penalized trust region, with an SCvx*-style penalty ramp):
  - dynamics: exact CW STM equality constraints (zero discretization
    error, as in Layer 0);
  - KOZ linearized about the current reference: rhat_k . r_k >= R - s_k
    with virtual buffers s_k >= 0 penalized by w_pen (feasibility must
    be earned: converged solutions have s == 0);
  - quadratic trust region ||X - X_ref||^2 weighted by w_tr keeps each
    subproblem honest about where the linearization is valid;
  - SCvx*-style update: if the true violation stalls while buffers are
    still active, w_pen ramps by gamma (bounded ramp count).

Continuous-time constraint satisfaction (CTCS), exact-flow form:
between burns the trajectory follows the exact linear flow
p(tau) = S Phi(tau) (x_k + B dv_k), which is LINEAR in the decision
variables. So instead of integral augmentation, we densely sample the
exact flow of each iterate, and every violation becomes an exact linear
cut at its worst time — iterated until a dense check is clean. For
coast-arc RPOD this is CTCS without approximation.

Sandbox layer (cvxpy); the flight-side product is the impulse schedule.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import cvxpy as cp
import numpy as np
from numpy.typing import NDArray

from podium.core import cw

F64 = NDArray[np.float64]

_B = np.zeros((6, 3))
_B[3:6, :] = np.eye(3)
_S = np.zeros((3, 6))
_S[:, 0:3] = np.eye(3)


@dataclass
class ScpResult:
    times: F64
    dvs: F64  # (K+1, 3)
    states: F64  # (K+1, 6) pre-burn
    objective: float  # true fuel cost (sum of burn norms)
    status: str  # "converged" | "max_iters" | solver status on failure
    iterations: int
    slack_final: float  # virtual-buffer norm at exit (0 => feasible)
    dense_violation: float  # worst continuous-time KOZ violation [m]
    n_cuts: int
    w_pen_final: float
    history: list[dict] = field(default_factory=list)

    def total_dv(self) -> float:
        return float(np.sum(np.linalg.norm(self.dvs, axis=1)))


@dataclass(frozen=True)
class EventuallyBoxSpec:
    """STL timed-reach: eventually within [t_lo, t_hi], the position is
    inside the axis-aligned box (center +- half). Robustness
    rho = max_{k in window} min_i (half_i - |r_k,i - c_i|); encoded in
    the SCP via the log-sum-exp smooth max, consumed CONSERVATIVELY:
    LSE/tau >= eps + ln(K)/tau  =>  true rho >= eps."""

    t_lo: float
    t_hi: float
    center: tuple[float, float, float]
    half: tuple[float, float, float]
    eps: float = 1.0  # required true-robustness margin [m]
    # smoothing sharpness [1/m]: conservatism is ln(K)/tau, so tau must
    # keep that well under the box half-width (tau=0.5, K=5 -> 3.2 m)
    tau: float = 0.5


class PtrDockingPlanner:
    """PTR loop for impulsive transfers with a true keep-out sphere."""

    def __init__(
        self,
        times: F64,
        koz_radius: float,
        objective: str = "l2",
        dv_max: float | None = None,
        ctcs_samples: int = 25,
        max_iters: int = 20,
        w_tr: float = 1e-4,
        w_pen0: float = 1.0,
        pen_gamma: float = 8.0,
        tol_feas: float = 1e-4,
        tol_step: float = 1e-3,
        stl_reach: EventuallyBoxSpec | None = None,
    ) -> None:
        self.times = np.asarray(times, dtype=np.float64)
        self.k = len(self.times) - 1
        self.r_koz = koz_radius
        self.objective = objective
        self.dv_max = dv_max
        self.ctcs_samples = ctcs_samples
        self.max_iters = max_iters
        self.w_tr = w_tr
        self.w_pen0 = w_pen0
        self.pen_gamma = pen_gamma
        self.tol_feas = tol_feas
        self.tol_step = tol_step
        self.stl_reach = stl_reach
        self.stl_nodes: list[int] = []
        if stl_reach is not None:
            self.stl_nodes = [
                i for i, t in enumerate(self.times)
                if stl_reach.t_lo <= t <= stl_reach.t_hi and 0 < i <= self.k
            ]
            if not self.stl_nodes:
                raise ValueError("STL window contains no trajectory nodes")

    # -- STL smooth robustness --------------------------------------------
    def _stl_margins(self, states: F64) -> F64:
        """Per-window-node box margins m_k = min_i(half_i - |r_i - c_i|)."""
        sp = self.stl_reach
        assert sp is not None
        out = np.empty(len(self.stl_nodes))
        for j, k in enumerate(self.stl_nodes):
            d = states[k, 0:3] - np.asarray(sp.center)
            out[j] = float(np.min(np.asarray(sp.half) - np.abs(d)))
        return out

    def stl_true_robustness(self, states: F64) -> float:
        """Exact (non-smooth) STL robustness of the node trajectory."""
        return float(np.max(self._stl_margins(states)))

    def _stl_lse(self, states: F64) -> tuple[float, F64, F64]:
        """(LSE/tau smoothed max, softmax weights, reference margins).

        Soundness structure: node margins enter the subproblem EXACTLY
        via hypograph variables (m is concave: m <= each affine face
        expression), and only the LSE — convex in the margins — is
        tangent-linearized. A convex function dominates its tangent, so
        tangent >= target implies LSE >= target implies true max-margin
        >= eps. Linearizing the concave margins directly instead is a
        RELAXATION the optimizer exploits by exiting through a face the
        subgradient didn't pick (observed as a period-2 oscillation)."""
        sp = self.stl_reach
        assert sp is not None
        m = self._stl_margins(states)
        mmax = float(np.max(m))
        w = np.exp(sp.tau * (m - mmax))
        lse = mmax + math.log(float(np.sum(w))) / sp.tau
        w = w / float(np.sum(w))
        return lse, w, m

    # -- exact flow helpers ---------------------------------------------
    def _flow_position(self, n: float, xk: F64, dvk: F64, tau: float) -> F64:
        out: F64 = _S @ (cw.stm(n, tau) @ (xk + _B @ dvk))
        return out

    def _dense_worst(self, n: float, states: F64, dvs: F64) -> tuple[float, list]:
        """Worst continuous-time KOZ violation over all coast arcs.

        Returns (worst_violation_m, cut_list) where each cut is
        (node k, tau, unit direction at the violation point)."""
        worst = 0.0
        cuts = []
        for k in range(self.k):
            dt = self.times[k + 1] - self.times[k]
            arc_worst, arc_tau = 0.0, 0.0
            for i in range(1, self.ctcs_samples):
                tau = dt * i / self.ctcs_samples
                p = self._flow_position(n, states[k], dvs[k], tau)
                dist = float(np.linalg.norm(p))
                v = self.r_koz - dist
                if v > arc_worst and dist > 1e-9:
                    arc_worst, arc_tau = v, tau
            if arc_worst > 0.0:
                worst = max(worst, arc_worst)
                cuts.append((k, arc_tau))
        return worst, cuts

    def _node_worst(self, states: F64) -> float:
        worst = 0.0
        for k in range(1, self.k):
            worst = max(worst, self.r_koz - float(np.linalg.norm(states[k, 0:3])))
        return worst

    # -- one convex subproblem -------------------------------------------
    def _subproblem(self, n: float, x0: F64, xf: F64, ref: F64,
                    ref_dv: F64, cuts: list, w_pen: float,
                    w_tr: float | None = None) -> "tuple | None":
        x = cp.Variable((6, self.k + 1))
        v = cp.Variable((3, self.k + 1))
        n_node = self.k - 1
        n_stl = 1 if self.stl_reach is not None else 0
        s = cp.Variable(n_node + len(cuts) + n_stl, nonneg=True) \
            if (n_node + len(cuts) + n_stl) else None
        cons = [x[:, 0] == x0]
        phi = cw.stm(n, float(self.times[1] - self.times[0]))
        for i in range(self.k):
            dt = float(self.times[i + 1] - self.times[i])
            phi_i = phi if abs(dt - (self.times[1] - self.times[0])) < 1e-12 \
                else cw.stm(n, dt)
            cons.append(x[:, i + 1] == phi_i @ x[:, i] + (phi_i @ _B) @ v[:, i])
        cons.append(x[:, self.k] + _B @ v[:, self.k] == xf)
        if self.dv_max is not None:
            for i in range(self.k + 1):
                cons.append(cp.norm(v[:, i], 2) <= self.dv_max)
        # node KOZ, linearized about the reference
        for j, k in enumerate(range(1, self.k)):
            r_ref = ref[k, 0:3]
            mag = float(np.linalg.norm(r_ref))
            rhat = r_ref / mag if mag > 1e-9 else np.array([1.0, 0.0, 0.0])
            assert s is not None
            cons.append(rhat @ x[0:3, k] >= self.r_koz - s[j])
        # continuous-time cuts: times are persistent, directions are
        # RE-LINEARIZED at the current reference every iteration (stale
        # fixed-direction cuts can permanently exclude the optimum and
        # keep virtual buffers active forever)
        for j, (k, tau) in enumerate(cuts):
            p_ref = self._flow_position(n, ref[k], ref_dv[k], float(tau))
            mag = float(np.linalg.norm(p_ref))
            rhat = p_ref / mag if mag > 1e-9 else np.array([1.0, 0.0, 0.0])
            phi_tau = cw.stm(n, float(tau))
            row = rhat @ (_S @ phi_tau)
            assert s is not None
            cons.append(row @ (x[:, k] + _B @ v[:, k])
                        >= self.r_koz - s[n_node + j])
        # STL timed-reach: exact hypograph margins + LSE tangent (see
        # _stl_lse for why this split is the sound one)
        if self.stl_reach is not None:
            sp = self.stl_reach
            lse, w, m_ref = self._stl_lse(ref)
            target = sp.eps + math.log(len(self.stl_nodes)) / sp.tau
            mk = cp.Variable(len(self.stl_nodes))
            for j, k in enumerate(self.stl_nodes):
                for i in range(3):
                    d_i = x[i, k] - sp.center[i]
                    cons.append(mk[j] <= sp.half[i] - d_i)
                    cons.append(mk[j] <= sp.half[i] + d_i)
            expr = lse + cp.sum(cp.multiply(w, mk - m_ref))
            assert s is not None
            cons.append(expr >= target - s[n_node + len(cuts)])

        fuel = cp.sum([cp.norm(v[:, i], 1 if self.objective == "l1" else 2)
                       for i in range(self.k + 1)])
        wt = self.w_tr if w_tr is None else w_tr
        obj = fuel + wt * (cp.sum_squares(x - ref.T)
                           + cp.sum_squares(v - ref_dv.T))
        if s is not None:
            obj = obj + w_pen * cp.sum(s)
        prob = cp.Problem(cp.Minimize(obj), cons)
        prob.solve(solver=cp.CLARABEL)
        if x.value is None:
            return None
        slack = float(np.sum(s.value)) if s is not None and s.value is not None \
            else 0.0
        return x.value.T.copy(), v.value.T.copy(), slack, str(prob.status)

    # -- the PTR loop ------------------------------------------------------
    def solve(self, x0: F64, xf: F64, n: float) -> ScpResult:
        ref = np.zeros((self.k + 1, 6))
        for i in range(self.k + 1):
            f = i / self.k
            ref[i, 0:3] = (1 - f) * x0[0:3] + f * xf[0:3]
        ref_dv = np.zeros((self.k + 1, 3))
        cuts: list = []
        w_pen = self.w_pen0
        w_tr = self.w_tr
        history: list[dict] = []
        prev_viol = math.inf
        prev_fuel = math.inf
        status = "max_iters"
        states, dvs, slack = ref, ref_dv, math.inf
        it = 0
        for it in range(1, self.max_iters + 1):
            out = self._subproblem(n, x0, xf, ref, ref_dv, cuts, w_pen, w_tr)
            if out is None:
                return ScpResult(self.times.copy(), ref_dv, ref, math.inf,
                                 "subproblem_failed", it, math.inf, math.inf,
                                 len(cuts), w_pen, history)
            states, dvs, slack, _ = out
            step = float(np.max(np.abs(states - ref)))
            fuel = float(np.sum(np.linalg.norm(dvs, axis=1)))
            node_viol = self._node_worst(states)
            dense_viol, new_cuts = self._dense_worst(n, states, dvs)
            stl_short = 0.0
            if self.stl_reach is not None:
                stl_short = max(0.0, self.stl_reach.eps
                                - self.stl_true_robustness(states))
            true_viol = max(node_viol, dense_viol, stl_short)
            history.append({"iter": it, "step": step, "slack": slack,
                            "violation": true_viol, "w_pen": w_pen,
                            "w_tr": w_tr, "cuts": len(cuts), "fuel": fuel})
            # SCvx*-style updates:
            #  - infeasibility stalling with buffers active -> penalty ramp;
            #  - feasible iterate whose linearization VALIDATED (true
            #    violation clean) -> expand the trust region so nearly-flat
            #    fuel valleys don't force a 1-(mu/w_tr) crawl.
            if slack > self.tol_feas and true_viol >= 0.9 * prev_viol:
                w_pen *= self.pen_gamma
            elif slack <= self.tol_feas and true_viol <= self.tol_feas:
                # bounded expansion: enough to break flat-valley crawl,
                # floored so linearized constraints keep a restoring pull
                # (w_tr -> 0 makes the loop oscillate, observed)
                w_tr = max(self.w_tr * 0.04, w_tr * 0.2)
            prev_viol = true_viol
            ref, ref_dv = states, dvs
            if new_cuts and dense_viol > self.tol_feas:
                # register at most one cut per arc per iteration
                known = {(k, round(t, 6)) for k, t in cuts}
                for c in new_cuts:
                    if (c[0], round(c[1], 6)) not in known:
                        cuts.append(c)
                prev_fuel = fuel
                continue  # constraint set changed: keep iterating
            # Converged when feasible (buffers empty, continuous check
            # clean) and stationary. Min-fuel problems have flat valleys
            # (many equally-cheap burn splits) along which the iterate can
            # drift without the cost changing — fuel stationarity is the
            # honest criterion there, small step the strict one.
            fuel_stat = abs(fuel - prev_fuel) <= 1e-5 * max(0.1, fuel)
            prev_fuel = fuel
            if slack <= self.tol_feas and true_viol <= self.tol_feas \
                    and (step <= self.tol_step or fuel_stat):
                status = "converged"
                break
        fuel = float(np.sum(np.linalg.norm(dvs, axis=1)))
        dense_viol, _ = self._dense_worst(n, states, dvs)
        return ScpResult(self.times.copy(), dvs, states, fuel, status, it,
                         slack, dense_viol, len(cuts), w_pen, history)
