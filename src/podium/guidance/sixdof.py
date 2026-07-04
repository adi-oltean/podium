"""6-DOF attitude-coupled PTR: joint translation + attitude planning.

The coupling is physical, not decorative: the chaser has ONE body-fixed
thruster along +x_body, so the translational acceleration is
(T/m) R(q) e1 — you cannot brake without slewing, and the planner must
discover the slew. Torque is 3-axis (wheels/RCS pairs), bounded.

State (13): r(3) LVLH, v(3), q(4) body->LVLH, w(3) body rates.
Control (4): throttle T in [0, t_max] (thrust magnitude, N) and
torque tau(3), |tau_i| <= tau_max.

Method — PTR exactly as in guidance/scp.py, adapted to the coupled
dynamics: propagate the nonlinear reference with RK4 (quaternion
renormalized at node boundaries), build discrete Jacobians by central
finite differences (17 states + 4 controls at ~10 nodes is cheap),
solve the convex subproblem with virtual-control slack (penalty ramps
when stalled) and a floored trust region, iterate until defects and
slack vanish. The quaternion norm is not constrained convexly; the
reference renormalization keeps drift ~1e-9 per iteration (measured in
the receipts via nonlinear replay).

Frame note: CW translational dynamics (circular target), so this is a
terminal-approach planner (tens of meters), where CW is mm-accurate
over the horizon.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from podium.core import quat

F64 = NDArray[np.float64]

NX = 13
NU = 4


def _deriv(x: F64, u: F64, n: float, mass: float, inertia: F64,
           inertia_inv: F64) -> F64:
    r, v, q, w = x[0:3], x[3:6], x[6:10], x[10:13]
    thrust, tau = u[0], u[1:4]
    a_cw = np.array([
        3.0 * n * n * r[0] + 2.0 * n * v[1],
        -2.0 * n * v[0],
        -n * n * r[2],
    ])
    a_thr = (thrust / mass) * quat.rotate(q, np.array([1.0, 0.0, 0.0]))
    out = np.empty(NX)
    out[0:3] = v
    out[3:6] = a_cw + a_thr
    out[6:10] = quat.deriv(q, w)
    out[10:13] = inertia_inv @ (tau - np.cross(w, inertia @ w))
    return out


@dataclass
class SixDofPlan:
    status: str
    times: F64
    states: F64  # (K, 13) nonlinear-propagated reference at convergence
    controls: F64  # (K-1, 4) zero-order-hold
    iterations: int
    slack: float
    defect: float
    fuel: float  # integral of |T| dt  [N s]
    history: list = field(default_factory=list)


class SixDofPlanner:
    """PTR for the body-fixed-thruster terminal approach."""

    def __init__(
        self,
        times: F64,
        mass: float = 500.0,
        inertia: F64 | None = None,
        t_max: float = 10.0,
        tau_max: float = 0.2,
        n_substeps: int = 4,
    ) -> None:
        self.times = np.asarray(times, dtype=np.float64)
        self.mass = mass
        self.inertia = (np.diag([120.0, 90.0, 60.0])
                        if inertia is None else inertia)
        self.inertia_inv = np.linalg.inv(self.inertia)
        self.t_max = t_max
        self.tau_max = tau_max
        self.n_substeps = n_substeps
        self.w_pen = 1e3
        self.w_tr = 1.0

    # -- nonlinear propagation over one interval -----------------------
    def _step(self, x: F64, u: F64, dt: float, n: float) -> F64:
        h = dt / self.n_substeps
        y = x.copy()
        for _ in range(self.n_substeps):
            k1 = _deriv(y, u, n, self.mass, self.inertia, self.inertia_inv)
            k2 = _deriv(y + 0.5 * h * k1, u, n, self.mass, self.inertia,
                        self.inertia_inv)
            k3 = _deriv(y + 0.5 * h * k2, u, n, self.mass, self.inertia,
                        self.inertia_inv)
            k4 = _deriv(y + h * k3, u, n, self.mass, self.inertia,
                        self.inertia_inv)
            y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        y[6:10] = quat.normalize(y[6:10])
        return y

    def _jacobians(self, x: F64, u: F64, dt: float, n: float
                   ) -> tuple[F64, F64, F64]:
        """Central-difference discrete Jacobians and the defect base."""
        x1 = self._step(x, u, dt, n)
        a = np.zeros((NX, NX))
        b = np.zeros((NX, NU))
        eps_x = 1e-6
        eps_u = 1e-6
        for i in range(NX):
            dxp = x.copy()
            dxm = x.copy()
            dxp[i] += eps_x
            dxm[i] -= eps_x
            a[:, i] = (self._step(dxp, u, dt, n)
                       - self._step(dxm, u, dt, n)) / (2.0 * eps_x)
        for i in range(NU):
            dup = u.copy()
            dum = u.copy()
            dup[i] += eps_u
            dum[i] -= eps_u
            b[:, i] = (self._step(x, dup, dt, n)
                       - self._step(x, dum, dt, n)) / (2.0 * eps_u)
        return a, b, x1

    # -- initial reference ----------------------------------------------
    def _initial_reference(self, x0: F64, xf: F64) -> tuple[F64, F64]:
        """Straight-line translation + SLERP attitude + a small thrust
        seed. The seeds are load-bearing: at a zero-thrust constant-
        attitude reference the thrust-attitude bilinearity d(T R(q)e1)
        vanishes (dT enters only along the fixed reference axis, dq
        enters times T=0), the linear model cannot reach any terminal
        needing off-axis force, and PTR wedges at irreducible slack —
        observed before this fix, pinned by the receipts."""
        k = len(self.times)
        xs = np.zeros((k, NX))
        us = np.zeros((k - 1, NU))
        q0, qf = x0[6:10], xf[6:10]
        dot = float(np.dot(q0, qf))
        qf_s = qf if dot >= 0.0 else -qf
        ang = float(np.arccos(np.clip(abs(dot), -1.0, 1.0)))
        for i in range(k):
            s = i / (k - 1)
            xs[i, 0:3] = (1 - s) * x0[0:3] + s * xf[0:3]
            xs[i, 3:6] = (1 - s) * x0[3:6] + s * xf[3:6]
            if ang < 1e-9:
                qi = q0.copy()
            else:
                qi = (math.sin((1 - s) * ang) * q0
                      + math.sin(s * ang) * qf_s) / math.sin(ang)
            xs[i, 6:10] = quat.normalize(qi)
            xs[i, 10:13] = 0.0
        us[:, 0] = 0.2 * self.t_max  # thrust seed keeps dq coupling alive
        return xs, us

    def solve(
        self,
        x0: F64,
        xf: F64,
        n: float,
        max_iter: int = 40,
        tol_defect: float = 1e-4,
        tol_slack: float = 1e-6,
    ) -> SixDofPlan:
        import cvxpy as cp

        k = len(self.times)
        dts = np.diff(self.times)
        xs, us = self._initial_reference(x0, xf)
        w_pen, w_tr = self.w_pen, self.w_tr
        history = []
        status = "max_iter"
        best = None
        stall = 0
        # hard trust-region boxes per state block (position, velocity,
        # quaternion, rate) and per control — soft quadratic weights
        # alone let the subproblem take steps that blow up the
        # nonlinear attitude propagation (observed: NaN reference)
        tr0 = {"pos": 8.0, "vel": 0.3, "q": 0.4, "w": 0.03,
               "T": 0.6 * self.t_max, "tau": self.tau_max}
        tr = dict(tr0)
        prev_defect = float("inf")

        for it in range(max_iter):
            a_list, b_list, base = [], [], []
            for i in range(k - 1):
                a, b, x1 = self._jacobians(xs[i], us[i], dts[i], n)
                a_list.append(a)
                b_list.append(b)
                base.append(x1)

            dx = cp.Variable((k, NX))
            du = cp.Variable((k - 1, NU))
            nu = cp.Variable((k - 1, NX))
            cons = [dx[0] == 0]
            for i in range(k - 1):
                cons.append(
                    xs[i + 1] + dx[i + 1] ==
                    base[i] + a_list[i] @ dx[i] + b_list[i] @ du[i]
                    + nu[i])
            tvar = us[:, 0] + du[:, 0]
            cons += [tvar >= 0, tvar <= self.t_max]
            for j in range(1, 4):
                cons += [cp.abs(us[:, j] + du[:, j]) <= self.tau_max]
            # hard trust region
            cons += [cp.abs(dx[:, 0:3]) <= tr["pos"],
                     cp.abs(dx[:, 3:6]) <= tr["vel"],
                     cp.abs(dx[:, 6:10]) <= tr["q"],
                     cp.abs(dx[:, 10:13]) <= tr["w"],
                     cp.abs(du[:, 0]) <= tr["T"],
                     cp.abs(du[:, 1:4]) <= tr["tau"]]
            # terminal boxes (position/velocity exact, attitude/rate box)
            cons.append(xs[-1, 0:6] + dx[-1, 0:6] == xf[0:6])
            cons.append(cp.abs(xs[-1, 6:10] + dx[-1, 6:10] - xf[6:10])
                        <= 0.02)
            cons.append(cp.abs(xs[-1, 10:13] + dx[-1, 10:13]) <= 0.002)

            fuel = cp.sum(cp.multiply(dts, tvar))
            obj = (fuel
                   + w_pen * cp.sum(cp.abs(nu))
                   + w_tr * (cp.sum_squares(dx) + cp.sum_squares(du)))
            prob = cp.Problem(cp.Minimize(obj), cons)
            try:
                prob.solve(solver=cp.CLARABEL)
            except cp.error.SolverError:
                prob.solve(solver=cp.SCS, eps=1e-8, max_iters=200_000)
            if prob.status not in ("optimal", "optimal_inaccurate"):
                status = f"subproblem_{prob.status}"
                break

            xs_new = xs + dx.value
            us_new = us + du.value
            for i in range(k):
                xs_new[i, 6:10] = quat.normalize(xs_new[i, 6:10])

            # validate on the NONLINEAR dynamics
            defect = 0.0
            y = xs_new[0].copy()
            traj = [y.copy()]
            for i in range(k - 1):
                y = self._step(y, us_new[i], dts[i], n)
                defect = max(defect,
                             float(np.max(np.abs(y - xs_new[i + 1]))))
                traj.append(y.copy())
            slack = float(np.sum(np.abs(nu.value)))
            step = float(np.max(np.abs(dx.value)))
            history.append((it, slack, defect, step,
                            float(np.sum(dts * us_new[:, 0]))))

            # accept/reject: a step whose nonlinear replay explodes or
            # sharply regresses shrinks the trust region and retries
            if not np.all(np.isfinite(defect)) or not np.isfinite(
                    float(np.max(np.abs(xs_new)))) \
                    or defect > max(10.0 * prev_defect, 1.0):
                for key in tr:
                    tr[key] *= 0.5
                continue
            # once feasible, CONTRACT the trust region whenever defect
            # stops improving — the soft-weight-decay rule from scp.py
            # is backwards here (bigger steps near the solution) and
            # produced a period-2 limit cycle at defect 0.465, measured
            if slack <= tol_slack and defect > 0.7 * prev_defect:
                for key in tr:
                    tr[key] = max(tr[key] * 0.5, tr0[key] * 1e-4)
            elif defect <= 0.7 * prev_defect:
                for key in tr:
                    tr[key] = min(tr0[key], tr[key] * 1.3)
            prev_defect = min(prev_defect, defect)
            xs, us = xs_new, us_new
            if slack <= tol_slack and defect <= tol_defect:
                status = "converged"
                best = (np.array(traj), us.copy())
                break
            if slack <= tol_slack and step < 1e-5:
                stall += 1
                if stall >= 2:
                    status = "stalled_feasible" \
                        if defect <= 10 * tol_defect else "stalled"
                    best = (np.array(traj), us.copy())
                    break
            else:
                stall = 0
            # ramp the slack penalty while infeasible
            if slack > tol_slack:
                w_pen = min(w_pen * 4.0, 1e6)

        if best is None:
            y = xs[0].copy()
            traj = [y.copy()]
            for i in range(k - 1):
                y = self._step(y, us[i], dts[i], n)
                traj.append(y.copy())
            best = (np.array(traj), us.copy())
        states, controls = best
        return SixDofPlan(
            status=status,
            times=self.times.copy(),
            states=states,
            controls=controls,
            iterations=len(history),
            slack=history[-1][1] if history else float("nan"),
            defect=float(np.max(np.abs(states[-1, 0:6] - xf[0:6]))),
            fuel=float(np.sum(dts * controls[:, 0])),
            history=history,
        )
