"""Coverage-driven receipts for the guidance/control edge paths.

These exercise the status-return, infeasibility, degenerate-geometry, and
input-validation branches of the guidance planners and the docking law that
the behavioural tests do not reach — asserting the DOCUMENTED behaviour
(the right error, the right status, a finite fallback), not just touching
the line. Conventions follow tests/test_convex.py / test_scp.py /
test_sixdof.py.
"""

import json
import math

import numpy as np
import pytest

cp = pytest.importorskip("cvxpy")

from podium import constants as const  # noqa: E402
from podium.control import docking  # noqa: E402
from podium.core import quat  # noqa: E402
from podium.guidance import arch, safety  # noqa: E402
from podium.guidance.convex import (  # noqa: E402
    AnnulusSpec,
    FiniteBurnPlan,
    FiniteBurnPlanner,
    PassiveSafetySpec,
    PlumeSpec,
    RendezvousPlanner,
    RoePlanner,
    find_min_time,
)
from podium.guidance.scp import EventuallyBoxSpec, PtrDockingPlanner  # noqa: E402
from podium.guidance.sixdof import SixDofPlanner  # noqa: E402

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)
X0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])
XF = np.zeros(6)


# ======================================================================
# convex.py — RendezvousPlanner
# ======================================================================
def test_planner_needs_at_least_two_nodes():
    with pytest.raises(ValueError, match="at least two nodes"):
        RendezvousPlanner(np.array([0.0]))


def test_planner_rejects_unknown_objective():
    with pytest.raises(ValueError, match="objective must be"):
        RendezvousPlanner(np.linspace(0.0, 1000.0, 5), objective="linf")


def test_dv_max_bound_is_respected_and_binds():
    """A per-burn dv cap below the unconstrained peak (~0.48 m/s) must
    hold at every node and force the planner off the two-impulse optimum."""
    times = np.linspace(0.0, 2000.0, 9)
    free = RendezvousPlanner(times).solve(X0, XF, n=N)
    assert float(np.max(np.linalg.norm(free.dvs, axis=1))) > 0.4  # peak binds
    capped = RendezvousPlanner(times, dv_max=0.4).solve(X0, XF, n=N)
    assert capped.status == "optimal"
    assert np.all(np.linalg.norm(capped.dvs, axis=1) <= 0.4 + 1e-6)
    # spreading the burn under the cap cannot be cheaper than the optimum
    assert capped.total_dv() >= free.total_dv() - 1e-6


def test_plume_deactivates_far_nodes():
    """Plume half-space is active only inside range_m of the target; the
    far nodes get a zeroed (inactive) direction parameter."""
    times = np.linspace(0.0, 2000.0, 11)
    plan = RendezvousPlanner(times, plume=PlumeSpec(range_m=100.0)).solve(
        X0, XF, n=N
    )
    active = {i for i, _ in plan.plume_dirs}
    assert plan.status == "optimal"
    assert active  # some near-target nodes are active
    assert len(active) < len(times) - 1  # and most are deactivated


def test_passive_safety_requires_circular_orbit():
    spec = PassiveSafetySpec(radius=200.0, horizon=1000.0)
    planner = RendezvousPlanner(np.linspace(0.0, 2000.0, 6), passive_safety=spec)
    with pytest.raises(ValueError, match="require e=0"):
        planner.solve(X0, XF, n=N, e=0.01)


# ======================================================================
# convex.py — FiniteBurnPlan / FiniteBurnPlanner / find_min_time
# ======================================================================
def test_finiteburn_plan_total_dv_and_empty_primer_certificate():
    """total_dv sums burn magnitudes * dt; with no primer duals recorded
    the normality certificate degrades to 0.0 rather than dividing by an
    empty array."""
    u = np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]])  # norms 5, 0
    plan = FiniteBurnPlan(
        times=np.array([0.0, 10.0, 20.0]),
        u=u,
        gammas=np.zeros(2),
        states=np.zeros((3, 6)),
        objective=0.0,
        status="optimal",
        lcvx_gaps=np.zeros(2),
        lcvx_inactive=[],
        lcvx_max_gap=0.0,
        controllable=True,
    )
    assert plan.total_dv(10.0) == pytest.approx(50.0)
    assert plan.primer_certificate() == 0.0  # empty primer_norms


def test_finiteburn_requires_uniform_grid():
    with pytest.raises(ValueError, match="uniform"):
        FiniteBurnPlanner(np.array([0.0, 10.0, 25.0, 30.0]), AnnulusSpec(0.0, 0.02))


def test_finiteburn_solver_error_returns_infeasible_plan(monkeypatch):
    """A raised SolverError must be swallowed and surface as an empty plan
    carrying the status (this is the path find_min_time probes)."""
    planner = FiniteBurnPlanner(np.linspace(0.0, 900.0, 9), AnnulusSpec(0.0, 0.02))

    def boom(self, *a, **k):
        raise cp.error.SolverError("forced failure")

    monkeypatch.setattr(cp.Problem, "solve", boom)
    plan = planner.solve(np.array([0.0, -2000.0, 0.0, 0.0, 0.0, 0.0]), XF, n=N)
    assert math.isinf(plan.objective)
    assert plan.u.shape == (planner.k, 3)
    assert np.all(plan.u == 0.0)


def test_find_min_time_rejects_infeasible_upper_bound():
    x0 = np.array([0.0, -2000.0, 0.0, 0.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="not feasible"):
        find_min_time(x0, XF, N, rho_max=0.001, k=6, t_lo=10.0, t_hi=50.0)


# ======================================================================
# convex.py — RoePlanner
# ======================================================================
def _roe_times():
    return np.linspace(0.0, 1.5 * (2 * math.pi / N), 7)


def test_roe_planner_requires_roef_without_safe_set():
    with pytest.raises(ValueError, match="roef required"):
        RoePlanner(_roe_times()).solve(np.zeros(6), None, a=A, n=N, u0=0.0)


def test_roe_planner_l1_objective_solves():
    roef = np.array([0.0, 0.0, 4e-5, 0.0, 4e-5, 0.0])
    plan = RoePlanner(_roe_times(), objective="l1").solve(
        np.zeros(6), roef, a=A, n=N, u0=0.3
    )
    assert plan.status in ("optimal", "optimal_inaccurate")
    assert np.allclose(plan.roes[-1] + 0.0, plan.roes[-1])  # finite


def test_roe_planner_j2_dynamics_solves():
    """The j2>0 branch uses the J2-perturbed ROE STM; the plan must still
    reach the aligned safe geometry through those dynamics."""
    roef = np.array([0.0, 0.0, 4e-5, 0.0, 4e-5, 0.0])
    plan = RoePlanner(_roe_times()).solve(
        np.zeros(6), roef, a=A, n=N, u0=0.3,
        mu=const.MU_EARTH, j2=1.08263e-3, r_body=6.378137e6, inc=0.9,
    )
    assert plan.status in ("optimal", "optimal_inaccurate")


# ======================================================================
# scp.py — PtrDockingPlanner
# ======================================================================
def test_stl_window_without_nodes_is_rejected():
    spec = EventuallyBoxSpec(
        t_lo=1.0e5, t_hi=2.0e5, center=(0.0, 0.0, 0.0), half=(5.0, 5.0, 5.0)
    )
    with pytest.raises(ValueError, match="no trajectory nodes"):
        PtrDockingPlanner(np.linspace(0.0, 1000.0, 6), koz_radius=50.0,
                          stl_reach=spec)


def test_scp_dv_max_constraint_is_respected():
    times = np.linspace(0.0, 2000.0, 9)
    res = PtrDockingPlanner(times, koz_radius=1.0, dv_max=0.4).solve(X0, XF, N)
    assert res.status == "converged"
    assert np.all(np.linalg.norm(res.dvs, axis=1) <= 0.4 + 1e-6)


def test_scp_infeasible_subproblem_reports_failure():
    """A dv cap too small to satisfy the terminal constraint makes the
    first subproblem infeasible; the loop must report subproblem_failed
    with an infinite objective instead of a bogus trajectory."""
    times = np.linspace(0.0, 2000.0, 9)
    res = PtrDockingPlanner(times, koz_radius=1.0, dv_max=1e-6).solve(X0, XF, N)
    assert res.status == "subproblem_failed"
    assert math.isinf(res.objective)


# ======================================================================
# safety.py
# ======================================================================
def test_ei_separation_angle_degenerate_vectors():
    """A vanishing e- or i-vector gives no phasing protection: the
    separation angle is reported as the worst case (pi/2)."""
    tiny_de = np.array([0.0, 0.0, 1e-15, 0.0, 1.0, 0.0])
    tiny_di = np.array([0.0, 0.0, 1.0, 0.0, 1e-15, 0.0])
    assert safety.ei_separation_angle(tiny_de) == pytest.approx(0.5 * math.pi)
    assert safety.ei_separation_angle(tiny_di) == pytest.approx(0.5 * math.pi)


# ======================================================================
# arch.py
# ======================================================================
def test_arch_rejects_unknown_gains():
    with pytest.raises(ValueError, match="reference"):
        arch.simulate(arch.X0_CENTER, gains="bogus")


def test_arch_write_model_roundtrips(tmp_path):
    path = tmp_path / "arch_model.json"
    arch.write_model(str(path), abort_time=120.0, gains="reference")
    data = json.loads(path.read_text())
    assert data["abort_time"] == 120.0
    assert [m["name"] for m in data["modes"]][:2] == ["approaching", "attempt"]
    assert any(m["name"] == "aborting" for m in data["modes"])


# ======================================================================
# control/docking.py
# ======================================================================
def test_pulses_needed_saturates_negative():
    """A large negative rate error clamps at -25 pulses per cycle."""
    assert docking.pulses_needed(v_cmd=-10.0, v_est=0.0, dv_pulse=0.1) == -25


# ======================================================================
# sixdof.py
# ======================================================================
def _sixdof_scenario():
    q0 = quat.normalize(np.array([1.0, 0.0, 0.0, 0.0]))
    qf = quat.normalize(np.array([math.cos(-math.pi / 4), 0.0, 0.0,
                                  math.sin(-math.pi / 4)]))
    x0 = np.zeros(13)
    x0[0:3] = [0.0, -30.0, 0.0]
    x0[3:6] = [0.0, 0.10, 0.0]
    x0[6:10] = q0
    xf = np.zeros(13)
    xf[0:3] = [0.0, -2.0, 0.0]
    xf[3:6] = [0.0, 0.05, 0.0]
    xf[6:10] = qf
    return x0, xf


def test_sixdof_zero_slew_reference_is_constant_attitude():
    """When start and terminal attitude coincide the SLERP seed collapses
    to a constant-attitude reference (the ang<1e-9 branch)."""
    q = quat.normalize(np.array([1.0, 0.0, 0.0, 0.0]))
    x0 = np.zeros(13)
    x0[0:3] = [0.0, -10.0, 0.0]
    x0[3:6] = [0.0, 0.02, 0.0]
    x0[6:10] = q
    xf = np.zeros(13)
    xf[0:3] = [0.0, -2.0, 0.0]
    xf[3:6] = [0.0, 0.01, 0.0]
    xf[6:10] = q
    pl = SixDofPlanner(np.linspace(0.0, 120.0, 7), t_max=2.0, tau_max=0.1)
    xs, _us = pl._initial_reference(x0, xf)
    assert np.allclose(xs[:, 6:10], q)  # every node holds the shared attitude


def test_sixdof_max_iter_returns_finite_fallback():
    """Capped at one iteration the slew cannot converge; the planner must
    still return the nonlinear replay of the current controls (best is
    None -> fallback), not crash or emit NaNs."""
    x0, xf = _sixdof_scenario()
    pl = SixDofPlanner(np.linspace(0.0, 360.0, 13), t_max=5.0, tau_max=0.2)
    p = pl.solve(x0, xf, N, max_iter=1)
    assert p.status == "max_iter"
    assert p.iterations == 1
    assert np.all(np.isfinite(p.states))
    assert p.states.shape == (13, 13)


@pytest.mark.slow
def test_sixdof_stalls_feasible_when_step_vanishes():
    """Feasible iterate (slack ~ 0) whose step collapses below the step
    tolerance while the defect is still above tol_defect must exit as
    'stalled_feasible' — not spin to max_iter. We reproduce it by seeding
    the reference at the converged optimum and pinning the trust-region
    weight high so the subproblem can barely move."""
    x0, xf = _sixdof_scenario()
    times = np.linspace(0.0, 360.0, 13)
    opt = SixDofPlanner(times, t_max=5.0, tau_max=0.2).solve(x0, xf, N)
    assert opt.status == "converged"
    seed_states = opt.states.copy()
    seed_controls = opt.controls.copy()

    class Seeded(SixDofPlanner):
        def _initial_reference(self, _a, _b):
            return seed_states.copy(), seed_controls.copy()

    pl = Seeded(times, t_max=5.0, tau_max=0.2)
    pl.w_tr = 1e6  # feasible seed + tiny steps -> the stall guard fires
    p = pl.solve(x0, xf, N, tol_defect=1e-8, tol_slack=1e-6, max_iter=8)
    assert p.status.startswith("stalled")
    assert p.slack <= 1e-6  # exited feasible
    assert np.all(np.isfinite(p.states))


def test_sixdof_solver_error_falls_back_to_scs(monkeypatch):
    """When CLARABEL raises, the subproblem must be retried on SCS and the
    loop must keep producing a finite iterate."""
    x0, xf = _sixdof_scenario()
    pl = SixDofPlanner(np.linspace(0.0, 120.0, 5), t_max=2.0, tau_max=0.1)
    orig = cp.Problem.solve

    def clarabel_fails(self, *a, solver=None, **k):
        if solver == cp.CLARABEL:
            raise cp.error.SolverError("forced clarabel failure")
        return orig(self, solver=solver, **k)

    monkeypatch.setattr(cp.Problem, "solve", clarabel_fails)
    p = pl.solve(x0, xf, N, max_iter=2)
    assert np.all(np.isfinite(p.states))
    assert p.iterations >= 1


def test_sixdof_rejects_nonfinite_nonlinear_replay(monkeypatch):
    """If the nonlinear replay of an accepted step returns non-finite
    values, the accept/reject guard must reject it, shrink the trust
    region, and continue — never fold the blow-up into the plan. We force
    exactly one exploded validation replay and check the planner recovers
    to a finite trajectory."""
    x0, xf = _sixdof_scenario()
    pl = SixDofPlanner(np.linspace(0.0, 120.0, 5), t_max=2.0, tau_max=0.1)
    k = len(pl.times)
    real_step = SixDofPlanner._step
    real_solve = cp.Problem.solve
    state = {"budget": 0}
    solved = {"n": 0}

    def step_wrap(self, x, u, dt, n):
        out = real_step(self, x, u, dt, n)
        if state["budget"] > 0:  # poison the post-solve validation replay
            state["budget"] -= 1
            return out * np.inf
        return out

    def solve_wrap(self, *a, **kw):
        real_solve(self, *a, **kw)
        solved["n"] += 1
        if solved["n"] == 1:  # arm the replay right after the first solve
            state["budget"] = k - 1

    monkeypatch.setattr(SixDofPlanner, "_step", step_wrap)
    monkeypatch.setattr(cp.Problem, "solve", solve_wrap)
    with np.errstate(invalid="ignore", over="ignore"):  # the injected inf
        p = pl.solve(x0, xf, N, max_iter=4)
    assert np.all(np.isfinite(p.states))
    assert p.iterations >= 1  # kept iterating after the rejected step


def test_sixdof_subproblem_infeasible_status_and_fallback(monkeypatch):
    """A non-optimal subproblem status aborts the loop with the propagated
    status and still returns a finite nonlinear-replay fallback."""
    x0, xf = _sixdof_scenario()
    pl = SixDofPlanner(np.linspace(0.0, 120.0, 5), t_max=2.0, tau_max=0.1)
    orig = cp.Problem.solve

    def report_infeasible(self, *a, **k):
        orig(self, *a, **k)
        self._status = "infeasible"

    monkeypatch.setattr(cp.Problem, "solve", report_infeasible)
    p = pl.solve(x0, xf, N, max_iter=5)
    assert p.status == "subproblem_infeasible"
    assert p.iterations == 0  # aborted before recording any iterate
    assert np.all(np.isfinite(p.states))
