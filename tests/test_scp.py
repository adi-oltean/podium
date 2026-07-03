"""PTR/SCvx* receipts: convex reduction, true nonconvex feasibility,
exact-flow CTCS cuts, penalty ramp, and the engine flight."""

import math

import numpy as np
import pytest

cp = pytest.importorskip("cvxpy")

from podium import constants as const  # noqa: E402
from podium.core import cw  # noqa: E402
from podium.guidance.convex import (  # noqa: E402
    KozSpec,
    RendezvousPlanner,
    plan_to_controller,
)
from podium.guidance.scp import PtrDockingPlanner, ScpResult  # noqa: E402
from podium.sim import Scenario, circular_target, run  # noqa: E402

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)


def dense_min_range(res: ScpResult, n: float, samples: int = 1000) -> float:
    """Independent dense check: exact flow, fine grid, min ||r||."""
    worst = math.inf
    for k in range(len(res.times) - 1):
        dt = res.times[k + 1] - res.times[k]
        xk = res.states[k] + np.concatenate([np.zeros(3), res.dvs[k]])
        for i in range(samples // (len(res.times) - 1) + 1):
            tau = dt * i / (samples // (len(res.times) - 1) + 1)
            p = (cw.stm(n, tau) @ xk)[0:3]
            worst = min(worst, float(np.linalg.norm(p)))
    return worst


def test_convex_problem_reduces_to_layer0():
    """No active KOZ: PTR must converge fast to the Layer-0 optimum."""
    times = np.linspace(0.0, 2000.0, 9)
    x0 = np.array([0.0, -800.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, -300.0, 0.0, 0.0, 0.0, 0.0])
    base = RendezvousPlanner(times).solve(x0, xf, n=N)
    ptr = PtrDockingPlanner(times, koz_radius=50.0).solve(x0, xf, N)
    assert ptr.status == "converged"
    assert ptr.slack_final <= 1e-6
    assert ptr.total_dv() <= base.total_dv() + 1e-4
    assert ptr.iterations <= 6


def test_nonconvex_passage_beats_hyperplane_heuristic():
    """Far-side passage: PTR satisfies the TRUE sphere constraint and
    costs no more than Layer-0's two-pass hyperplane approximation."""
    times = np.linspace(0.0, 600.0, 13)
    x0 = np.array([0.0, -300.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, 300.0, 0.0, 0.0, 0.0, 0.0])
    r = 200.0
    heur = RendezvousPlanner(times, koz=KozSpec(radius=r)).solve(x0, xf, n=N)
    ptr = PtrDockingPlanner(times, koz_radius=r).solve(x0, xf, N)
    assert ptr.status == "converged"
    assert ptr.slack_final <= 1e-6
    # true nonconvex constraint at nodes
    for k in range(1, len(times) - 1):
        assert float(np.linalg.norm(ptr.states[k, 0:3])) >= r - 1e-3
    # within trust-region regularization bias of the heuristic's cost
    # (measured gap 0.009%); the point is the TRUE constraint holds
    assert ptr.total_dv() <= heur.total_dv() * 1.001
    # and continuously, by independent dense check
    assert dense_min_range(ptr, N) >= r - 0.5


def test_ctcs_cuts_catch_intersample_dip():
    """Coarse grid (few nodes): node constraints alone let the coast arc
    dip inside the sphere between nodes; the exact-flow cuts must catch
    and eliminate it to a clean independent dense check."""
    times = np.linspace(0.0, 700.0, 5)  # very coarse: 4 arcs
    x0 = np.array([0.0, -350.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, 350.0, 0.0, 0.0, 0.0, 0.0])
    r = 250.0
    ptr = PtrDockingPlanner(times, koz_radius=r, ctcs_samples=40)
    res = ptr.solve(x0, xf, N)
    assert res.status == "converged"
    assert res.n_cuts > 0  # the dip genuinely happened and was cut
    assert res.dense_violation <= 1e-3
    assert dense_min_range(res, N) >= r - 0.5


def test_scvxstar_penalty_ramp():
    """Start with a uselessly small penalty: the SCvx*-style ramp must
    grow it and still converge to a feasible solution."""
    times = np.linspace(0.0, 600.0, 13)
    x0 = np.array([0.0, -300.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, 300.0, 0.0, 0.0, 0.0, 0.0])
    ptr = PtrDockingPlanner(times, koz_radius=200.0, w_pen0=1e-6)
    res = ptr.solve(x0, xf, N)
    assert res.status == "converged"
    assert res.w_pen_final > 1e-6  # ramp engaged
    assert res.slack_final <= 1e-6


def test_deterministic():
    times = np.linspace(0.0, 600.0, 13)
    x0 = np.array([0.0, -300.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, 300.0, 0.0, 0.0, 0.0, 0.0])
    r1 = PtrDockingPlanner(times, koz_radius=200.0).solve(x0, xf, N)
    r2 = PtrDockingPlanner(times, koz_radius=200.0).solve(x0, xf, N)
    assert np.array_equal(r1.dvs, r2.dvs)
    assert r1.iterations == r2.iterations


@pytest.mark.slow
def test_ptr_plan_flies_through_engine():
    """The SCP plan holds up against the nonlinear truth: keep-out
    respected within the linearization allowance, arrival achieved."""
    times = np.linspace(0.0, 600.0, 13)
    x0 = np.array([0.0, -300.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, 300.0, 0.0, 0.0, 0.0, 0.0])
    r = 200.0
    res = PtrDockingPlanner(times, koz_radius=r).solve(x0, xf, N)
    assert res.status == "converged"

    class PlanShim:
        times = res.times
        dvs = res.dvs

    sc = Scenario(duration=620.0, rv_target0=circular_target(A),
                  x_rel0=x0.copy(), dt_gnc=1.0, truth_substeps=4)
    tr = run(sc, plan_to_controller(PlanShim))
    rng = tr.channels()["range"]
    # CW-vs-truth error at 300 m scale over 600 s is well under 2 m
    assert float(np.min(rng[5:-5])) >= r - 2.0
    assert abs(float(rng[-1]) - 300.0) < 5.0