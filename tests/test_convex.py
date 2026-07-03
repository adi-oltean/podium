"""Layer-0 convex guidance receipts: analytic reproduction, transcription
exactness, constraint bite, DPP re-solve, ROE reconfiguration, and a
planned trajectory flown closed-loop through the engine."""

import math

import numpy as np
import pytest

cp = pytest.importorskip("cvxpy")

from podium import constants as const  # noqa: E402
from podium.core import cw, ya  # noqa: E402
from podium.core import roe as roe_mod  # noqa: E402
from podium.guidance import safety  # noqa: E402
from podium.guidance.convex import (  # noqa: E402
    ConeSpec,
    KozSpec,
    PlumeSpec,
    RendezvousPlanner,
    RoePlanner,
    plan_to_controller,
)
from podium.sim import Scenario, circular_target, mean_motion_of, run  # noqa: E402
from podium.sim import spec as sp  # noqa: E402

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)
X0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])
XF = np.zeros(6)


def replay(plan, n, e=0.0, theta0=0.0):
    """Propagate the planned burns through the same STM family."""
    x = plan.states[0].copy()
    for i in range(len(plan.times) - 1):
        x[3:6] += plan.dvs[i]
        dt = plan.times[i + 1] - plan.times[i]
        if e == 0.0:
            x = cw.stm(n, dt) @ x
        else:
            th = ya.propagate_true_anomaly(n, e, theta0, plan.times[i])
            x = ya.stm(n, e, th, dt) @ x
    x[3:6] += plan.dvs[-1]
    return x


def test_single_interval_reproduces_two_impulse():
    """K=1 with fixed endpoints has a unique feasible burn pair — it must
    equal the closed-form CW two-impulse solution."""
    tof = 1500.0
    planner = RendezvousPlanner(np.array([0.0, tof]))
    plan = planner.solve(X0, XF, n=N)
    dv1, dv2 = cw.two_impulse(X0, XF, N, tof)
    assert plan.status == "optimal"
    assert np.allclose(plan.dvs[0], dv1, atol=1e-6)
    assert np.allclose(plan.dvs[1], dv2, atol=1e-6)


def test_multinode_no_worse_than_two_impulse():
    tof = 3000.0
    dv1, dv2 = cw.two_impulse(X0, XF, N, tof)
    baseline = float(np.linalg.norm(dv1) + np.linalg.norm(dv2))
    planner = RendezvousPlanner(np.linspace(0.0, tof, 11))
    plan = planner.solve(X0, XF, n=N)
    assert plan.status == "optimal"
    assert plan.total_dv() <= baseline + 1e-6


def test_transcription_exactness_cw_and_ya():
    """Replaying the plan through the same STMs must hit the target at
    solver tolerance — the transcription itself has no discretization
    error."""
    times = np.linspace(0.0, 2500.0, 9)
    planner = RendezvousPlanner(times)
    for e, th0 in ((0.0, 0.0), (0.15, 1.1)):
        plan = planner.solve(X0, XF, n=N, e=e, theta0=th0)
        assert plan.status == "optimal"
        xf = replay(plan, N, e, th0)
        assert np.linalg.norm(xf[:3] - XF[:3]) < 1e-4
        assert np.linalg.norm(xf[3:] - XF[3:]) < 1e-7


def test_dpp_resolve_with_new_parameters():
    times = np.linspace(0.0, 2000.0, 6)
    planner = RendezvousPlanner(times)
    p1 = planner.solve(X0, XF, n=N)
    x0b = np.array([100.0, -600.0, 50.0, 0.0, 0.2, 0.0])
    p2 = planner.solve(x0b, XF, n=N)  # same compiled problem, new params
    assert p2.status == "optimal"
    xf = replay(p2, N)
    assert np.linalg.norm(xf[:3]) < 1e-4
    assert not np.allclose(p1.dvs, p2.dvs)


def test_approach_cone_bites_and_holds():
    """Start displaced off the corridor axis; the unconstrained plan cuts
    the corner, the constrained plan stays inside the cone."""
    times = np.linspace(0.0, 1800.0, 13)
    x0 = np.array([150.0, -800.0, 0.0, 0.0, 0.0, 0.0])
    cone = ConeSpec(apex=(0.0, 0.0, 0.0), axis=(0.0, -1.0, 0.0),
                    half_angle=math.radians(15.0), from_time=600.0)
    axis = np.array([0.0, -1.0, 0.0])
    proj = np.eye(3) - np.outer(axis, axis)
    tan_a = math.tan(cone.half_angle)

    def violation(plan):
        worst = -math.inf
        for i, t in enumerate(times):
            if t >= cone.from_time:
                r = plan.states[i, :3]
                worst = max(worst, np.linalg.norm(proj @ r) - tan_a * (axis @ r))
        return worst

    free = RendezvousPlanner(times).solve(x0, XF, n=N)
    assert violation(free) > 1.0  # the constraint is not vacuous
    coned = RendezvousPlanner(times, cone=cone).solve(x0, XF, n=N)
    assert coned.status == "optimal"
    assert violation(coned) < 1e-6
    assert coned.objective >= free.objective - 1e-9


def test_koz_hyperplanes_keep_distance():
    """Short pass-through transfer (nearly straight over the origin): with
    the KOZ every intermediate node stays outside the sphere (hyperplane
    implies distance since ||r|| >= n^T r), and fuel cost does not
    decrease."""
    times = np.linspace(0.0, 600.0, 13)
    x0 = np.array([0.0, -300.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, 300.0, 0.0, 0.0, 0.0, 0.0])
    koz = KozSpec(radius=200.0)

    free = RendezvousPlanner(times).solve(x0, xf, n=N)
    interior = range(1, len(times) - 1)
    min_free = min(np.linalg.norm(free.states[i, :3]) for i in interior)
    assert min_free < 100.0  # unconstrained path cuts through the sphere
    guarded = RendezvousPlanner(times, koz=koz).solve(x0, xf, n=N)
    assert guarded.status == "optimal"
    min_guard = min(np.linalg.norm(guarded.states[i, :3]) for i in interior)
    assert min_guard >= koz.radius - 1e-6
    assert guarded.objective >= free.objective - 1e-9


def test_plume_halfspace_blocks_target_pointing_burns():
    """Every plume constraint active in the final solve must hold; the
    arrival braking burn is exempt by design (Plan records the active
    (node, direction) pairs for exactly this transparency)."""
    times = np.linspace(0.0, 1500.0, 11)
    x0 = np.array([0.0, -400.0, 0.0, 0.0, 0.0, 0.0])
    plume = PlumeSpec(range_m=500.0)
    plan = RendezvousPlanner(times, plume=plume).solve(x0, XF, n=N)
    assert plan.status == "optimal"
    assert len(plan.plume_dirs) >= 5  # the zone genuinely covers nodes
    for i, d in plan.plume_dirs:
        assert i < len(times) - 1  # arrival node exempt
        assert float(d @ plan.dvs[i]) <= 1e-7


def test_l1_objective_solves_and_is_sparser_in_1norm():
    times = np.linspace(0.0, 2400.0, 9)
    l1 = RendezvousPlanner(times, objective="l1").solve(X0, XF, n=N)
    l2 = RendezvousPlanner(times, objective="l2").solve(X0, XF, n=N)
    assert l1.status == l2.status == "optimal"
    assert np.sum(np.abs(l1.dvs)) <= np.sum(np.abs(l2.dvs)) + 1e-9


@pytest.mark.slow
def test_planned_approach_flies_through_engine():
    """The whole Layer-0 story in one receipt: plan on the CW STM, fly
    closed-loop through the engine against the nonlinear truth, judge by
    specs. Terminal error budget = linearization error at 1 km scale
    (envelope tests bound it near C*sep^2/a ~ meters) + burn-time
    quantization on the 1 s grid."""
    tof = 2400.0
    times = np.linspace(0.0, tof, 9)
    hold = np.array([0.0, -30.0, 0.0, 0.0, 0.0, 0.0])
    planner = RendezvousPlanner(times, objective="l2")
    plan = planner.solve(X0, hold, n=N)

    specs = (
        sp.always_below("stay_inside", "range", 1100.0),
        sp.eventually_below("arrive", "range", 40.0),
        sp.final_between("terminal_speed", "speed", 0.0, 0.06),
    )
    sc = Scenario(
        duration=tof + 30.0,
        rv_target0=circular_target(A),
        x_rel0=X0.copy(),
        dt_gnc=1.0,
        truth_substeps=5,
    )
    assert abs(mean_motion_of(sc.rv_target0) - N) < 1e-9
    tr = run(sc, plan_to_controller(plan), specs)
    for name, margin in tr.spec_margins.items():
        assert margin > 0.0, f"{name}: {margin}"
    final = tr.channels()["range"][-1]
    assert abs(final - 30.0) < 8.0


def test_roe_planner_reconfigures_and_reports_safety():
    """Reconfigure to an e/i-aligned safe geometry; the plan must hit the
    target ROE through its own dynamics and the resulting geometry must
    show positive RN-plane margin."""
    a, n = A, N
    orbit = 2 * math.pi / n
    times = np.linspace(0.0, 1.5 * orbit, 7)
    roe0 = np.zeros(6)
    roef = np.array([0.0, 0.0, 4e-5, 0.0, 4e-5, 0.0])  # aligned de/di
    plan = RoePlanner(times).solve(roe0, roef, a=a, n=n, u0=0.3)
    assert plan.status == "optimal"
    # replay through the same dynamics
    r = plan.roes[0].copy()
    for i in range(len(times) - 1):
        u_i = 0.3 + n * plan.times[i]
        r = roe_mod.stm_keplerian(n, plan.times[i + 1] - plan.times[i]) @ (
            r + roe_mod.control_matrix(a, n, u_i) @ plan.dvs[i]
        )
    r = r + roe_mod.control_matrix(a, n, 0.3 + n * plan.times[-1]) @ plan.dvs[-1]
    # tolerance = solver feasibility tolerance (Clarabel default ~1e-8),
    # not machine precision: the replay checks the solver's constraints.
    assert np.allclose(r, roef, atol=1e-7)
    assert safety.rn_margin(roef, a, 200.0) > 0.0
    # dv sanity: dominant cost is the e-vector change, |d(de)|*n*a/2
    expected = np.linalg.norm(roef[2:4]) * n * a / 2.0
    assert plan.total_dv() < 3.0 * expected
