"""Layer-0 convex guidance receipts: analytic reproduction, transcription
exactness, constraint bite, DPP re-solve, ROE reconfiguration, and a
planned trajectory flown closed-loop through the engine."""

import math

import numpy as np
import pytest

cp = pytest.importorskip("cvxpy")

from podium import constants as const  # noqa: E402
from podium.control import lqr  # noqa: E402
from podium.core import cw, ya  # noqa: E402
from podium.core import roe as roe_mod  # noqa: E402
from podium.guidance import convex as ConvexMod  # noqa: E402
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


# --- Layer-0 follow-ups (#9) -------------------------------------------


def test_lcvx_annulus_lossless_and_exact():
    """Thrust-annulus finite burn: the relaxation must come back lossless
    (audited, not assumed — discrete-time LCvx bounds non-tight nodes by
    the state dimension), the annulus must genuinely bite, and replaying
    the ZOH profile through the exact discretization must hit the target
    at solver tolerance."""
    # Normal scenario: near-minimum-time transfer, so thrust is needed
    # throughout and the primer never vanishes.
    annulus = ConvexMod.AnnulusSpec(rho_min=0.008, rho_max=0.02)
    times = np.linspace(0.0, 900.0, 17)
    x0 = np.array([0.0, -2000.0, 0.0, 0.0, 0.0, 0.0])
    planner = ConvexMod.FiniteBurnPlanner(times, annulus)
    plan = planner.solve(x0, XF, n=N)
    assert plan.status == "optimal"
    assert plan.controllable
    mags = np.linalg.norm(plan.u, axis=1)
    assert np.all(mags >= annulus.rho_min - 1e-7)
    assert np.all(mags <= annulus.rho_max + 1e-7)
    # losslessness audit clean, and the lower bound genuinely binding
    # (measured: 12 of 16 nodes ride rho_min, max gap ~3e-9)
    assert len(plan.lcvx_inactive) == 0
    assert plan.lcvx_max_gap < 1e-5 * annulus.rho_max
    assert int(np.sum(mags < annulus.rho_min + 1e-6)) >= 8
    # replay through the exact ZOH dynamics
    ad, bd = lqr.cw_discrete(N, float(times[1] - times[0]))
    x = plan.states[0].copy()
    for i in range(len(times) - 1):
        x = ad @ x + bd @ plan.u[i]
    assert np.linalg.norm(x - XF) < 1e-5
    # the lower bound is not vacuous: without it the optimum coasts
    free = ConvexMod.FiniteBurnPlanner(
        times, ConvexMod.AnnulusSpec(rho_min=0.0, rho_max=0.02)
    ).solve(x0, XF, n=N)
    assert float(np.min(np.linalg.norm(free.u, axis=1))) < annulus.rho_min
    assert plan.objective >= free.objective - 1e-9


def test_lcvx_audit_catches_degenerate_relaxation():
    """Excess-capacity problem (forced minimum fuel far above what the
    transfer needs — the discrete analogue of coast arcs): the relaxation
    goes loose and the shipped audit must say so, loudly. This is the
    audit doing its job: such a solution is not a valid thrust profile."""
    annulus = ConvexMod.AnnulusSpec(rho_min=0.004, rho_max=0.02)
    times = np.linspace(0.0, 1500.0, 16)
    x0 = np.array([0.0, -2000.0, 0.0, 0.0, 0.0, 0.0])
    plan = ConvexMod.FiniteBurnPlanner(times, annulus).solve(x0, XF, n=N)
    assert plan.status == "optimal"  # the SOCP is fine; the audit is not
    assert len(plan.lcvx_inactive) > 6  # exceeds the n_x theory bound
    assert plan.lcvx_max_gap > 1e-4  # measured ~1.6e-3


def test_passive_safety_scenarios_protect_failure_drifts():
    """Breger-How: without the constraint some failure drift enters the
    keep-out sphere; with it, every protected (node, sample) drift stays
    outside — and since the hyperplane normal is unit, satisfaction
    implies true distance."""
    times = np.linspace(0.0, 2000.0, 11)
    x0 = np.array([0.0, -1200.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, -150.0, 0.0, 0.0, 0.0, 0.0])
    spec = ConvexMod.PassiveSafetySpec(
        radius=200.0, horizon=2000.0, n_samples=6,
        failure_nodes=tuple(range(1, 6)),
    )
    taus = np.linspace(0.0, spec.horizon, spec.n_samples + 1)[1:]

    def worst_drift(plan):
        worst = math.inf
        for j in spec.failure_nodes:
            for tau in taus:
                xd = cw.stm(N, float(tau)) @ plan.states[j]
                worst = min(worst, float(np.linalg.norm(xd[:3])))
        return worst

    free = RendezvousPlanner(times).solve(x0, xf, n=N)
    assert worst_drift(free) < spec.radius  # scenario genuinely dangerous
    guarded = RendezvousPlanner(times, passive_safety=spec).solve(x0, xf, n=N)
    assert guarded.status == "optimal"
    assert worst_drift(guarded) >= spec.radius - 1e-6
    assert guarded.objective >= free.objective - 1e-9


def test_qp_tracking_objective_dpp_resolve():
    """QP tracking: the compiled problem re-solves against different
    Parameter references, and tracking a reference beats ignoring it."""
    times = np.linspace(0.0, 1200.0, 9)
    planner = RendezvousPlanner(times, objective="qp_tracking",
                                track_state_weight=1e-4,
                                track_control_weight=1.0)
    x0 = np.array([100.0, -500.0, 0.0, 0.0, 0.0, 0.0])

    def straight_ref():
        ref = np.zeros((len(times), 6))
        for i in range(len(times)):
            frac = i / (len(times) - 1)
            ref[i, 0:3] = (1 - frac) * x0[0:3]
        return ref

    ref = straight_ref()
    p_track = planner.solve(x0, XF, n=N, x_ref=ref)
    p_zero = planner.solve(x0, XF, n=N)  # same compiled problem, ref=0
    assert p_track.status == p_zero.status == "optimal"

    def tracking_error(plan):
        return float(np.sum((plan.states - ref) ** 2))

    assert tracking_error(p_track) < tracking_error(p_zero)
    # both still satisfy the hard terminal condition
    for p in (p_track, p_zero):
        assert np.linalg.norm(replay(p, N)[:3]) < 1e-4


def test_roe_safe_set_terminal():
    """Convex e/i safe-set terminal: alignment cones + minimum magnitudes
    hold, |da| bounded, and the exact RN-plane scan confirms the achieved
    geometry clears the keep-out radius."""
    a, n = A, N
    orbit = 2 * math.pi / n
    times = np.linspace(0.0, 1.5 * orbit, 7)
    spec = ConvexMod.SafeSetSpec(direction=(1.0, 0.0), e_min=4e-5,
                                 i_min=4e-5, cone_angle=math.radians(10.0),
                                 da_tol=1e-7)
    plan = RoePlanner(times, safe_set=spec).solve(
        np.zeros(6), None, a=a, n=n, u0=0.3
    )
    assert plan.status == "optimal"
    # terminal post-burn ROE
    rT = plan.roes[-1] + roe_mod.control_matrix(
        a, n, 0.3 + n * plan.times[-1]
    ) @ plan.dvs[-1]
    de, di = rT[2:4], rT[4:6]
    uhat = np.array([1.0, 0.0])
    assert uhat @ de >= spec.e_min - 1e-9
    assert uhat @ di >= spec.i_min - 1e-9
    tan_c = math.tan(spec.cone_angle)
    assert np.linalg.norm(de - (uhat @ de) * uhat) <= tan_c * (uhat @ de) + 1e-9
    assert np.linalg.norm(di - (uhat @ di) * uhat) <= tan_c * (uhat @ di) + 1e-9
    assert abs(rT[0]) <= spec.da_tol + 1e-9
    # the receipt: exact scan of the achieved geometry
    assert safety.rn_margin(rT, a, 150.0) > 0.0
    expected = np.linalg.norm(de) * n * a / 2.0
    assert plan.total_dv() < 3.0 * expected


# --- Layer-0 residuals (#12) -------------------------------------------


def test_ya_discrete_reduces_and_composes():
    """YA ZOH: equals cw_discrete at e=0 and satisfies the two-interval
    composition identity Ad2 Ad1 / Ad2 Bd1 + Bd2 == full-interval maps."""
    dt = 120.0
    a_cw, b_cw = lqr.cw_discrete(N, dt)
    a_ya, b_ya = lqr.ya_discrete(N, 0.0, 0.7, dt)
    assert np.allclose(a_ya, a_cw, rtol=1e-9, atol=1e-12)
    assert np.allclose(b_ya, b_cw, rtol=1e-6, atol=1e-9)
    e, th0 = 0.2, 1.1
    a1, b1 = lqr.ya_discrete(N, e, th0, dt)
    th1 = ya.propagate_true_anomaly(N, e, th0, dt)
    a2, b2 = lqr.ya_discrete(N, e, th1, dt)
    a_full, b_full = lqr.ya_discrete(N, e, th0, 2 * dt)
    assert np.allclose(a2 @ a1, a_full, rtol=1e-9, atol=1e-9)
    assert np.allclose(a2 @ b1 + b2, b_full, rtol=1e-5, atol=1e-7)


def test_finite_burn_eccentric():
    """LCvx on YA dynamics (e=0.15): solves, and replaying the ZOH
    profile through the same per-interval maps hits the target."""
    annulus = ConvexMod.AnnulusSpec(rho_min=0.008, rho_max=0.02)
    times = np.linspace(0.0, 900.0, 13)
    x0 = np.array([0.0, -2000.0, 0.0, 0.0, 0.0, 0.0])
    e, th0 = 0.15, 0.9
    plan = ConvexMod.FiniteBurnPlanner(times, annulus).solve(
        x0, XF, n=N, e=e, theta0=th0
    )
    assert plan.status == "optimal"
    assert plan.controllable
    x = plan.states[0].copy()
    dt = float(times[1] - times[0])
    for i in range(len(times) - 1):
        th_i = ya.propagate_true_anomaly(N, e, th0, float(times[i]))
        ad, bd = lqr.ya_discrete(N, e, th_i, dt)
        x = ad @ x + bd @ plan.u[i]
    assert np.linalg.norm(x - XF) < 1e-5
    assert len(plan.lcvx_inactive) <= 6


def test_primer_certificate_separates_normal_from_degenerate():
    """Normality certificate from the duals: primer/dt stays O(0.1) on
    the normal problem and collapses (~1e-7) on the degenerate one —
    and at tight interior-slack nodes the primer equals dt exactly
    (Gamma-stationarity), which pins the dual convention."""
    x0 = np.array([0.0, -2000.0, 0.0, 0.0, 0.0, 0.0])
    normal = ConvexMod.FiniteBurnPlanner(
        np.linspace(0.0, 900.0, 17), ConvexMod.AnnulusSpec(0.008, 0.02)
    ).solve(x0, XF, n=N)
    degen = ConvexMod.FiniteBurnPlanner(
        np.linspace(0.0, 1500.0, 16), ConvexMod.AnnulusSpec(0.004, 0.02)
    ).solve(x0, XF, n=N)
    assert normal.primer_certificate() > 0.01
    assert degen.primer_certificate() < 1e-4
    # Gamma interior at the first node of the normal problem => primer == dt
    dt = float(normal.times[1] - normal.times[0])
    assert abs(normal.primer_norms[0] - dt) < 1e-3 * dt


def test_min_time_presolve_yields_normal_lcvx():
    """Coast-arc antidote: bisect the minimum feasible time, then LCvx at
    1.15x t_min passes the losslessness audit."""
    x0 = np.array([0.0, -2000.0, 0.0, 0.0, 0.0, 0.0])
    rho_max = 0.02
    t_min = ConvexMod.find_min_time(x0, XF, N, rho_max, k=12,
                                    t_lo=200.0, t_hi=2000.0)
    assert 300.0 < t_min < 1200.0
    # just below is infeasible (bisection bracket is real)
    tight = ConvexMod.FiniteBurnPlanner(
        np.linspace(0.0, 0.9 * t_min, 13), ConvexMod.AnnulusSpec(0.0, rho_max)
    ).solve(x0, XF, n=N)
    assert tight.status not in ("optimal",)
    plan = ConvexMod.FiniteBurnPlanner(
        np.linspace(0.0, 1.15 * t_min, 13),
        ConvexMod.AnnulusSpec(0.3 * rho_max, rho_max),
    ).solve(x0, XF, n=N)
    assert plan.status == "optimal"
    assert plan.lcvx_max_gap < 1e-5 * rho_max
    assert len(plan.lcvx_inactive) == 0


def test_passive_safety_dense_margins_reported():
    times = np.linspace(0.0, 2000.0, 11)
    x0 = np.array([0.0, -1200.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([0.0, -150.0, 0.0, 0.0, 0.0, 0.0])
    spec = ConvexMod.PassiveSafetySpec(
        radius=200.0, horizon=2000.0, n_samples=6,
        failure_nodes=tuple(range(1, 6)),
    )
    plan = RendezvousPlanner(times, passive_safety=spec).solve(x0, xf, n=N)
    assert set(plan.ps_margins) == set(spec.failure_nodes)
    # sparse samples guarantee >= radius at the samples; the dense report
    # may dip slightly between them — it must be reported, and small
    for j, margin in plan.ps_margins.items():
        assert margin > -10.0, (j, margin)


@pytest.mark.slow
def test_mib_quantized_plan_flies():
    """MIB bridge: sigma-delta quantization keeps the residual bounded by
    half a click per axis, and the quantized plan still arrives when
    flown through the engine against the nonlinear truth."""
    tof = 2400.0
    times = np.linspace(0.0, tof, 9)
    hold = np.array([0.0, -30.0, 0.0, 0.0, 0.0, 0.0])
    plan = RendezvousPlanner(times, objective="l2").solve(X0, hold, n=N)
    q = 0.002  # m/s per click (~1 N s on a 500 kg vehicle)
    qplan = ConvexMod.quantize_plan(plan, q)
    # every burn is an integer number of clicks
    assert np.allclose(np.round(qplan.dvs / q) * q, qplan.dvs, atol=1e-12)
    # sigma-delta: cumulative per-axis error bounded by half a click
    resid = np.sum(plan.dvs - qplan.dvs, axis=0)
    assert np.all(np.abs(resid) <= 0.5 * q + 1e-12)
    sc = Scenario(
        duration=tof + 30.0,
        rv_target0=circular_target(A),
        x_rel0=X0.copy(),
        dt_gnc=1.0,
        truth_substeps=5,
    )
    tr = run(sc, plan_to_controller(qplan))
    final = tr.channels()["range"][-1]
    assert abs(final - 30.0) < 12.0  # baseline test allows 8 m unquantized
