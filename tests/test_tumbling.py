"""Tumbling-dock study receipts: kinematic port matching through the
engine, rotating corridor held, dv-vs-rate monotonicity, and the
feasibility boundary."""

import math

import numpy as np
import pytest

cp = pytest.importorskip("cvxpy")

from podium import constants as const  # noqa: E402
from podium.guidance import tumbling  # noqa: E402
from podium.guidance.convex import plan_to_controller  # noqa: E402
from podium.sim import Scenario, circular_target, run  # noqa: E402

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)
X0 = np.array([0.0, -150.0, 0.0, 0.0, 0.0, 0.0])
TIMES = np.linspace(0.0, 480.0, 25)


def test_port_kinematics():
    tb = tumbling.Tumble(w_spin=0.01, rho_port=10.0, phase0=0.3)
    for t in (0.0, 37.5, 200.0):
        p = tb.port_state(t)
        assert abs(float(np.linalg.norm(p[0:3])) - 10.0) < 1e-12
        assert abs(float(np.linalg.norm(p[3:6])) - 0.1) < 1e-12
        # velocity is perpendicular to position (circular motion)
        assert abs(float(p[0:3] @ p[3:6])) < 1e-10


def test_plan_matches_rotating_port():
    tb = tumbling.Tumble(w_spin=0.005, rho_port=10.0)
    plan = tumbling.plan_tumbling_dock(TIMES, X0, tb, N)
    assert plan.status == "optimal"
    assert plan.terminal_pos_err < 1e-5
    assert plan.terminal_vel_err < 1e-6
    # corridor: final-half nodes inside the rotating cone
    tan_a = math.tan(math.radians(20.0))
    for i in range(len(TIMES) // 2 + 1, len(TIMES)):
        p = tb.port_state(float(TIMES[i]))[0:3]
        axis = p / np.linalg.norm(p)
        rel = plan.states[i, 0:3] - p
        perp = rel - (rel @ axis) * axis
        assert float(np.linalg.norm(perp)) <= tan_a * float(rel @ axis) + 1e-4


def test_envelope_findings_and_closure():
    """The phase-fixed envelope, with both study findings pinned:
    (1) a counterintuitive dip — matching a slow port beats nulling all
    motion; (2) strong monotone growth past the dip until the dv-capped
    envelope closes (0.05 rad/s infeasible at 0.35 m/s budget)."""
    rates = np.array([0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.2])
    table = tumbling.envelope_sweep(rates, TIMES, X0, N, rho_port=10.0,
                                    dv_max=0.35,
                                    fix_arrival_phase=-math.pi / 2)
    feas = [row for row in table if row["feasible"]]
    assert len(feas) >= 4
    # STUDY FINDING 2: even phase-fixed, dv(0.002) < dv(0) — matching a
    # small tangential port velocity is CHEAPER than nulling all motion,
    # because CW natural drift supplies velocity for free. The naive
    # "dv grows by rho*w" intuition is not a theorem for the full
    # trajectory problem; the honest claims are the measured table:
    assert feas[1]["dv"] < feas[0]["dv"]  # the counterintuitive dip
    # strong upward trend past the dip (measured 0.63 -> 2.14)
    assert feas[-1]["dv"] > 3.0 * feas[0]["dv"]
    for a, b in zip(feas[1:], feas[2:]):
        assert b["dv"] >= a["dv"], (a, b)  # monotone beyond the dip
    # regression pins (deterministic solver; loose windows)
    assert 0.55 < feas[0]["dv"] < 0.70
    assert 1.9 < feas[-1]["dv"] < 2.4
    # the envelope genuinely closes: 0.05 rad/s infeasible under the cap
    assert not table[-2]["feasible"]
    assert not table[-1]["feasible"]


def test_free_phase_is_nonmonotone_study_finding():
    """The confound itself, pinned as a receipt: with a FREE phase the
    fuel-vs-rate curve is non-monotone (arrival geometry dominates),
    which is why the envelope methodology fixes the arrival phase."""
    rates = np.array([0.0, 0.002, 0.005, 0.01, 0.02])
    table = tumbling.envelope_sweep(rates, TIMES, X0, N,
                                    rho_port=10.0, dv_max=0.35)
    feas = [row["dv"] for row in table if row["feasible"]]
    assert len(feas) >= 4
    diffs = np.diff(feas)
    assert float(np.min(diffs)) < -1e-3  # a genuine decrease exists


@pytest.mark.slow
def test_tumbling_plan_flies_and_arrives_on_port():
    """Engine receipt: fly the plan against the nonlinear truth; at t_f
    the chaser sits ON the (independently recomputed) rotating port with
    matched velocity, within the linearization budget at 150 m scale."""
    tb = tumbling.Tumble(w_spin=0.005, rho_port=10.0)
    plan = tumbling.plan_tumbling_dock(TIMES, X0, tb, N)
    assert plan.status == "optimal"

    class Shim:
        times = plan.times
        dvs = plan.dvs

    sc = Scenario(duration=float(TIMES[-1]), rv_target0=circular_target(A),
                  x_rel0=X0.copy(), dt_gnc=0.5, truth_substeps=2)
    tr = run(sc, plan_to_controller(Shim, tick_tol=0.25))
    xf = tr.x_rel[-1]
    port = tb.port_state(float(TIMES[-1]))
    assert float(np.linalg.norm(xf[0:3] - port[0:3])) < 1.0
    assert float(np.linalg.norm(xf[3:6] - port[3:6])) < 0.01