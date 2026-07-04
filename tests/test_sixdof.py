"""6-DOF attitude-coupled PTR receipts: the planner must discover the
braking slew (the thruster is body-fixed — coupling is physical), the
plan must replay on the nonlinear dynamics, and every bound must hold."""

import math

import numpy as np
import pytest

pytest.importorskip("cvxpy")

from podium import constants as const  # noqa: E402
from podium.core import quat  # noqa: E402
from podium.guidance.sixdof import SixDofPlanner, _deriv  # noqa: E402

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)


def _scenario():
    q0 = quat.normalize(np.array([1.0, 0.0, 0.0, 0.0]))
    # thruster (+x body) ends RETROGRADE (-y LVLH): the final braking
    # burn is physically possible; camera on -x body faces the target.
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


@pytest.fixture(scope="module")
def plan():
    x0, xf = _scenario()
    pl = SixDofPlanner(np.linspace(0.0, 360.0, 13), t_max=5.0,
                       tau_max=0.2)
    return pl, pl.solve(x0, xf, N), x0, xf


@pytest.mark.slow
def test_converges_and_replays_nonlinearly(plan):
    """Convergence with zero slack, and an INDEPENDENT nonlinear
    re-propagation of the planned controls hits the terminal box —
    dynamic feasibility, not just convex feasibility."""
    pl, p, x0, xf = plan
    assert p.status == "converged", p.status
    assert p.slack < 1e-9
    y = x0.copy()
    dts = np.diff(p.times)
    for i in range(len(dts)):
        y = pl._step(y, p.controls[i], dts[i], N)
    assert np.linalg.norm(y[0:3] - xf[0:3]) < 0.05
    assert np.linalg.norm(y[3:6] - xf[3:6]) < 0.005
    assert np.max(np.abs(y[6:10] - xf[6:10])) < 0.03
    assert np.max(np.abs(y[10:13])) < 0.003


@pytest.mark.slow
def test_planner_discovers_the_braking_slew(plan):
    """The coupling receipt: attitude swings ~90 deg over the horizon,
    and at the FINAL burn the body thrust axis points retrograde
    (-y LVLH) — braking. Nobody told the planner to slew; only the
    terminal attitude box and the body-fixed thruster did."""
    pl, p, _x0, _xf = plan
    ex = np.array([1.0, 0.0, 0.0])
    dir0 = quat.rotate(p.states[0, 6:10], ex)
    dirf = quat.rotate(p.states[-1, 6:10], ex)
    swing = math.degrees(math.acos(
        float(np.clip(np.dot(dir0, dirf), -1.0, 1.0))))
    assert swing > 60.0, swing
    # final interval with meaningful thrust: axis must oppose closing
    burn_nodes = np.flatnonzero(p.controls[:, 0] > 0.05)
    assert len(burn_nodes) >= 1
    k_last = int(burn_nodes[-1])
    axis = quat.rotate(p.states[k_last, 6:10], ex)
    assert axis[1] < -0.7, axis  # retrograde (-y): braking


@pytest.mark.slow
def test_bounds_and_dynamics_consistency(plan):
    pl, p, _x0, _xf = plan
    assert np.all(p.controls[:, 0] >= -1e-9)
    assert np.all(p.controls[:, 0] <= pl.t_max + 1e-9)
    assert np.all(np.abs(p.controls[:, 1:4]) <= pl.tau_max + 1e-9)
    # thrust acceleration in the dynamics is exactly (T/m) R(q) e1
    x = p.states[3]
    u = np.array([2.0, 0.0, 0.0, 0.0])
    d_on = _deriv(x, u, N, pl.mass, pl.inertia, pl.inertia_inv)
    d_off = _deriv(x, np.zeros(4), N, pl.mass, pl.inertia,
                   pl.inertia_inv)
    a_thr = d_on[3:6] - d_off[3:6]
    expect = (2.0 / pl.mass) * quat.rotate(x[6:10],
                                           np.array([1.0, 0.0, 0.0]))
    assert np.allclose(a_thr, expect, atol=1e-12)
    # quaternion norm preserved along the plan (renormalized RK4)
    norms = np.linalg.norm(p.states[:, 6:10], axis=1)
    assert np.max(np.abs(norms - 1.0)) < 1e-9