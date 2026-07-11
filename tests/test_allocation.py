"""Thruster-allocation receipts (#39): the 6-DOF guidance's commanded
body wrench is realized by push-only discrete thrusters.

The physics that makes this non-trivial is non-negativity — a plain
pseudoinverse returns negative "thrust" no thruster can produce. The
min-propellant LP respects u >= 0 and minimizes total thrust; the
least-squares fallback reports the closest achievable wrench when a
demand exceeds the cluster's authority.
"""

import math

import numpy as np
import pytest

pytest.importorskip("scipy")

from podium.control import allocation as al
from podium.control.allocation import allocate, standard_cluster


def test_cluster_has_full_wrench_authority():
    """B is rank 6: the cluster can produce any wrench direction."""
    cfg = standard_cluster()
    b = cfg.effectiveness()
    assert b.shape == (6, 24)
    assert np.linalg.matrix_rank(b, tol=1e-9) == 6


def test_directions_are_unit_and_torque_arm_correct():
    cfg = standard_cluster(half=1.5)
    assert np.allclose(np.linalg.norm(cfg.directions, axis=1), 1.0)
    b = cfg.effectiveness()
    # column i torque rows equal r_i x d_i exactly
    for i in range(cfg.n):
        expect = np.cross(cfg.positions[i], cfg.directions[i])
        assert np.allclose(b[3:6, i], expect)


def test_feasible_wrenches_reproduce_nonnegative():
    """Random wrenches inside the cluster's authority reproduce to 1e-9
    with strictly non-negative thrusts."""
    cfg = standard_cluster()
    b = cfg.effectiveness()
    rng = np.random.default_rng(39)
    for _ in range(200):
        # a wrench known feasible: a non-negative combination of columns
        u_true = rng.uniform(0.0, 1.0, cfg.n) * (rng.random(cfg.n) < 0.3)
        w = b @ u_true
        a = allocate(cfg, w)
        assert a.feasible, w
        assert np.all(a.u >= -1e-12)
        assert a.residual < 1e-9
        # propellant is minimal: no more than the witness combination
        assert a.propellant <= u_true.sum() + 1e-9


def test_pure_couple_is_realizable():
    """A pure torque with ZERO net force — the classic RCS requirement —
    is realized by opposing thrusters (net force cancels)."""
    cfg = standard_cluster()
    torque = np.array([0.0, 0.0, 0.0, 0.03, -0.02, 0.05])
    a = allocate(cfg, torque)
    assert a.feasible
    assert np.linalg.norm(a.realized[0:3]) < 1e-9   # no net force
    assert np.allclose(a.realized[3:6], torque[3:6], atol=1e-9)


def test_pseudoinverse_would_go_negative():
    """The witness that non-negativity matters: the minimum-norm
    (pseudoinverse) solution for a realizable wrench contains negative
    'thrust', which the LP avoids while still reproducing the wrench."""
    cfg = standard_cluster()
    b = cfg.effectiveness()
    w = np.array([0.4, 0.0, 0.0, 0.0, 0.1, 0.0])
    u_pinv = np.linalg.pinv(b) @ w
    assert u_pinv.min() < -1e-6            # pseudoinverse is infeasible
    a = allocate(cfg, w)
    assert a.feasible and a.u.min() >= -1e-12
    assert a.residual < 1e-9


def test_infeasible_demand_reports_closest():
    """A wrench beyond the per-thruster bound is flagged infeasible and
    the closest achievable wrench + residual are reported (no lying)."""
    cfg = standard_cluster()
    # only 4 thrusters point +x; capped at u=1 each, max Fx is 4, so a
    # 20 N +x demand cannot be met
    big = np.array([20.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    a = allocate(cfg, big, u_max=1.0)
    assert not a.feasible
    assert a.residual > 1.0
    assert np.all(a.u <= 1.0 + 1e-9)


@pytest.mark.slow
def test_sixdof_plan_wrenches_allocate():
    """Every node of a 6-DOF PTR plan's (thrust, torque) command
    allocates feasibly onto the cluster — the guidance output is
    hardware-realizable."""
    pytest.importorskip("cvxpy")
    from podium import constants as const  # noqa: PLC0415
    from podium.core import quat  # noqa: PLC0415
    from podium.guidance.sixdof import SixDofPlanner  # noqa: PLC0415

    n = math.sqrt(const.MU_EARTH / 6_778_137.0**3)
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
    pl = SixDofPlanner(np.linspace(0.0, 360.0, 13), t_max=5.0,
                       tau_max=0.2)
    plan = pl.solve(x0, xf, n)
    assert plan.status == "converged"

    # cluster sized so t_max thrust + tau_max torque are within authority
    cfg = standard_cluster(half=2.0)
    for k in range(len(plan.controls)):
        thrust, tau = plan.controls[k, 0], plan.controls[k, 1:4]
        wrench = np.concatenate([[thrust, 0.0, 0.0], tau])
        a = allocate(cfg, wrench, u_max=50.0)
        assert a.feasible, (k, wrench)
        assert a.residual < 1e-7


def test_bounded_lsq_is_no_worse_than_clipped_nnls():
    """The infeasible fallback is a TRUE box-bounded least-squares fit
    (0 <= u <= u_max), not an unbounded NNLS clipped to u_max -- so it stays
    within bounds and its residual is never worse than the clip-after-NNLS it
    replaced."""
    from scipy.optimize import nnls  # noqa: PLC0415

    cfg = standard_cluster()
    b = cfg.effectiveness()
    w = np.array([20.0, 3.0, -2.0, 0.5, 0.0, 0.0])   # beyond authority at u_max=1
    a = allocate(cfg, w, u_max=1.0)
    assert not a.feasible
    assert np.all(a.u >= -1e-12) and np.all(a.u <= 1.0 + 1e-9)
    u_nnls, _ = nnls(b, w)
    r_clip = float(np.linalg.norm(b @ np.minimum(u_nnls, 1.0) - w))
    assert a.residual <= r_clip + 1e-9


def test_config_validation():
    with pytest.raises(ValueError, match="both be"):
        al.ThrusterConfig(positions=np.zeros((4, 3)),
                          directions=np.zeros((3, 3)))


def test_zero_direction_rejected():
    """A [0, 0, 0] thrust direction cannot be normalized — it would inject a
    NaN column into the effectiveness matrix — so it is rejected up front."""
    with pytest.raises(ValueError, match="nonzero"):
        al.ThrusterConfig(positions=np.zeros((2, 3)),
                          directions=np.array([[1.0, 0.0, 0.0],
                                               [0.0, 0.0, 0.0]]))
