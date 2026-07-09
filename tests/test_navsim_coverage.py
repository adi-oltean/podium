"""Coverage-closing receipts for nav / sim / core / dynamics.

Each test pins the documented BEHAVIOUR of an otherwise-unexercised
branch (a guard clause, a rarely-taken sign correction, an optional
measurement mode), not merely the line: a near-singular / corrupted /
boundary input is asserted to be handled as the module promises, never
to produce a NaN or an exception where finiteness is guaranteed.
"""

import math

import numpy as np
import pytest

from podium import constants as const
from podium.core import cw as cw_mod
from podium.core import integrators, ya
from podium.dynamics import nonlinear as nl
from podium.nav import ekf, sensors
from podium.sim import circular_target, engine
from podium.sim import spec as spec_mod

A = 6_778_137.0
N = math.sqrt(const.MU_EARTH / A**3)


# --------------------------------------------------------------------------
# core/integrators.py:31 — euler_step (the cheap onboard predictor)
# --------------------------------------------------------------------------
def test_euler_step_is_the_tangent_line_predictor():
    """One explicit-Euler step is exactly the tangent-line extrapolation
    x + dt f(t, x); on a constant field it is even exact, and it is
    first-order (not RK4) so it carries the expected O(dt^2) local error."""
    # constant derivative field: Euler is exact
    const_f = lambda t, x: np.array([2.0, -3.0])  # noqa: E731
    x0 = np.array([1.0, 1.0])
    assert np.allclose(integrators.euler_step(const_f, 0.0, x0, 0.5),
                       x0 + 0.5 * np.array([2.0, -3.0]))

    # linear decay x' = -x: one step is the tangent line, and it differs
    # from the true e^{-dt} decay by the first-order truncation term
    decay = lambda t, x: -x  # noqa: E731
    dt = 0.1
    stepped = integrators.euler_step(decay, 0.0, np.array([1.0]), dt)
    assert np.allclose(stepped, np.array([1.0 - dt]))
    true = math.exp(-dt)
    assert stepped[0] < true  # tangent falls below the convex decay curve
    assert abs(stepped[0] - true) < dt * dt  # local error is O(dt^2)


# --------------------------------------------------------------------------
# nav/sensors.py:74 — RelGnss.start zero-bias branch
# --------------------------------------------------------------------------
def test_relgnss_zero_bias_budget_draws_no_bias():
    """With the constant-bias budget at its default 0, start() draws the
    zero vector (no per-run offset) rather than sampling the rng — so the
    measurement stream is unbiased about truth."""
    gnss = sensors.RelGnss(pos_std=1.0, vel_std=0.01)  # bias_pos_std == 0.0
    rng = np.random.default_rng(0)
    bias = gnss.start(rng)
    assert np.array_equal(bias, np.zeros(3))

    x = np.array([5.0, -400.0, 3.0, 0.1, 0.0, -0.2])
    zs = np.array([gnss.measure(x, rng, bias) for _ in range(4000)])
    # unbiased: sample mean of the position error collapses toward 0
    assert np.allclose(np.mean(zs[:, 0:3] - x[0:3], axis=0), 0.0, atol=0.1)


# --------------------------------------------------------------------------
# nav/ekf.py:183-184 — position+velocity measurement configuration
# --------------------------------------------------------------------------
def test_ekf_position_velocity_measurement_mode():
    """Passing r_vel selects the full 6-state measurement (H = I6) with a
    per-block diagonal R; a velocity measurement then corrects the
    velocity estimate, which the position-only mode cannot do."""
    f = ekf.RelNavEkf(N, dt=2.0, q_accel=1e-8, r_pos=5.0, r_vel=0.05)
    assert np.array_equal(f.h, ekf.H_POSVEL)
    assert np.allclose(np.diag(f.r), [25.0, 25.0, 25.0, 0.0025, 0.0025, 0.0025])
    assert f.last_nu.shape == (6,) and f.last_s.shape == (6, 6)

    # a 6-vector measurement that disagrees in velocity pulls the velocity
    # estimate toward it (impossible with the position-only H)
    z = np.array([0.0, 0.0, 0.0, 1.0, -2.0, 0.5])
    est = f.step(z)
    assert np.all(np.abs(est[3:6]) > 0.0)
    assert np.sign(est[3]) == 1.0 and np.sign(est[4]) == -1.0


# --------------------------------------------------------------------------
# nav/ekf.py:202 — commanded dv fed through step()'s prediction
# --------------------------------------------------------------------------
def test_ekf_step_feeds_commanded_burn_through_prediction():
    """A known impulsive burn passed to step() is added to the velocity
    BEFORE the time update, so the maneuver propagates the estimate
    instead of being mistaken for state error. The returned (post-update,
    pre-prediction) estimate is unchanged by dv; the next predicted state
    differs by exactly Phi applied to the velocity impulse."""
    z = np.array([12.0, -5.0, 3.0])
    dv = np.array([0.1, -0.2, 0.05])
    f_dv = ekf.RelNavEkf(N, dt=3.0, q_accel=1e-9, r_pos=4.0,
                         x0=np.array([10.0, -4.0, 2.0, 0.0, 0.0, 0.0]))
    f_no = ekf.RelNavEkf(N, dt=3.0, q_accel=1e-9, r_pos=4.0,
                         x0=np.array([10.0, -4.0, 2.0, 0.0, 0.0, 0.0]))
    est_dv = f_dv.step(z, dv=dv)
    est_no = f_no.step(z, dv=None)
    assert np.allclose(est_dv, est_no)  # dv does not touch the estimate
    expected = f_no.phi @ np.concatenate([np.zeros(3), dv])
    assert np.allclose(f_dv.x - f_no.x, expected)


# --------------------------------------------------------------------------
# nav/ekf.py:92 — covariance-repair clamp in the sequential update
# --------------------------------------------------------------------------
def test_sequential_update_clamp_survives_corrupted_covariance():
    """The scalar innovation variance s = P[m,m] + r_var is clamped up to
    r_var so the division is provably safe even for a NON-PSD covariance
    (a negative diagonal that a valid filter never has). The update then
    stays finite and symmetric instead of dividing by a non-positive
    number and blowing up."""
    x = np.zeros(6)
    p = np.eye(6)
    p[0, 0] = -1.0  # corrupt: invalid covariance, P[0,0] + r_var < r_var
    r_var = 0.01
    assert p[0, 0] + r_var < r_var  # the clamp branch is genuinely taken
    xs, ps = ekf.update_sequential(x, p, np.array([1.0, 0.0, 0.0]), r_var)
    assert np.all(np.isfinite(xs))
    assert np.all(np.isfinite(ps))
    assert np.allclose(ps, ps.T)  # symmetrised despite the bad input


# --------------------------------------------------------------------------
# sim/spec.py:59,63,91 — eventually_above semantics + unknown-kind guard
# --------------------------------------------------------------------------
def test_spec_eventually_above_margin_and_unknown_kind():
    """F_[window](channel >= lo) is satisfied by the BEST sample in the
    window, so its robustness is max(signal - lo); and an unrecognised
    spec kind is rejected loudly rather than silently scored."""
    t = np.array([0.0, 1.0, 2.0, 3.0])
    sig = np.array([1.0, 4.0, 2.0, 3.0])
    sp = spec_mod.eventually_above("reach", "ch", lo=3.5)
    assert sp.kind == "eventually_above"
    # best sample is 4.0 -> margin 0.5 (> 0 satisfied)
    assert math.isclose(sp.margin(t, sig), 0.5)
    # never-reached threshold -> negative robustness
    assert spec_mod.eventually_above("x", "ch", lo=10.0).margin(t, sig) < 0.0

    bogus = spec_mod.Spec("bad", "ch", kind="does_not_exist")
    with pytest.raises(ValueError, match="unknown spec kind"):
        bogus.margin(t, sig)


# --------------------------------------------------------------------------
# sim/plots.py:55 — single-channel stacked plot
# --------------------------------------------------------------------------
def test_plot_channels_single_channel_wraps_axes():
    """fig.subplots(1, 1) returns a bare Axes, so the single-channel case
    must wrap it in a list to stay iterable; the resulting figure has one
    populated axis labelled with the channel name."""
    pytest.importorskip("matplotlib")
    from podium.sim import plots  # noqa: PLC0415 — optional extra, guarded above

    times = np.linspace(0.0, 4.0, 5)
    x_rel = np.zeros((5, 6))
    x_rel[:, 0] = np.array([10.0, 8.0, 6.0, 4.0, 2.0])
    tr = engine.Trace(times, x_rel, np.zeros((5, 6)), [], {})
    fig = plots.plot_channels(tr, channels=("range",))
    assert len(fig.axes) == 1
    ax = fig.axes[0]
    assert ax.get_ylabel() == "range"
    assert len(ax.lines) == 1
    assert np.allclose(ax.lines[0].get_ydata(), x_rel[:, 0])


# --------------------------------------------------------------------------
# sim/engine.py:132 — hemisphere continuity of the viewer quaternions
# --------------------------------------------------------------------------
def test_frame_quaternions_are_hemisphere_continuous():
    """Shepperd's method returns a raw sign per sample, so the LVLH->ECI
    quaternion sequence flips hemisphere as the frame sweeps an orbit.
    frame_quaternions() negates flipped samples so a viewer can slerp
    sample-to-sample: every adjacent pair then has a non-negative dot."""
    rv_t0 = circular_target(A, inc=0.9, raan=0.3, argp=0.5, nu=0.0)
    period = 2 * math.pi / N
    _t, x_rel, rv_t = nl.propagate_relative(
        rv_t0, np.array([10.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        tof=1.2 * period, dt=period / 60.0,
    )
    # confirm the raw (uncorrected) sequence really does flip — otherwise
    # the continuity branch would be vacuously "covered"
    raw = np.array([
        engine._quat_from_matrix(nl.lvlh_rotation(rv_t[i, 0:3], rv_t[i, 3:6]).T)
        for i in range(len(rv_t))
    ])
    assert np.any(np.sum(raw[1:] * raw[:-1], axis=1) < 0.0)

    tr = engine.Trace(_t, x_rel, rv_t, [], {})
    q = tr.frame_quaternions()
    assert np.allclose(np.linalg.norm(q, axis=1), 1.0)  # still unit
    assert np.all(np.sum(q[1:] * q[:-1], axis=1) >= 0.0)  # continuous


# --------------------------------------------------------------------------
# sim/engine.py:115 — crossing exactly on a sample
# --------------------------------------------------------------------------
def test_crossing_times_reports_exact_sample_hit():
    """When a channel equals the threshold exactly at a sample, that
    sample time is reported directly (no interpolation), distinct from the
    interpolated sign-change case."""
    times = np.array([0.0, 1.0, 2.0])
    x_rel = np.zeros((3, 6))
    x_rel[:, 0] = np.array([10.0, 5.0, 2.0])  # range hits 5.0 exactly at t=1
    tr = engine.Trace(times, x_rel, np.zeros((3, 6)), [], {})
    assert tr.crossing_times("range", 5.0) == [1.0]

    # contrast: a threshold BETWEEN samples is linearly interpolated
    got = tr.crossing_times("range", 7.5)  # between 10 and 5 -> t=0.5
    assert len(got) == 1 and math.isclose(got[0], 0.5)


# --------------------------------------------------------------------------
# dynamics/nonlinear.py:165 — argp quadrant fix for a southern periapsis
# --------------------------------------------------------------------------
def test_elements_from_rv_recovers_southern_periapsis_argp():
    """When the eccentricity vector points below the equatorial plane
    (e_vec_z < 0) the argument of perigee is in (pi, 2*pi); acos alone
    returns the wrong quadrant, so elements_from_rv reflects it. A
    round-trip through elements_to_rv recovers argp exactly."""
    argp_in = 4.2  # > pi -> periapsis in the southern hemisphere
    r, v = nl.elements_to_rv(A, 0.1, 0.6, 0.4, argp_in, 0.3, const.MU_EARTH)
    e_vec = ((float(np.dot(v, v)) - const.MU_EARTH / float(np.linalg.norm(r))) * r
             - float(np.dot(r, v)) * v) / const.MU_EARTH
    assert e_vec[2] < 0.0  # the branch precondition
    el = nl.elements_from_rv(r, v, const.MU_EARTH)
    assert math.isclose(float(el[4]), argp_in, abs_tol=1e-9)
    # and it is the reflected quadrant, not the raw acos value
    assert el[4] > math.pi


# --------------------------------------------------------------------------
# core/ya.py:83 — documents WHY the negative-theta wrap is unreachable
# --------------------------------------------------------------------------
def test_propagate_true_anomaly_is_always_in_canonical_range():
    """propagate_true_anomaly always returns a true anomaly in [0, 2*pi):
    the mean anomaly is wrapped into [0, 2*pi) and true_from_eccentric is
    monotone and non-negative there, so the belt-and-braces
    `theta1 += 2*pi` guard (ya.py line 83) is provably dead code. This
    test pins the invariant that makes it dead rather than the line."""
    rng = np.random.default_rng(7)
    for _ in range(3000):
        e = float(rng.uniform(0.0, 0.9))
        theta0 = float(rng.uniform(-2.0, 8.0))
        n = float(rng.uniform(1e-4, 1e-2))
        dt = float(rng.uniform(-30000.0, 30000.0))
        th = ya.propagate_true_anomaly(n, e, theta0, dt)
        assert 0.0 <= th < 2.0 * math.pi + 1e-12
    # sanity: at e=0 it reduces to the linear mean-anomaly advance
    assert math.isclose(ya.propagate_true_anomaly(N, 0.0, 0.0, 10.0),
                        (N * 10.0) % (2.0 * math.pi), abs_tol=1e-9)
    # and consistency with the CW limit is unaffected (e=0 STM is finite)
    assert np.all(np.isfinite(cw_mod.stm(N, 100.0)))
