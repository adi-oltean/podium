"""ROE module receipts.

The J2 STM is pinned entrywise by a central-difference Jacobian of the
exact secular flow; the Keplerian STM and near-circular map/control
matrix are validated against the nonlinear ECI truth model; safety
metrics against brute-force separation scans.
"""

import math

import numpy as np
import pytest

from podium import constants as const
from podium.core import cw, integrators, roe, ya
from podium.dynamics import nonlinear as nl
from podium.guidance import safety

MU = const.MU_EARTH
J2 = const.J2_EARTH
RE = const.R_EARTH


def secular_rates(a, e, inc, mu, j2, re):
    """Exact secular J2 rates (raan_dot, argp_dot, m_dot) — textbook."""
    n = math.sqrt(mu / a**3)
    eta = math.sqrt(1 - e * e)
    kappa = 0.75 * j2 * math.sqrt(mu) * re * re / (a**3.5 * eta**4)
    ci = math.cos(inc)
    raan_dot = -2.0 * kappa * ci
    argp_dot = kappa * (5.0 * ci * ci - 1.0)
    m_dot = n + kappa * eta * (3.0 * ci * ci - 1.0)
    return raan_dot, argp_dot, m_dot


def flow_mean_j2(el, dt, mu=MU, j2=J2, re=RE):
    """Exact secular flow of mean elements under Keplerian + J2."""
    raan_dot, argp_dot, m_dot = secular_rates(el[0], el[1], el[2], mu, j2, re)
    out = el.copy()
    out[3] = el[3] + raan_dot * dt
    out[4] = el[4] + argp_dot * dt
    out[5] = el[5] + m_dot * dt
    return out


CHIEF = np.array([7_100_000.0, 0.08, math.radians(52.0), 0.7, 1.1, 0.4])


def test_roe_elements_roundtrip():
    r0 = np.array([2e-4, -8e-4, 1e-4, -2e-4, 5e-5, 1.2e-4])
    deputy = roe.elements_from_roe(CHIEF, r0)
    back = roe.roe_from_elements(CHIEF, deputy)
    assert np.allclose(back, r0, rtol=0, atol=1e-15)


def test_j2_stm_matches_flow_jacobian():
    """Central-difference Jacobian of the exact secular flow == stm_j2.

    This pins every entry of the closed-form matrix; tolerance reflects
    O(h^2) differencing error against entries scaled by dt."""
    dt = 20_000.0
    phi = roe.stm_j2(MU, J2, RE, CHIEF[0], CHIEF[1], CHIEF[2], CHIEF[4], dt)
    h = 1e-8
    jac = np.zeros((6, 6))
    for k in range(6):
        dp = np.zeros(6)
        dp[k] = h
        el_p = roe.elements_from_roe(CHIEF, dp)
        el_m = roe.elements_from_roe(CHIEF, -dp)
        roe_p = roe.roe_from_elements(flow_mean_j2(CHIEF, dt), flow_mean_j2(el_p, dt))
        roe_m = roe.roe_from_elements(flow_mean_j2(CHIEF, dt), flow_mean_j2(el_m, dt))
        jac[:, k] = (roe_p - roe_m) / (2 * h)
    assert np.allclose(phi, jac, rtol=2e-5, atol=2e-5)


def test_j2_stm_reduces_to_keplerian():
    dt = 5_000.0
    phi = roe.stm_j2(MU, 0.0, RE, CHIEF[0], CHIEF[1], CHIEF[2], CHIEF[4], dt)
    n = math.sqrt(MU / CHIEF[0] ** 3)
    assert np.allclose(phi, roe.stm_keplerian(n, dt), atol=1e-15)


@pytest.mark.slow
def test_keplerian_stm_vs_truth():
    """Two-body ECI truth: osculating == mean, so the ROE history from the
    truth model checks the Keplerian STM directly. The only STM
    approximation is n_d - n ~ -1.5 n da, so the residual must scale
    quadratically with da."""
    n = math.sqrt(MU / CHIEF[0] ** 3)
    errs = []
    for scale in (1.0, 0.5):
        r0 = scale * np.array([3e-5, -2e-4, 8e-5, -6e-5, 4e-5, 7e-5])
        dep = roe.elements_from_roe(CHIEF, r0)
        nu_c = _true_from_mean(CHIEF[5], CHIEF[1])
        nu_d = _true_from_mean(dep[5], dep[1])
        rv_c = np.concatenate(nl.elements_to_rv(CHIEF[0], CHIEF[1], CHIEF[2],
                                                CHIEF[3], CHIEF[4], nu_c, MU))
        rv_d = np.concatenate(nl.elements_to_rv(dep[0], dep[1], dep[2],
                                                dep[3], dep[4], nu_d, MU))
        tof = 3 * 2 * math.pi / n
        y = np.concatenate([rv_c, rv_d])
        f = nl._deriv(nl.ForceConfig(), 100.0, 100.0)
        steps = 3000
        h = tof / steps
        for i in range(steps):
            y = integrators.rk4_step(f, i * h, y, h)
        el_c = nl.elements_from_rv(y[0:3], y[3:6], MU)
        el_d = nl.elements_from_rv(y[6:9], y[9:12], MU)
        roe_truth = roe.roe_from_elements(el_c, el_d)
        roe_stm = roe.stm_keplerian(n, tof) @ r0
        errs.append(np.linalg.norm(roe_truth - roe_stm))
        # absolute: residual is the second-order drift term ~ 3.4*n*t*da^2
        assert errs[-1] < 6.0 * n * tof * (scale * 3e-5) ** 2 + 1e-9
    assert errs[0] / errs[1] > 3.0  # quadratic in the ROE scale


def _true_from_mean(m, e):
    ecc = ya.kepler_eccentric(m, e)
    return ya.true_from_eccentric(ecc, e)


def test_map_matches_cw_near_circular():
    """Near-circular chief: map -> Keplerian STM -> map must agree with the
    CW propagation of the mapped initial state."""
    a = 6_900_000.0
    n = math.sqrt(MU / a**3)
    u0 = 0.9
    r0 = np.array([1e-5, -6e-5, 2e-5, -1.5e-5, 1.2e-5, 2.5e-5])
    x0 = roe.map_roe_to_lvlh(r0, a, n, u0)
    for tof in (500.0, 3000.0):
        u1 = u0 + n * tof
        x_roe = roe.map_roe_to_lvlh(roe.stm_keplerian(n, tof) @ r0, a, n, u1)
        x_cw = cw.stm(n, tof) @ x0
        # identical linear models expressed in different coordinates
        assert np.allclose(x_roe[:3], x_cw[:3], atol=1e-6)
        assert np.allclose(x_roe[3:], x_cw[3:], atol=1e-9)


def test_map_roundtrip():
    a, n, u = 6_900_000.0, 1.1e-3, 2.2
    r0 = np.array([3e-5, -1e-4, 5e-5, 2e-5, -4e-5, 6e-5])
    x = roe.map_roe_to_lvlh(r0, a, n, u)
    assert np.allclose(roe.map_lvlh_to_roe(x, a, n, u), r0, atol=1e-18)


def test_control_matrix_vs_truth_impulses():
    """Finite-difference impulse receipts: apply dv along each LVLH axis in
    the truth model (instantaneous velocity change; no propagation) and
    compare the osculating-ROE jump against Gamma. Chief e = 0.01, so the
    near-circular Gamma is accurate to O(e) ~ 1%."""
    chief = np.array([6_900_000.0, 0.01, math.radians(51.0), 0.4, 0.7, 1.9])
    a = chief[0]
    n = math.sqrt(MU / a**3)
    nu = _true_from_mean(chief[5], chief[1])
    r_eci, v_eci = nl.elements_to_rv(chief[0], chief[1], chief[2],
                                     chief[3], chief[4], nu, MU)
    rot = nl.lvlh_rotation(r_eci, v_eci)  # rows: R, T-ish, N
    el0 = nl.elements_from_rv(r_eci, v_eci, MU)
    # true argument of latitude enters Gamma; mean vs true differs at O(e)
    u = chief[4] + nu
    gamma = roe.control_matrix(a, n, u)
    dv = 0.02  # m/s: small enough that O(dv^2) << O(e) tolerance
    for axis in range(3):
        v_new = v_eci + rot[axis] * dv
        el1 = nl.elements_from_rv(r_eci, v_new, MU)
        droe = roe.roe_from_elements(el0, el1)
        expected = gamma[:, axis] * dv
        scalefree = np.abs(droe - expected) * (n * a / dv)
        # entries are O(1) scale-free; near-circular truncation is O(e)=0.01,
        # allow 3x margin
        assert np.all(scalefree < 0.03), f"axis {axis}: {scalefree}"


@pytest.mark.slow
def test_j2_stm_vs_truth_secular_drift():
    """15-orbit ECI truth with J2: compare secular trends (mean-longitude
    drift and e-vector rotation) via ROE built from per-orbit means of
    osculating elements. Orbit-averaging suppresses short-period J2 terms
    to O(J2^2); tolerance dominated by the residual averaging error."""
    chief = np.array([7_000_000.0, 0.05, math.radians(50.0), 0.3, 0.8, 0.0])
    r0 = np.array([0.0, 2e-4, 1e-4, -5e-5, 8e-5, -6e-5])
    dep = roe.elements_from_roe(chief, r0)
    n = math.sqrt(MU / chief[0] ** 3)
    period = 2 * math.pi / n
    orbits = 15
    cfg = nl.ForceConfig(j2=J2)

    nu_c = _true_from_mean(chief[5], chief[1])
    nu_d = _true_from_mean(dep[5], dep[1])
    rv_c = np.concatenate(nl.elements_to_rv(chief[0], chief[1], chief[2],
                                            chief[3], chief[4], nu_c, MU))
    rv_d = np.concatenate(nl.elements_to_rv(dep[0], dep[1], dep[2],
                                            dep[3], dep[4], nu_d, MU))
    f = nl._deriv(cfg, 100.0, 100.0)
    y = np.concatenate([rv_c, rv_d])
    dt = 10.0
    steps_per_orbit = int(period / dt)

    def orbit_mean_roe(y):
        """Average osculating ROE over one orbit, advancing y in place."""
        acc = np.zeros(6)
        nonlocal_y = y
        for i in range(steps_per_orbit):
            nonlocal_y = integrators.rk4_step(f, 0.0, nonlocal_y, dt)
            el_c = nl.elements_from_rv(nonlocal_y[0:3], nonlocal_y[3:6], MU)
            el_d = nl.elements_from_rv(nonlocal_y[6:9], nonlocal_y[9:12], MU)
            acc += roe.roe_from_elements(el_c, el_d)
        return acc / steps_per_orbit, nonlocal_y

    roe_first, y = orbit_mean_roe(y)
    for _ in range(orbits - 2):
        _, y = orbit_mean_roe(y)
    roe_last, y = orbit_mean_roe(y)

    span = (orbits - 1) * period
    phi = roe.stm_j2(MU, J2, RE, chief[0], chief[1], chief[2], chief[4], span)
    roe_pred = phi @ roe_first
    # dl secular drift and e-vector rotation are the observable trends;
    # osculating-vs-mean and averaging residuals are O(J2 * (Re/a)^2) ~ 1e-3
    # of the elements themselves — compare drifts, not absolutes.
    drift_truth = roe_last - roe_first
    drift_pred = roe_pred - roe_first
    assert abs(drift_pred[1] - drift_truth[1]) < 0.12 * abs(drift_truth[1])
    rot_truth = math.atan2(roe_last[3], roe_last[2]) - math.atan2(roe_first[3], roe_first[2])
    rot_pred = math.atan2(roe_pred[3], roe_pred[2]) - math.atan2(roe_first[3], roe_first[2])
    assert abs(rot_pred - rot_truth) < 0.25 * abs(rot_truth) + 2e-4


def test_ei_separation_safety_ordering():
    a = 7_000_000.0
    aligned = np.array([0.0, 0.0, 2e-5, 0.0, 2e-5, 0.0])       # de || di
    perp = np.array([0.0, 0.0, 2e-5, 0.0, 0.0, 2e-5])          # de ⊥ di
    assert safety.ei_separation_angle(aligned) < 1e-9
    assert abs(safety.ei_separation_angle(perp) - math.pi / 2) < 1e-9
    # parallel vectors keep the RN trajectory off the origin; perpendicular
    # phasing lets radial and cross-track zeros coincide
    assert safety.min_rn_separation(aligned, a) > 100.0
    assert safety.min_rn_separation(perp, a) < 10.0


def test_min_rn_separation_matches_bruteforce():
    a = 7_000_000.0
    r = np.array([1e-5, 0.0, 3e-5, -1e-5, 2e-5, 2.5e-5])
    n = 1.1e-3
    best = math.inf
    for k in range(3600):
        u = 2 * math.pi * k / 3600
        x = roe.map_roe_to_lvlh(r, a, n, u)
        best = min(best, math.hypot(x[0], x[2]))
    assert abs(safety.min_rn_separation(r, a) - best) < 0.02 * best + 0.01


def test_rn_margin_sign():
    a = 7_000_000.0
    r = np.array([0.0, 0.0, 4e-5, 0.0, 4e-5, 0.0])  # aligned, ~280 m RN
    assert safety.rn_margin(r, a, 200.0) > 0.0
    assert safety.rn_margin(r, a, 400.0) < 0.0


# --- #6 follow-ups: drag STM, eccentric map, analytic screening --------


def _int_pow(a0, rate, m, t):
    """Integral of (a0 + rate*s)^m ds over [0, t], series in xi = rate*t/a0.

    The naive closed form ((a0+rt)^(m+1) - a0^(m+1))/(r(m+1)) is
    catastrophically cancelling for the tiny rates the FD check uses (the
    dadot-sensitive part sits below the subtraction's rounding noise);
    the series is exact to ~xi^4 ~ 1e-24 here."""
    xi = rate * t / a0
    return a0**m * t * (
        1.0
        + m * xi / 2.0
        + m * (m - 1.0) * xi * xi / 6.0
        + m * (m - 1.0) * (m - 2.0) * xi**3 / 24.0
    )


def flow_mean_j2_drag(el, dadot, a_chief, dt, mu=MU, j2=J2, re=RE):
    """Exact secular flow of deputy mean elements under Keplerian + J2
    with a(s) = a0 + a_chief*dadot*s (density-model-free differential
    drag: constant relative sma decay, e and inc held). Closed-form time
    integrals — the oracle for the augmented STM."""
    a0, e, inc = el[0], el[1], el[2]
    rate = a_chief * dadot
    eta = math.sqrt(1.0 - e * e)
    kfac = 0.75 * j2 * math.sqrt(mu) * re * re / eta**4
    ci = math.cos(inc)
    i_m35 = _int_pow(a0, rate, -3.5, dt)
    i_m15 = _int_pow(a0, rate, -1.5, dt)
    out = el.copy()
    out[0] = a0 + rate * dt
    out[3] = el[3] - 2.0 * kfac * ci * i_m35
    out[4] = el[4] + kfac * (5.0 * ci * ci - 1.0) * i_m35
    out[5] = el[5] + math.sqrt(mu) * i_m15 + kfac * eta * (3.0 * ci * ci - 1.0) * i_m35
    return out


def test_j2_drag_stm_matches_flow_jacobian():
    """Central-difference Jacobian of the exact augmented flow pins every
    entry of the 7x7 density-model-free J2+drag STM."""
    dt = 20_000.0
    phi = roe.stm_j2_drag(MU, J2, RE, CHIEF[0], CHIEF[1], CHIEF[2], CHIEF[4], dt)
    h7 = np.array([1e-8] * 6 + [1e-10])  # dadot in 1/s: scale h to the units
    jac = np.zeros((7, 7))
    for k in range(7):
        dp = np.zeros(7)
        dp[k] = h7[k]
        el_p = roe.elements_from_roe(CHIEF, dp[:6])
        el_m = roe.elements_from_roe(CHIEF, -dp[:6])
        chief_t = flow_mean_j2_drag(CHIEF, 0.0, CHIEF[0], dt)
        dep_p = flow_mean_j2_drag(el_p, dp[6], CHIEF[0], dt)
        dep_m = flow_mean_j2_drag(el_m, -dp[6], CHIEF[0], dt)
        roe_p = np.concatenate([roe.roe_from_elements(chief_t, dep_p), [dp[6]]])
        roe_m = np.concatenate([roe.roe_from_elements(chief_t, dep_m), [-dp[6]]])
        jac[:, k] = (roe_p - roe_m) / (2.0 * h7[k])
    assert np.allclose(phi, jac, rtol=3e-5, atol=3e-5)


def test_j2_drag_stm_reduces_to_j2():
    dt = 8_000.0
    phi7 = roe.stm_j2_drag(MU, J2, RE, CHIEF[0], CHIEF[1], CHIEF[2], CHIEF[4], dt)
    phi6 = roe.stm_j2(MU, J2, RE, CHIEF[0], CHIEF[1], CHIEF[2], CHIEF[4], dt)
    assert np.allclose(phi7[:6, :6], phi6, atol=1e-18)
    assert phi7[0, 6] == dt
    assert np.allclose(phi7[6, :6], 0.0)


@pytest.mark.slow
def test_j2_drag_stm_vs_truth_differential_drag():
    """Dual-ECI truth with differential ballistic coefficients: extract
    dadot from the first orbit, then the augmented STM must predict the
    quadratic mean-longitude runaway over 8 orbits. Tolerance reflects
    drag decay not being exactly linear and orbit-averaging residue."""
    a = 6_778_137.0
    n = math.sqrt(MU / a**3)
    period = 2 * math.pi / n
    orbits = 8
    cfg = nl.ForceConfig(j2=J2, drag=nl.DragConfig(rho0=1e-11))
    rv0 = np.concatenate(nl.elements_to_rv(a, 0.001, 0.9, 0.5, 1.2, 0.3, MU))
    times, x_rel, rv_t = nl.propagate_relative(
        rv0, np.zeros(6), orbits * period, dt=5.0, cfg=cfg,
        bc_target=200.0, bc_chaser=40.0,
    )
    # orbit-averaged osculating ROE per orbit (co-located start: roe0 = 0)
    steps = len(times) - 1
    per_orbit = steps // orbits
    roe_hist = np.zeros((orbits, 6))
    for k in range(orbits):
        sl = slice(k * per_orbit, (k + 1) * per_orbit)
        acc = np.zeros(6)
        for i in range(sl.start, sl.stop):
            el_c = nl.elements_from_rv(rv_t[i, 0:3], rv_t[i, 3:6], MU)
            rv_c_full = nl.lvlh_to_eci(rv_t[i], x_rel[i], cfg, 200.0)
            el_d = nl.elements_from_rv(rv_c_full[0:3], rv_c_full[3:6], MU)
            acc += roe.roe_from_elements(el_c, el_d)
        roe_hist[k] = acc / per_orbit
    # dadot from the first-orbit da slope
    dadot = (roe_hist[1, 0] - roe_hist[0, 0]) / period
    t_span = (orbits - 1) * period
    phi = roe.stm_j2_drag(MU, J2, RE, a, 0.001, 0.9, 1.2, t_span)
    roe7_0 = np.concatenate([roe_hist[0], [dadot]])
    pred = phi @ roe7_0
    dl_truth = roe_hist[-1, 1] - roe_hist[0, 1]
    dl_pred = pred[1] - roe_hist[0, 1]
    assert abs(dl_truth) * a > 100.0  # the effect is macroscopic (>100 m)
    assert abs(dl_pred - dl_truth) < 0.25 * abs(dl_truth)
    # and the da prediction tracks the decay trend
    assert abs(pred[0] - roe_hist[-1, 0]) < 0.35 * abs(roe_hist[-1, 0] - roe_hist[0, 0])


def test_eccentric_map_reduces_to_near_circular():
    a = 6_900_000.0
    n = math.sqrt(MU / a**3)
    chief = np.array([a, 1e-8, math.radians(51.0), 0.4, 0.7, 1.3])
    r0 = np.array([1e-5, -6e-5, 2e-5, -1.5e-5, 1.2e-5, 2.5e-5])
    u = chief[4] + chief[5]  # = argp + M ~ argp + nu at e ~ 0
    x_ecc = nl.map_roe_to_lvlh_eccentric(r0, chief, MU)
    x_circ = roe.map_roe_to_lvlh(r0, a, n, u)
    # near-circular map is itself O(e)-accurate; at e=1e-8 they must agree
    assert np.allclose(x_ecc[:3], x_circ[:3], atol=2e-3)
    assert np.allclose(x_ecc[3:], x_circ[3:], atol=1e-6)


def test_eccentric_map_first_order_exactness():
    """The Jacobian map's error against the exact nonlinear separation
    must scale quadratically with ROE size — first-order in ROE at any
    eccentricity (here e = 0.3, far outside the near-circular map)."""
    chief = np.array([7_200_000.0, 0.3, math.radians(60.0), 1.0, 0.8, 2.1])
    r_base = np.array([2e-5, -8e-5, 4e-5, -3e-5, 2e-5, 5e-5])
    errs = []
    for scale in (1.0, 0.5):
        r0 = scale * r_base
        dep = roe.elements_from_roe(chief, r0)
        x_exact = nl.eci_to_lvlh(
            nl._rv_from_mean_elements(chief, MU),
            nl._rv_from_mean_elements(dep, MU),
            nl.ForceConfig(mu=MU), 1.0,
        )
        x_map = nl.map_roe_to_lvlh_eccentric(r0, chief, MU)
        errs.append(float(np.linalg.norm(x_map[:3] - x_exact[:3])))
    assert errs[0] / errs[1] > 3.0  # quadratic in the ROE scale
    sep = float(np.linalg.norm(x_exact[:3]))
    assert errs[1] < 0.02 * sep  # small in absolute terms too


def test_analytic_min_rn_matches_dense_scan():
    """The quartic-root minimum equals a 200k-point scan to high
    precision across seeded random geometries — the scan is the oracle,
    the analytic form is the product."""
    rng = np.random.default_rng(123)
    a = 7_000_000.0
    us = np.linspace(0.0, 2 * math.pi, 200_000, endpoint=False)
    cu, su = np.cos(us), np.sin(us)
    for _ in range(60):
        r = rng.normal(0.0, 3e-5, 6)
        x = r[0] - r[2] * cu - r[3] * su
        z = r[4] * su - r[5] * cu
        dense = a * float(np.sqrt(np.min(x * x + z * z)))
        analytic = safety.min_rn_separation_analytic(r, a)
        assert analytic <= dense + 1e-6  # true min is <= any sampled value
        assert abs(analytic - dense) < 1e-6 * a * 3e-5 + 1e-4 * max(dense, 1e-9)
    # degenerate geometries
    assert safety.min_rn_separation_analytic(np.zeros(6), a) == 0.0
    circ = np.array([2e-5, 0.0, 0.0, 0.0, 0.0, 0.0])  # pure radial offset
    assert abs(safety.min_rn_separation_analytic(circ, a) - a * 2e-5) < 1e-6
