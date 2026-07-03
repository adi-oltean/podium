"""Nonlinear relative-motion truth model for LEO/MEO RPOD.

Both spacecraft are integrated in ECI with identical force models (two-body
+ optional J2 + optional exponential-atmosphere drag with per-spacecraft
ballistic coefficients), and the relative state is formed by differencing
in the target LVLH frame. This is exact — no relative-dynamics
approximation — and float64 differencing error at LEO radii (~1e-9 m) is
far below truth-model needs.

LVLH convention as everywhere in Podium: x radial (zenith), y along-track,
z cross-track (orbit normal). The frame angular velocity includes the
out-of-plane component (r/h)(a_pert . z_hat) from orbit-plane precession
under perturbations — dropping it corrupts relative velocities at the
mm/s level under J2, which matters for docking-rate budgets.

This is a *sandbox-side* module (full Python allowed): it is the reference
that flight-side linearized models in ``podium.core`` are validated
against, never flight code itself.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from podium import constants as const
from podium.core import integrators, ya
from podium.core import roe as roe_mod

F64 = NDArray[np.float64]


class DensityPerturbation:
    """Seeded mean-reverting (Ornstein-Uhlenbeck) log-density factor.

    The perturbation is discretized EXACTLY (p_{k+1} = phi p_k + noise
    with phi = exp(-dt/tau)) on a fixed grid at construction and linearly
    interpolated afterwards, so density(t) is a deterministic function of
    time: RK4 propagation and bit-identical replay are preserved — the
    seed is the only source of randomness.

    Calibration: the multiplicative factor is exp(p) with stationary
    std sigma_log. The default sigma_log = 0.35 puts the +2-sigma
    excursion at exp(0.7) ~ 2.0x (+100%), inside the +50-125% band
    observed at 200-400 km in May-2024-class storms; tau defaults to
    6 h (storm-driver correlation time). Beyond the grid the last value
    holds (clamped), loudly documented rather than silently wrapped.
    """

    def __init__(
        self,
        seed: int,
        duration: float,
        dt: float = 60.0,
        sigma_log: float = 0.35,
        tau: float = 21_600.0,
    ) -> None:
        self.dt = dt
        self.sigma_log = sigma_log
        self.tau = tau
        rng = np.random.default_rng(seed)
        m = int(math.ceil(duration / dt)) + 2
        phi = math.exp(-dt / tau)
        innov = sigma_log * math.sqrt(max(0.0, 1.0 - phi * phi))
        p = np.zeros(m)
        p[0] = sigma_log * rng.standard_normal()  # stationary start
        xi = rng.standard_normal(m - 1)
        for k in range(m - 1):
            p[k + 1] = phi * p[k] + innov * xi[k]
        self._p = p

    def log_factor(self, t: float) -> float:
        u = t / self.dt
        if u <= 0.0:
            return float(self._p[0])
        i = int(u)
        if i >= len(self._p) - 1:
            return float(self._p[-1])
        frac = u - i
        return float((1.0 - frac) * self._p[i] + frac * self._p[i + 1])

    def factor(self, t: float) -> float:
        return math.exp(self.log_factor(t))


@dataclass
class DragConfig:
    """Exponential atmosphere, co-rotating with Earth.

    Defaults are representative of ~400 km altitude at moderate solar
    activity; density there varies by an order of magnitude over the solar
    cycle, so treat rho0 as a scenario parameter, not a constant. An
    optional DensityPerturbation multiplies the baseline (seeded,
    deterministic in time).
    """

    rho0: float = 3.0e-12  # density at h0 [kg/m^3]
    h0: float = 400e3  # reference altitude [m]
    scale_height: float = 60e3  # [m]
    perturbation: DensityPerturbation | None = None

    def density(self, altitude: float, t: float = 0.0) -> float:
        rho = self.rho0 * math.exp(-(altitude - self.h0) / self.scale_height)
        if self.perturbation is not None:
            rho *= self.perturbation.factor(t)
        return rho


@dataclass(frozen=True)
class ForceConfig:
    """Force-model switches for the truth propagation."""

    mu: float = const.MU_EARTH
    j2: float = 0.0  # set to const.J2_EARTH to enable
    r_body: float = const.R_EARTH
    drag: DragConfig | None = None
    omega_earth: float = const.OMEGA_EARTH


def elements_to_rv(
    a: float, e: float, inc: float, raan: float, argp: float, nu: float, mu: float
) -> tuple[F64, F64]:
    """Classical orbital elements to ECI position/velocity."""
    p = a * (1.0 - e * e)
    r_mag = p / (1.0 + e * math.cos(nu))
    r_pf = np.array([r_mag * math.cos(nu), r_mag * math.sin(nu), 0.0])
    v_pf = math.sqrt(mu / p) * np.array([-math.sin(nu), e + math.cos(nu), 0.0])

    cO, sO = math.cos(raan), math.sin(raan)
    ci, si = math.cos(inc), math.sin(inc)
    cw, sw = math.cos(argp), math.sin(argp)
    rot = np.array(
        [
            [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
            [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
            [sw * si, cw * si, ci],
        ]
    )
    return rot @ r_pf, rot @ v_pf


def elements_from_rv(r: F64, v: F64, mu: float) -> F64:
    """Osculating classical elements [a, e, inc, raan, argp, M] from ECI r, v.

    Sandbox utility (inverse of elements_to_rv) used by tests and analysis.
    Angles in radians, wrapped to [0, 2*pi). Requires a bound, non-degenerate
    orbit with e > 0 and 0 < inc < pi (argp/raan ill-conditioned otherwise —
    callers in tests keep e >= 0.01).
    """
    rn = float(np.linalg.norm(r))
    v2 = float(np.dot(v, v))
    h_vec = np.cross(r, v)
    node = np.cross(np.array([0.0, 0.0, 1.0]), h_vec)
    node_n = float(np.linalg.norm(node))
    e_vec = ((v2 - mu / rn) * r - float(np.dot(r, v)) * v) / mu
    e = float(np.linalg.norm(e_vec))
    a = 1.0 / (2.0 / rn - v2 / mu)
    inc = math.acos(float(h_vec[2]) / float(np.linalg.norm(h_vec)))
    raan = math.atan2(float(node[1]), float(node[0])) % (2.0 * math.pi)
    argp = math.acos(min(1.0, max(-1.0, float(np.dot(node, e_vec)) / (node_n * e))))
    if e_vec[2] < 0.0:
        argp = 2.0 * math.pi - argp
    nu = math.acos(min(1.0, max(-1.0, float(np.dot(e_vec, r)) / (e * rn))))
    if float(np.dot(r, v)) < 0.0:
        nu = 2.0 * math.pi - nu
    ecc_anom = 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(0.5 * nu), math.sqrt(1.0 + e) * math.cos(0.5 * nu)
    )
    mean_anom = (ecc_anom - e * math.sin(ecc_anom)) % (2.0 * math.pi)
    return np.array([a, e, inc, raan, argp, mean_anom])


def perturb_accel(r: F64, v: F64, cfg: ForceConfig, bc: float, t: float = 0.0) -> F64:
    """Perturbing (non-central) acceleration in ECI: J2 + drag."""
    a = np.zeros(3)
    rn = float(np.linalg.norm(r))
    if cfg.j2 != 0.0:
        # Standard J2 acceleration, z = Earth polar axis.
        k = -1.5 * cfg.j2 * cfg.mu * cfg.r_body**2 / rn**5
        z2r2 = (r[2] / rn) ** 2
        a[0] += k * r[0] * (1.0 - 5.0 * z2r2)
        a[1] += k * r[1] * (1.0 - 5.0 * z2r2)
        a[2] += k * r[2] * (3.0 - 5.0 * z2r2)
    if cfg.drag is not None:
        rho = cfg.drag.density(rn - cfg.r_body, t)
        # Atmosphere co-rotates: v_rel = v - omega_e x r.
        v_rel = v - np.array(
            [-cfg.omega_earth * r[1], cfg.omega_earth * r[0], 0.0]
        )
        a -= 0.5 * rho * float(np.linalg.norm(v_rel)) / bc * v_rel
    return a


def total_accel(r: F64, v: F64, cfg: ForceConfig, bc: float, t: float = 0.0) -> F64:
    rn = float(np.linalg.norm(r))
    return -cfg.mu / rn**3 * r + perturb_accel(r, v, cfg, bc, t)


def lvlh_rotation(r: F64, v: F64) -> F64:
    """Rows are the LVLH basis vectors in ECI: R @ x_eci = x_lvlh."""
    x_hat = r / np.linalg.norm(r)
    h_vec = np.cross(r, v)
    z_hat = h_vec / np.linalg.norm(h_vec)
    y_hat = np.cross(z_hat, x_hat)
    return np.array([x_hat, y_hat, z_hat])


def lvlh_omega(r: F64, v: F64, a_pert: F64) -> F64:
    """LVLH frame angular velocity, expressed in LVLH coordinates.

    omega_z = h / r^2 (orbital rate); omega_x = (r/h) a_z tracks the
    rotation of the orbit normal under out-of-plane perturbing
    acceleration a_z (LVLH z-component); omega_y = 0 identically because
    the radial direction stays in the instantaneous orbit plane.
    """
    h_vec = np.cross(r, v)
    h = float(np.linalg.norm(h_vec))
    rn = float(np.linalg.norm(r))
    a_z = float(np.dot(a_pert, h_vec / h))
    return np.array([rn / h * a_z, 0.0, h / (rn * rn)])


def eci_to_lvlh(rv_target: F64, rv_chaser: F64, cfg: ForceConfig, bc_target: float) -> F64:
    """Relative state of chaser w.r.t. target in the target LVLH frame."""
    rt, vt = rv_target[0:3], rv_target[3:6]
    rot = lvlh_rotation(rt, vt)
    omega = lvlh_omega(rt, vt, perturb_accel(rt, vt, cfg, bc_target))
    rho = rot @ (rv_chaser[0:3] - rt)
    drho = rot @ (rv_chaser[3:6] - vt) - np.cross(omega, rho)
    return np.concatenate([rho, drho])


def lvlh_to_eci(rv_target: F64, x_lvlh: F64, cfg: ForceConfig, bc_target: float) -> F64:
    """Chaser ECI state from a relative LVLH state (inverse of eci_to_lvlh)."""
    rt, vt = rv_target[0:3], rv_target[3:6]
    rot = lvlh_rotation(rt, vt)
    omega = lvlh_omega(rt, vt, perturb_accel(rt, vt, cfg, bc_target))
    rc = rt + rot.T @ x_lvlh[0:3]
    vc = vt + rot.T @ (x_lvlh[3:6] + np.cross(omega, x_lvlh[0:3]))
    return np.concatenate([rc, vc])


def _rv_from_mean_elements(el: F64, mu: float) -> F64:
    """ECI state from mean elements [a, e, inc, raan, argp, M] (two-body:
    mean == osculating)."""
    ecc = ya.kepler_eccentric(float(el[5]), float(el[1]))
    nu = ya.true_from_eccentric(ecc, float(el[1]))
    r, v = elements_to_rv(el[0], el[1], el[2], el[3], el[4], nu, mu)
    return np.concatenate([r, v])


def roe_lvlh_jacobian(chief: F64, mu: float, h: float = 1e-7) -> F64:
    """Eccentric-valid first-order ROE -> LVLH map (6x6 Jacobian).

    Central-difference Jacobian of the EXACT nonlinear chain
    roe -> deputy elements -> ECI -> LVLH about the chief's mean elements
    [a, e, inc, raan, argp, M] — valid at any eccentricity the element
    conversions support, unlike the O(e) near-circular map in
    podium.core.roe. Sandbox-side by design (the differencing is not
    static-subset); the truncation is O(|roe|^2) like any first-order
    map, quantified in the validity-envelope tests.
    """
    rv_c = _rv_from_mean_elements(chief, mu)
    cfg = ForceConfig(mu=mu)
    jac = np.zeros((6, 6))
    for k in range(6):
        dp = np.zeros(6)
        dp[k] = h
        el_p = roe_mod.elements_from_roe(chief, dp)
        el_m = roe_mod.elements_from_roe(chief, -dp)
        x_p = eci_to_lvlh(rv_c, _rv_from_mean_elements(el_p, mu), cfg, 1.0)
        x_m = eci_to_lvlh(rv_c, _rv_from_mean_elements(el_m, mu), cfg, 1.0)
        jac[:, k] = (x_p - x_m) / (2.0 * h)
    return jac


def map_roe_to_lvlh_eccentric(roe: F64, chief: F64, mu: float) -> F64:
    """First-order eccentric-valid ROE -> LVLH state via roe_lvlh_jacobian."""
    out: F64 = roe_lvlh_jacobian(chief, mu) @ roe
    return out


def _deriv(cfg: ForceConfig, bc_target: float, bc_chaser: float) -> integrators.Deriv:
    def f(t: float, y: F64) -> F64:
        out = np.empty(12)
        out[0:3] = y[3:6]
        out[3:6] = total_accel(y[0:3], y[3:6], cfg, bc_target, t)
        out[6:9] = y[9:12]
        out[9:12] = total_accel(y[6:9], y[9:12], cfg, bc_chaser, t)
        return out

    return f


def propagate_relative(
    rv_target0: F64,
    x_lvlh0: F64,
    tof: float,
    dt: float,
    cfg: ForceConfig | None = None,
    bc_target: float = 100.0,
    bc_chaser: float = 100.0,
) -> tuple[F64, F64, F64]:
    """Propagate chaser relative to target; deterministic fixed-step RK4.

    Parameters
    ----------
    rv_target0 : (6,) target ECI state at t=0.
    x_lvlh0 : (6,) chaser state relative to target, LVLH.
    tof, dt : total time and step [s]; the last step is shortened to land
        exactly on tof.
    bc_* : ballistic coefficients m/(Cd*A) [kg/m^2] (used only with drag).

    Returns
    -------
    times : (N+1,)
    x_rel : (N+1, 6) relative LVLH history (includes t=0).
    rv_target : (N+1, 6) target ECI history.
    """
    if cfg is None:
        cfg = ForceConfig()
    rv_chaser0 = lvlh_to_eci(rv_target0, x_lvlh0, cfg, bc_target)
    y = np.concatenate([rv_target0, rv_chaser0])
    f = _deriv(cfg, bc_target, bc_chaser)

    n_full = int(tof / dt)
    steps = [dt] * n_full
    rem = tof - n_full * dt
    if rem > 1e-9:
        steps.append(rem)

    times = np.zeros(len(steps) + 1)
    x_rel = np.zeros((len(steps) + 1, 6))
    rv_t = np.zeros((len(steps) + 1, 6))
    x_rel[0] = eci_to_lvlh(y[0:6], y[6:12], cfg, bc_target)
    rv_t[0] = y[0:6]

    t = 0.0
    for i, h in enumerate(steps):
        y = integrators.rk4_step(f, t, y, h)
        t += h
        times[i + 1] = t
        x_rel[i + 1] = eci_to_lvlh(y[0:6], y[6:12], cfg, bc_target)
        rv_t[i + 1] = y[0:6]
    return times, x_rel, rv_t
