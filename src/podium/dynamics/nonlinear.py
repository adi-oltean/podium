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
from podium.core import integrators

F64 = NDArray[np.float64]


@dataclass(frozen=True)
class DragConfig:
    """Exponential atmosphere, co-rotating with Earth.

    Defaults are representative of ~400 km altitude at moderate solar
    activity; density there varies by an order of magnitude over the solar
    cycle, so treat rho0 as a scenario parameter, not a constant.
    """

    rho0: float = 3.0e-12  # density at h0 [kg/m^3]
    h0: float = 400e3  # reference altitude [m]
    scale_height: float = 60e3  # [m]

    def density(self, altitude: float) -> float:
        return self.rho0 * math.exp(-(altitude - self.h0) / self.scale_height)


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


def perturb_accel(r: F64, v: F64, cfg: ForceConfig, bc: float) -> F64:
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
        rho = cfg.drag.density(rn - cfg.r_body)
        # Atmosphere co-rotates: v_rel = v - omega_e x r.
        v_rel = v - np.array(
            [-cfg.omega_earth * r[1], cfg.omega_earth * r[0], 0.0]
        )
        a -= 0.5 * rho * float(np.linalg.norm(v_rel)) / bc * v_rel
    return a


def total_accel(r: F64, v: F64, cfg: ForceConfig, bc: float) -> F64:
    rn = float(np.linalg.norm(r))
    return -cfg.mu / rn**3 * r + perturb_accel(r, v, cfg, bc)


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


def _deriv(cfg: ForceConfig, bc_target: float, bc_chaser: float) -> integrators.Deriv:
    def f(_t: float, y: F64) -> F64:
        out = np.empty(12)
        out[0:3] = y[3:6]
        out[3:6] = total_accel(y[0:3], y[3:6], cfg, bc_target)
        out[6:9] = y[9:12]
        out[9:12] = total_accel(y[6:9], y[9:12], cfg, bc_chaser)
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
