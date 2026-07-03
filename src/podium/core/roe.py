"""Relative orbital elements: state, maps, Koenig STMs, control input.

Quasi-nonsingular ROE (D'Amico convention) of a deputy relative to a chief
with mean elements (a, e, inc, raan, argp, M), u = argp + M:

    da  = (a_d - a) / a                 relative semi-major axis
    dl  = (u_d - u) + (raan_d - raan) cos(inc)   relative mean longitude
    dex = e_d cos(argp_d) - e cos(argp)
    dey = e_d sin(argp_d) - e sin(argp)
    dix = inc_d - inc
    diy = (raan_d - raan) sin(inc)

State ordering: [da, dl, dex, dey, dix, diy] (dimensionless; multiply by a
for meters). The STMs propagate MEAN elements: secular Keplerian and J2
flows only. Osculating-mean conversion is a truth-side concern; receipts
in tests/test_roe.py validate the STMs as linearizations of the exact
secular flow (finite-difference Jacobian check pins every entry) and
against the nonlinear ECI truth model.

References: D'Amico & Montenbruck, JGCD 29(3), 2006 (e/i-vector
separation); Koenig, Guffanti & D'Amico, JGCD 40(7), 2017,
doi:10.2514/1.G002409 (closed-form STMs); Gaias & Lovera, JGCD 44, 2021
(ROE <-> Cartesian bridge). The near-circular LVLH map and impulsive
control matrix follow the standard first-order forms and are validated
empirically against the truth model (finite-difference impulses).

Static-subset compliant: fixed shapes, closed-form only, bounded loops.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

from podium.verify import Interval, contract

F64 = NDArray[np.float64]

_TWO_PI = 2.0 * math.pi


def _wrap_pi(x: float) -> float:
    """Wrap angle to (-pi, pi]."""
    return math.atan2(math.sin(x), math.cos(x))


def roe_from_elements(chief: F64, deputy: F64) -> F64:
    """Quasi-nonsingular ROE of deputy w.r.t. chief.

    Both element sets are [a, e, inc, raan, argp, M] (mean elements, rad).
    Angle differences are wrapped to (-pi, pi], so separations must be
    well inside half an orbit along-track.
    """
    a_c, e_c, i_c = chief[0], chief[1], chief[2]
    a_d, e_d = deputy[0], deputy[1]
    du = _wrap_pi((deputy[4] + deputy[5]) - (chief[4] + chief[5]))
    draan = _wrap_pi(deputy[3] - chief[3])
    out = np.empty(6)
    out[0] = (a_d - a_c) / a_c
    out[1] = du + draan * math.cos(i_c)
    out[2] = e_d * math.cos(deputy[4]) - e_c * math.cos(chief[4])
    out[3] = e_d * math.sin(deputy[4]) - e_c * math.sin(chief[4])
    out[4] = deputy[2] - i_c
    out[5] = draan * math.sin(i_c)
    return out


def elements_from_roe(chief: F64, roe: F64) -> F64:
    """Deputy mean elements [a, e, inc, raan, argp, M] from chief + ROE.

    Exact algebraic inverse of roe_from_elements (no linearization).
    """
    a_c, e_c, i_c = chief[0], chief[1], chief[2]
    ex = e_c * math.cos(chief[4]) + roe[2]
    ey = e_c * math.sin(chief[4]) + roe[3]
    e_d = math.sqrt(ex * ex + ey * ey)
    argp_d = math.atan2(ey, ex)
    draan = roe[5] / math.sin(i_c)
    u_d = (chief[4] + chief[5]) + roe[1] - draan * math.cos(i_c)
    out = np.empty(6)
    out[0] = a_c * (1.0 + roe[0])
    out[1] = e_d
    out[2] = i_c + roe[4]
    out[3] = chief[3] + draan
    out[4] = argp_d
    out[5] = u_d - argp_d
    return out


@contract(n=Interval(1e-5, 1e-2))
def stm_keplerian(n: float, dt: float) -> F64:
    """Keplerian ROE STM: only the mean longitude drifts, dl' = -1.5 n da."""
    phi = np.eye(6)
    phi[1, 0] = -1.5 * n * dt
    return phi


@contract(e=Interval(0.0, 0.9), inc=Interval(1e-3, 3.14))
def stm_j2(
    mu: float, j2: float, r_body: float, a: float, e: float, inc: float, argp: float, dt: float
) -> F64:
    """Koenig closed-form J2 (+ Keplerian) ROE STM, quasi-nonsingular.

    Linearizes the secular J2 flow (standard rates: raan' = -2*kappa*cos i,
    argp' = kappa*Q, M' = n + kappa*eta*P) about the chief; valid for
    arbitrary eccentricity below the contract bound. Every entry is pinned
    in tests by a central-difference Jacobian of the exact secular flow.
    """
    n = math.sqrt(mu / (a * a * a))
    eta = math.sqrt(1.0 - e * e)
    kappa = 0.75 * j2 * math.sqrt(mu) * r_body * r_body / (
        math.pow(a, 3.5) * math.pow(eta, 4.0)
    )
    ci = math.cos(inc)
    p_i = 3.0 * ci * ci - 1.0
    q_i = 5.0 * ci * ci - 1.0
    s_i = math.sin(2.0 * inc)
    t_i = math.sin(inc) ** 2
    e_f = 1.0 + eta
    f_f = 4.0 + 3.0 * eta
    g_f = 1.0 / (eta * eta)

    exi = e * math.cos(argp)
    eyi = e * math.sin(argp)
    argp_dot = kappa * q_i
    argp_f = argp + argp_dot * dt
    exf = e * math.cos(argp_f)
    eyf = e * math.sin(argp_f)
    cwt = math.cos(argp_dot * dt)
    swt = math.sin(argp_dot * dt)

    phi = np.eye(6)
    # relative mean longitude row
    phi[1, 0] = -(1.5 * n + 3.5 * kappa * e_f * p_i) * dt
    phi[1, 2] = kappa * exi * f_f * g_f * p_i * dt
    phi[1, 3] = kappa * eyi * f_f * g_f * p_i * dt
    phi[1, 4] = -kappa * f_f * s_i * dt
    # relative eccentricity vector rows (rotation by argp drift + couplings)
    phi[2, 0] = 3.5 * kappa * eyf * q_i * dt
    phi[2, 2] = cwt - 4.0 * kappa * exi * eyf * g_f * q_i * dt
    phi[2, 3] = -swt - 4.0 * kappa * eyi * eyf * g_f * q_i * dt
    phi[2, 4] = 5.0 * kappa * eyf * s_i * dt
    phi[3, 0] = -3.5 * kappa * exf * q_i * dt
    phi[3, 2] = swt + 4.0 * kappa * exi * exf * g_f * q_i * dt
    phi[3, 3] = cwt + 4.0 * kappa * eyi * exf * g_f * q_i * dt
    phi[3, 4] = -5.0 * kappa * exf * s_i * dt
    # relative inclination vector rows (dix constant; diy couples)
    phi[5, 0] = 3.5 * kappa * s_i * dt
    phi[5, 2] = -4.0 * kappa * exi * g_f * s_i * dt
    phi[5, 3] = -4.0 * kappa * eyi * g_f * s_i * dt
    phi[5, 4] = 2.0 * kappa * t_i * dt
    return phi


@contract(e=Interval(0.0, 0.9), inc=Interval(1e-3, 3.14))
def stm_j2_drag(
    mu: float, j2: float, r_body: float, a: float, e: float, inc: float, argp: float, dt: float
) -> F64:
    """Koenig density-model-free J2+drag STM on the augmented 7-state
    [da, dl, dex, dey, dix, diy, dadot], dadot = d(da)/dt.

    Differential drag is modeled as a CONSTANT relative semi-major-axis
    decay rate (estimated by the filter, not derived from a density
    model — hence density-model-free; near-circular differential-drag
    assumption: e-vector decay neglected). da(t) = da0 + dadot*t, so
    every secular rate that couples to da accumulates dadot with weight
    t^2/2 instead of t: the new column is the J2 STM's da column scaled
    by dt/2, plus the identity coupling da <- dadot*dt. Pinned entrywise
    in tests by an FD Jacobian of the exact augmented flow (closed-form
    time integrals of the rates along a(t)).
    """
    phi = np.eye(7)
    phi6 = stm_j2(mu, j2, r_body, a, e, inc, argp, dt)
    phi[0:6, 0:6] = phi6
    phi[0, 6] = dt
    for i in range(1, 6):
        phi[i, 6] = phi6[i, 0] * (0.5 * dt)
    return phi


@contract(n=Interval(1e-5, 1e-2))
def map_roe_to_lvlh(roe: F64, a: float, n: float, u: float) -> F64:
    """First-order near-circular map: ROE -> LVLH state [m, m/s].

    LVLH as everywhere in Podium: x radial, y along-track, z cross-track.
    Valid to O(e) and O(|roe|^2); u is the chief mean argument of latitude.
    """
    su, cu = math.sin(u), math.cos(u)
    out = np.empty(6)
    out[0] = a * (roe[0] - roe[2] * cu - roe[3] * su)
    out[1] = a * (roe[1] + 2.0 * roe[2] * su - 2.0 * roe[3] * cu)
    out[2] = a * (roe[4] * su - roe[5] * cu)
    out[3] = a * n * (roe[2] * su - roe[3] * cu)
    out[4] = a * n * (-1.5 * roe[0] + 2.0 * roe[2] * cu + 2.0 * roe[3] * su)
    out[5] = a * n * (roe[4] * cu + roe[5] * su)
    return out


@contract(n=Interval(1e-5, 1e-2))
def map_lvlh_to_roe(x: F64, a: float, n: float, u: float) -> F64:
    """Exact linear inverse of map_roe_to_lvlh."""
    su, cu = math.sin(u), math.cos(u)
    xa = x[0] / a
    ya = x[1] / a
    za = x[2] / a
    vxa = x[3] / (a * n)
    vya = x[4] / (a * n)
    vza = x[5] / (a * n)
    out = np.empty(6)
    out[0] = 4.0 * xa + 2.0 * vya
    out[1] = ya - 2.0 * vxa
    c1 = out[0] - xa  # = dex*cos u + dey*sin u
    out[2] = c1 * cu + vxa * su
    out[3] = c1 * su - vxa * cu
    out[4] = za * su + vza * cu
    out[5] = -za * cu + vza * su
    return out


@contract(n=Interval(1e-5, 1e-2))
def control_matrix(a: float, n: float, u: float) -> F64:
    """Impulsive control-input matrix Gamma (6x3), near-circular.

    d(roe) = Gamma @ dv with dv = [dv_R, dv_T, dv_N] in m/s (radial,
    along-track, cross-track at chief argument of latitude u). The
    N-coupling into dl cancels exactly between du and draan*cos(inc) —
    a known feature of the mean-longitude choice.
    """
    su, cu = math.sin(u), math.cos(u)
    g = np.zeros((6, 3))
    s = 1.0 / (n * a)
    g[0, 1] = 2.0 * s
    g[1, 0] = -2.0 * s
    g[2, 0] = su * s
    g[2, 1] = 2.0 * cu * s
    g[3, 0] = -cu * s
    g[3, 1] = 2.0 * su * s
    g[4, 2] = cu * s
    g[5, 2] = su * s
    return g
