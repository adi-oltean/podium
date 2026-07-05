"""Rigid-body rotational dynamics (truth model).

State: quaternion q (scalar-first, body -> reference) + body angular
velocity w [rad/s]. Euler's equations with a full inertia tensor and an
external body-frame torque; RK4 with quaternion renormalization per step
(the receipt suite pins energy and angular-momentum conservation, which
is where unnormalized quaternions and sloppy integrators show up).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from podium.core import quat

F64 = NDArray[np.float64]


def euler_deriv(w: F64, inertia: F64, torque_body: F64) -> F64:
    """Body angular acceleration: I w' = tau - w x (I w)."""
    out: F64 = np.linalg.solve(inertia, torque_body - np.cross(w, inertia @ w))
    return out


def step(q: F64, w: F64, inertia: F64, torque_body: F64, dt: float) -> tuple[F64, F64]:
    """One RK4 step of the coupled kinematics + Euler equations."""
    def f(qq: F64, ww: F64) -> tuple[F64, F64]:
        return quat.deriv(qq, ww), euler_deriv(ww, inertia, torque_body)

    k1q, k1w = f(q, w)
    k2q, k2w = f(q + 0.5 * dt * k1q, w + 0.5 * dt * k1w)
    k3q, k3w = f(q + 0.5 * dt * k2q, w + 0.5 * dt * k2w)
    k4q, k4w = f(q + dt * k3q, w + dt * k3w)
    q_out = quat.normalize(q + (dt / 6.0) * (k1q + 2 * k2q + 2 * k3q + k4q))
    w_out = w + (dt / 6.0) * (k1w + 2 * k2w + 2 * k3w + k4w)
    return q_out, w_out


def gravity_gradient_torque(nadir_body: F64, inertia: F64, n: float) -> F64:
    """Gravity-gradient torque in the body frame for a circular orbit,

        tau = 3 n^2 (o_hat x I o_hat),

    where o_hat is the unit NADIR direction (toward the primary's
    centre) expressed in body coordinates and n is the orbital mean
    motion. Vanishes when a principal axis points at nadir (the
    gravity-gradient equilibrium); to first order it restores the body
    toward the local-vertical/local-horizontal orientation, the basis
    of passive gravity-gradient stabilization."""
    o = np.asarray(nadir_body, dtype=np.float64)
    o = o / np.linalg.norm(o)
    tau: F64 = 3.0 * n * n * np.cross(o, inertia @ o)
    return tau


def aerodynamic_torque(v_rel_body: F64, rho: float, cd_area: float,
                       r_cp: F64) -> F64:
    """Aerodynamic disturbance torque in the body frame,

        F_drag = -1/2 rho (Cd A) |v_rel| v_rel,   tau = r_cp x F_drag,

    where v_rel_body is the atmosphere-relative velocity in body
    coordinates, cd_area = Cd*A, and r_cp is the center-of-pressure
    offset from the center of mass (body frame). This is the second
    dominant LEO attitude disturbance after gravity gradient. With the
    center of pressure BEHIND the center of mass (downstream) the torque
    weathervanes the body toward the flow — passive aerodynamic
    stabilization. The force is consistent with the truth model's drag
    (podium.dynamics.nonlinear): F = m*a_drag with cd_area = m/bc."""
    v = np.asarray(v_rel_body, dtype=np.float64)
    f_drag: F64 = -0.5 * rho * cd_area * float(np.linalg.norm(v)) * v
    tau: F64 = np.cross(np.asarray(r_cp, dtype=np.float64), f_drag)
    return tau


def srp_torque(sun_dir_body: F64, area: float, cr: float, r_cp: F64,
               pressure: float = 4.5606e-6, illuminated: bool = True) -> F64:
    """Solar-radiation-pressure disturbance torque in the body frame,

        F_srp = -P Cr A s_hat,   tau = r_cp x F_srp,

    where s_hat is the unit direction TO the Sun in body coordinates
    (the force pushes away from the Sun, hence the minus), P is the SRP
    constant at the spacecraft's heliocentric distance (default 1 AU),
    Cr in [1, 2] the reflectivity (1 absorbing, 2 perfect specular),
    A the projected illuminated area, and r_cp the center-of-pressure
    offset from the c.m. In eclipse (illuminated=False) the torque is
    zero. SRP is the dominant environmental attitude disturbance at GEO;
    negligible but nonzero in LEO except in Earth's shadow."""
    if not illuminated:
        return np.zeros(3)
    s = np.asarray(sun_dir_body, dtype=np.float64)
    s = s / np.linalg.norm(s)
    f_srp: F64 = -pressure * cr * area * s
    tau: F64 = np.cross(np.asarray(r_cp, dtype=np.float64), f_srp)
    return tau


def dipole_field(r_eci: F64, dipole_moment: float = 7.94e22,
                 axis: F64 | None = None) -> F64:
    """Geomagnetic field [T] of a centered dipole at position r_eci,

        B = (mu0 m / 4 pi r^3) [3 (m_hat . r_hat) r_hat - m_hat],

    with mu0/4pi = 1e-7. `axis` is the dipole unit vector in ECI
    (default -z, i.e. the geomagnetic south pole near the north
    geographic pole, aligned-dipole approximation). Magnitude at the
    equatorial surface is ~3.1e-5 T and falls as 1/r^3."""
    r = np.asarray(r_eci, dtype=np.float64)
    rn = float(np.linalg.norm(r))
    r_hat = r / rn
    m_hat = (np.array([0.0, 0.0, -1.0]) if axis is None
             else np.asarray(axis, dtype=np.float64)
             / np.linalg.norm(axis))
    coeff = 1.0e-7 * dipole_moment / rn**3
    b: F64 = coeff * (3.0 * float(np.dot(m_hat, r_hat)) * r_hat - m_hat)
    return b


def magnetic_torque(dipole_body: F64, b_body: F64) -> F64:
    """Magnetic disturbance torque tau = m x B, where m is the
    spacecraft's residual (or commanded magnetorquer) magnetic dipole
    [A m^2] and B the local geomagnetic field [T], both in the body
    frame. The fourth classic environmental attitude disturbance; also
    the actuation model for magnetic torquers."""
    tau: F64 = np.cross(np.asarray(dipole_body, dtype=np.float64),
                        np.asarray(b_body, dtype=np.float64))
    return tau


def kinetic_energy(w: F64, inertia: F64) -> float:
    return 0.5 * float(w @ (inertia @ w))


def momentum_inertial(q: F64, w: F64, inertia: F64) -> F64:
    """Angular momentum rotated into the reference frame (conserved when
    the external torque is zero)."""
    return quat.rotate(q, inertia @ w)
