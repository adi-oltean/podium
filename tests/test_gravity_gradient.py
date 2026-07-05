"""Gravity-gradient torque receipts (#44): the dominant LEO
environmental torque, validated against exact analytic results.

Equilibrium (torque zero at nadir-aligned principal axes) and the
restoring sign are pure algebra. The dynamic locking and the classic
pitch-libration frequency omega = n*sqrt(3(I_R - I_Y)/I_P) are checked
in an orbit-coupled propagation, and the gravity-gradient stability
boundary is exhibited (the unstable inertia ordering diverges).

Axis convention: body 1 = roll (along-track), 2 = pitch (orbit
normal), 3 = yaw (nadir). Stable gravity-gradient config has the
minimum moment along nadir (I_R > I_Y).
"""

import numpy as np
import pytest

from podium.core import quat
from podium.dynamics import attitude as att
from podium.sim.engine import _quat_from_matrix

N = 0.0011  # orbital mean motion [rad/s]

# LVLH-aligned equilibrium: body roll->along-track, pitch->normal,
# yaw->nadir. Columns of R are the body axes expressed in inertial at
# t=0 (nadir=(-1,0,0), along-track=(0,1,0), normal=(0,0,1)).
_R_EQ = np.array([[0.0, 0.0, -1.0],
                  [1.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0]])


def _q_eq():
    return _quat_from_matrix(_R_EQ.T)


def _nadir_inertial(t):
    return np.array([-np.cos(N * t), -np.sin(N * t), 0.0])


def _propagate(inertia, q, w, steps, dt=1.0):
    """Orbit-coupled gravity-gradient propagation; returns the signed
    pitch-libration series (angle of the body nadir axis from true
    nadir, about the orbit normal)."""
    pitch = []
    for k in range(steps):
        nadir_in = _nadir_inertial(k * dt)
        nb = quat.rotate(quat.conjugate(q), nadir_in)
        tau = att.gravity_gradient_torque(nb, inertia, N)
        q, w = att.step(q, w, inertia, tau, dt)
        b3 = quat.rotate(q, np.array([0.0, 0.0, 1.0]))
        pitch.append(float(np.arcsin(np.clip(
            np.cross(nadir_in, b3)[2], -1.0, 1.0))))
    return np.array(pitch)


def test_torque_zero_at_nadir_aligned_equilibrium():
    """When a principal axis points at nadir the torque is exactly
    zero — the gravity-gradient equilibrium."""
    inertia = np.diag([180.0, 140.0, 100.0])
    for axis in np.eye(3):
        tau = att.gravity_gradient_torque(axis, inertia, N)
        assert np.allclose(tau, 0.0, atol=1e-15), axis


def test_pitch_torque_structure():
    """A pitch tilt (nadir moved by eps in the roll-yaw plane) yields a
    torque purely about the pitch axis of exactly
    3 n^2 (I_R - I_Y) sin(eps) cos(eps) — the term whose sign sets
    gravity-gradient stability (the dynamic restoring is shown by the
    libration and instability tests below)."""
    i_r, i_y = 180.0, 100.0
    inertia = np.diag([i_r, 140.0, i_y])
    eps = 0.01
    o = np.array([np.sin(eps), 0.0, np.cos(eps)])
    tau = att.gravity_gradient_torque(o, inertia, N)
    assert abs(tau[0]) < 1e-15 and abs(tau[2]) < 1e-15  # pure pitch
    expect = 3 * N**2 * (i_r - i_y) * np.sin(eps) * np.cos(eps)
    assert abs(tau[1] - expect) < 1e-12


def test_equilibrium_locks_onto_nadir():
    """From the LVLH-aligned state, co-rotating at the orbital rate,
    the body stays nadir-locked to <0.1 deg over a full orbit — the
    torque + integrator reproduce gravity-gradient capture."""
    inertia = np.diag([180.0, 140.0, 100.0])
    q = _q_eq()
    w = np.array([0.0, -N, 0.0])           # co-rotate with the orbit
    steps = int(2 * np.pi / N)
    pitch = _propagate(inertia, q, w, steps)
    assert np.max(np.abs(pitch)) < np.radians(0.1)


@pytest.mark.slow
def test_pitch_libration_frequency_matches_analytic():
    """A tiny pitch-rate kick launches a small libration whose measured
    frequency matches omega = n*sqrt(3(I_R - I_Y)/I_P) — the classic
    gravity-gradient result — to better than 1%."""
    i_r, i_p, i_y = 180.0, 140.0, 100.0
    inertia = np.diag([i_r, i_p, i_y])
    q = _q_eq()
    w = np.array([0.0, -N, 0.0]) + np.array([0.0, 1.0e-5, 0.0])
    dt = 1.0
    pitch = _propagate(inertia, q, w, int(3 * 2 * np.pi / N), dt)
    s = pitch - pitch.mean()
    crossings = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    assert len(crossings) >= 3
    w_meas = 2 * np.pi / (np.mean(np.diff(crossings)) * dt)
    w_analytic = N * np.sqrt(3 * (i_r - i_y) / i_p)
    assert abs(w_meas - w_analytic) < 0.01 * w_analytic


@pytest.mark.slow
def test_unstable_ordering_diverges():
    """Gravity-gradient stability boundary: with the nadir moment
    LARGEST (I_Y > I_R) the equilibrium is unstable and a tiny kick
    diverges past 30 deg, whereas the stable ordering stays small."""
    unstable = np.diag([100.0, 140.0, 180.0])   # nadir (I3) largest
    q = _q_eq()
    w = np.array([0.0, -N, 0.0]) + np.array([0.0, 1.0e-5, 0.0])
    pitch = _propagate(unstable, q, w, int(2 * 2 * np.pi / N))
    assert np.max(np.abs(pitch)) > np.radians(30.0)
