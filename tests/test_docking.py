"""Docking pulse-control laws: closed-loop click-dynamics simulations.

The plant mirrors the SpaceX ISS-sim control model: per axis, velocity
changes ONLY in fixed increments (one per 'click'), control replans at
1 Hz, dynamics are a pure integrator (drift-free over these timescales).
Success gates are the sim's: |lateral| < 0.2 m, |angle| < 0.2 deg, closing
rate < 0.2 m/s at contact.
"""

import pytest

from podium.control import docking

DV_TRANS = 0.05  # m/s per click (sim fine-mode increment scale)
DV_ROT = 0.05  # deg/s per click


def run_axial(range0, dt=1.0, dv=DV_TRANS):
    """Closing axis: range decreases at rate v; control tracks the profile."""
    r, v, t, clicks = range0, 0.0, 0.0, 0
    while r > 0.0 and t < 3600.0:
        v_cmd = docking.approach_rate(r)
        n = docking.pulses_needed(v_cmd, v, dv)
        v += n * dv
        clicks += abs(n)
        r -= v * dt
        t += dt
    return r, v, t, clicks


def run_lateral(pos0, v0=0.0, tmax=400.0, dt=1.0, dv=DV_TRANS):
    p, v, t, clicks = pos0, v0, 0.0, 0
    while t < tmax:
        v_cmd = docking.lateral_rate_cmd(p)
        n = docking.pulses_needed(v_cmd, v, dv)
        v += n * dv
        clicks += abs(n)
        p += v * dt
        t += dt
    return p, v, clicks


@pytest.mark.parametrize("range0", [20.0, 75.0, 200.0])
def test_axial_contacts_within_rate_gate(range0):
    r, v, t, _ = run_axial(range0)
    assert r <= 0.0  # reached the port
    assert 0.0 < v < 0.2  # closing, under the 0.2 m/s gate
    assert t < 1800.0  # and in reasonable time


def test_axial_rate_tracks_profile_far_field():
    # From 200 m the commanded profile immediately ramps toward the ceiling;
    # after 60 s the vehicle should be closing at multiple pulse increments.
    r, v = 200.0, 0.0
    for _ in range(60):
        n = docking.pulses_needed(docking.approach_rate(r), v, DV_TRANS)
        v += n * DV_TRANS
        r -= v
    assert v > 1.0
    assert r < 130.0


@pytest.mark.parametrize(("pos0", "v0"), [(12.0, 0.0), (-8.0, 0.2), (0.5, -0.1)])
def test_lateral_nulls_within_gate(pos0, v0):
    p, v, _ = run_lateral(pos0, v0)
    assert abs(p) < 0.2  # sim lateral success gate
    assert abs(v) < DV_TRANS  # residual below one pulse


@pytest.mark.parametrize("ang0", [24.9, -12.0, 1.5])
def test_attitude_nulls_within_gate(ang0):
    a, w = ang0, 0.0
    for _ in range(400):
        w_cmd = docking.angle_rate_cmd(a)
        n = docking.pulses_needed(w_cmd, w, DV_ROT)
        w += n * DV_ROT
        a += w
    assert abs(a) < 0.2  # sim angular success gate
    assert abs(w) <= DV_ROT + 1e-12


def test_no_overshoot_through_port_side():
    # Lateral convergence must not ring: once inside 0.2 m it stays inside.
    p, v = 12.0, 0.0
    inside = False
    for _ in range(400):
        v_cmd = docking.lateral_rate_cmd(p)
        v += docking.pulses_needed(v_cmd, v, DV_TRANS) * DV_TRANS
        p += v
        if inside:
            assert abs(p) < 0.25
        if abs(p) < 0.15 and abs(v) <= DV_TRANS:
            inside = True
    assert inside


def test_pulses_deadband_prevents_chatter():
    # Within half a pulse of the command, no clicks are issued.
    assert docking.pulses_needed(0.10, 0.08, 0.05) == 0
    assert docking.pulses_needed(0.10, 0.0, 0.05) == 2
    assert docking.pulses_needed(-0.10, 0.0, 0.05) == -2
    assert docking.pulses_needed(0.0, 0.26, 0.05) == -5


def test_approach_rate_profile_shape():
    assert docking.approach_rate(0.0) == docking.CONTACT_RATE
    assert docking.approach_rate(5.0) == docking.CONTACT_RATE
    assert docking.approach_rate(60.0) == pytest.approx(1.0)
    assert docking.approach_rate(1e5) == docking.FAR_RATE
