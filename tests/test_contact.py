"""MuJoCo contact/capture receipts: nominal capture, the IDSS-box
tie-in (acceptance box conditions physically capture), funnel misses
beyond the mouth, bounce at excessive rate, force-vs-rate scaling,
determinism, and the capture envelope."""

import numpy as np
import pytest

pytest.importorskip("mujoco")

from podium.sim import contact, idss  # noqa: E402


@pytest.mark.slow
def test_centered_slow_approach_captures():
    o = contact.simulate_contact(closing_rate=0.08)
    assert o.captured and not o.bounced
    assert o.peak_force == 0.0  # straight into the throat, no wall


@pytest.mark.slow
def test_idss_box_corners_capture():
    """The tie-in receipt: contact conditions at the IDSS translation
    box corners physically capture in the contact sim (20 N docking
    thrust — probe-drogue practice keeps pushing through capture)."""
    box = idss.IdssBox()
    for closing in (box.closing_min, box.closing_max):
        for off in (0.0, box.lateral_offset_max):
            for lrate in (0.0, box.lateral_rate_max):
                o = contact.simulate_contact(
                    closing_rate=closing, lateral_offset=off,
                    lateral_rate=lrate, thrust=20.0)
                assert o.captured, (closing, off, lrate, o)
    # the nominal corner also captures fully ballistically
    o = contact.simulate_contact(closing_rate=0.10, lateral_offset=0.10,
                                 lateral_rate=0.04)
    assert o.captured


@pytest.mark.slow
def test_miss_beyond_mouth_and_bounce_at_speed():
    miss = contact.simulate_contact(closing_rate=0.08, lateral_offset=0.45,
                                    thrust=20.0)
    assert not miss.captured
    fast = contact.simulate_contact(closing_rate=2.0)
    assert not fast.captured
    assert fast.bounced  # entered and came back out


@pytest.mark.slow
def test_peak_force_grows_with_closing_rate():
    forces = []
    for cr in (0.1, 0.3, 0.6, 1.2):
        o = contact.simulate_contact(closing_rate=cr, lateral_offset=0.15)
        forces.append(o.peak_force)
    assert all(b > a for a, b in zip(forces, forces[1:])), forces


@pytest.mark.slow
def test_deterministic():
    a = contact.simulate_contact(closing_rate=0.1, lateral_offset=0.12,
                                 lateral_rate=0.03)
    b = contact.simulate_contact(closing_rate=0.1, lateral_offset=0.12,
                                 lateral_rate=0.03)
    assert (a.captured, a.max_tip_x, a.seat_time, a.peak_force) \
        == (b.captured, b.max_tip_x, b.seat_time, b.peak_force)


@pytest.mark.slow
def test_capture_envelope_boundary():
    """Envelope over lateral offset at nominal closing rate: everything
    comfortably inside the mouth captures, everything beyond misses,
    and the boundary sits between."""
    table = contact.capture_envelope(
        offsets=np.array([0.0, 0.1, 0.2, 0.35, 0.45]),
        closing_rates=np.array([0.08]),
    )
    by_off = {row["offset"]: row["captured"] for row in table}
    assert by_off[0.0] and by_off[0.1] and by_off[0.2]
    assert not by_off[0.45]