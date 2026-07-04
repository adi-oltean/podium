"""Reference-mission receipts: the full stack in one run — PTR/CTCS
replanning on EKF estimates through imperfect actuators, corridor-gated
terminal feedback into the IDSS box, MuJoCo capture, attitude hold, the
exact-rational abort certificate, and a byte-deterministic audit
bundle. Dispersed campaign: every mission captures."""

import math

import numpy as np
import pytest

pytest.importorskip("cvxpy")
pytest.importorskip("mujoco")

from podium.sim import mission  # noqa: E402


@pytest.fixture(scope="module")
def ref():
    return mission.fly(seed=7)


@pytest.mark.slow
def test_reference_mission_docks_and_captures(ref):
    assert ref.captured
    assert 2400.0 < ref.contact_time < 4200.0
    for k, v in ref.idss_translation.items():
        assert v > 0.0, (k, v)
    for k, v in ref.idss_rotation.items():
        assert v > 0.0, (k, v)
    for k, v in ref.spec_margins.items():
        assert v > 0.0, (k, v)
    # measured ~12.4 m/s: 2 km transfer + continuous terminal feedback
    # chatter (noise-fed rate commands each tick); a deadband would trim
    # it — recorded as polish, not hidden by a loose bound
    assert ref.dv_total < 15.0
    assert ref.barrier_ok  # abort-safety certificate verified exactly


@pytest.mark.slow
def test_audit_bundle_deterministic(ref):
    b1 = mission.audit_bundle(ref, 7)
    res2 = mission.fly(seed=7)
    b2 = mission.audit_bundle(res2, 7)
    assert b1 == b2  # byte-for-byte
    assert '"captured": true' in b1
    assert '"verified_exact_rational": true' in b1


@pytest.mark.slow
def test_dispersed_campaign_all_capture():
    for seed in (1, 2, 3):
        res = mission.fly(seed=seed, dispersed=True)
        assert res.captured, seed
        assert res.idss_translation["lateral_offset"] > 0.0, seed
        assert math.isfinite(res.contact_time)


def test_initial_states_stay_in_certified_set():
    """Dispersed starts must remain inside the barrier-certified
    ellipsoid (the abort-safety fact must actually apply to them)."""
    c = np.array([float(v) for v in mission.SAFE_CASE.center])
    r = np.array([float(v) for v in mission.SAFE_CASE.radii])
    for seed in range(20):
        x = mission._initial_state(np.random.default_rng(seed))
        u = np.concatenate([x[0:3], x[3:6] / mission.N_REF])
        s = np.sum(((u - c) / r) ** 2)
        assert s <= 1.0 + 1e-9, (seed, s)