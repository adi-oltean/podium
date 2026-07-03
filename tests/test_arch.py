"""ARCH rendezvous benchmark receipts (python side).

The reachability proof itself runs in the Julia gate (tools/reach); here
we pin the model to Podium's kernel, sanity-check stability, validate the
export schema, and run falsification-style spec checks from the initial
set — necessary conditions that catch model transcription errors before
the expensive proof ever runs.
"""

import json
import math

import numpy as np

from podium.control import lqr
from podium.core import cw
from podium.guidance import arch


def test_abort_mode_is_planar_cw():
    """The abort matrix must equal Podium's CW linearization at the GEO
    mean motion implied by its own 2n entry — ties the benchmark to
    podium.core.cw rather than trusting transcribed constants."""
    n = arch.N_RAD_MIN
    # build the planar sub-block of cw_deriv as a matrix (per-minute units)
    a_cw = np.zeros((4, 4))
    for j, e in enumerate(np.eye(4)):
        full = np.array([e[0], e[1], 0.0, e[2], e[3], 0.0])
        d = cw.cw_deriv(full, n, np.zeros(3))
        a_cw[:, j] = [d[0], d[1], d[3], d[4]]
    assert np.allclose(arch.A_ABORT[:4, :4], a_cw, atol=1e-12)


def test_controlled_modes_are_stable():
    for a in (arch.A_APPROACH, arch.A_ATTEMPT):
        eig = np.linalg.eigvals(a[:4, :4])
        assert np.all(eig.real < 0.0)


def test_export_schema_and_roundtrip():
    model = arch.export_model(abort_time=120.0)
    js = json.loads(json.dumps(model))
    assert [m["name"] for m in js["modes"]] == list(arch.MODE_NAMES)
    assert len(js["transitions"]) == 3
    assert js["initial"]["center"][0] == -900.0
    a2 = np.array(js["modes"][1]["A"])
    assert np.allclose(a2, arch.A_ATTEMPT)
    # every halfspace has a 5-vector and a bound
    for m in js["modes"]:
        for h in m["invariant"]:
            assert len(h["a"]) == 5
    # no-abort variant has two modes only
    assert len(arch.export_model(-1.0)["modes"]) == 2


def test_simulation_specs_srna01():
    """No-abort scenario: every corner of the initial set must reach the
    attempt mode and satisfy LOS + velocity throughout — a necessary
    condition for the set-based proof."""
    for x0 in arch.initial_corners():
        times, states, modes = arch.simulate(x0, abort_time=-1.0)
        assert modes[-1] == 2  # reached and stayed in attempt mode
        m = arch.spec_margins(states, modes)
        assert m["los_cone"] > 0.0
        assert m["velocity"] > 0.0
        # progress: ends close to the origin
        assert np.hypot(states[-1, 0], states[-1, 1]) < 5.0


def test_simulation_specs_sra01():
    """Abort at t=120 min: all three properties must hold from every
    corner, and the abort mode must actually be entered."""
    for x0 in arch.initial_corners():
        times, states, modes = arch.simulate(x0, abort_time=120.0)
        assert modes[-1] == 3
        m = arch.spec_margins(states, modes)
        assert m["los_cone"] > 0.0
        assert m["velocity"] > 0.0
        assert m["abort_avoidance"] > 0.0


def test_clock_is_time():
    _, states, _ = arch.simulate(arch.X0_CENTER, dt=0.02, horizon=10.0)
    assert math.isclose(states[-1, 4], 10.0, rel_tol=1e-9)


# --- Podium-synthesized controller variant (#11) -----------------------


def test_care_solves_riccati():
    """Residual of A'P + PA - P B R^-1 B' P + Q at machine precision."""
    a = arch.cw_planar(arch.N_RAD_MIN)
    for q in (arch.Q_APPROACH, arch.Q_ATTEMPT):
        p = lqr.care(a, arch._B_ACCEL, q, arch.R_CTRL)
        res = (
            a.T @ p + p @ a
            - p @ arch._B_ACCEL @ np.linalg.inv(arch.R_CTRL) @ arch._B_ACCEL.T @ p
            + q
        )
        assert np.max(np.abs(res)) < 1e-10 * max(1.0, float(np.max(np.abs(q))))
        assert np.allclose(p, p.T)
        assert np.all(np.linalg.eigvalsh(p) > 0.0)


def test_podium_gains_are_genuinely_ours():
    k1, k2 = arch.podium_gains()
    r1, r2 = arch.implied_reference_gains()
    # same regime (the Q/R were chosen for that), but different numbers
    assert not np.allclose(k1, r1, rtol=1e-3)
    assert not np.allclose(k2, r2, rtol=1e-3)
    # closed loops Hurwitz
    a1, a2, _ = arch.podium_mode_matrices()
    for a in (a1, a2):
        assert np.all(np.linalg.eigvals(a[:4, :4]).real < 0.0)


def test_podium_variant_simulation_specs():
    """All three properties from every initial corner, both scenarios —
    the necessary condition before the reachability proof."""
    for abort in (-1.0, 120.0):
        for x0 in arch.initial_corners():
            _, states, modes = arch.simulate(x0, abort_time=abort,
                                             gains="podium")
            m = arch.spec_margins(states, modes)
            assert m["los_cone"] > 0.0, (abort, x0[:2], m)
            assert m["velocity"] > 0.0, (abort, x0[:2], m)
            if abort >= 0.0:
                assert modes[-1] == 3
                assert m["abort_avoidance"] > 0.0, (abort, x0[:2], m)
            else:
                assert modes[-1] == 2


def test_podium_variant_export():
    js = arch.export_model(120.0, gains="podium")
    assert "podium" in str(js["name"])
    a1_ref = np.array(arch.export_model(120.0)["modes"][0]["A"])  # type: ignore[index]
    a1_pod = np.array(js["modes"][0]["A"])  # type: ignore[index]
    assert not np.allclose(a1_ref, a1_pod)
