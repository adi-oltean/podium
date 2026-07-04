"""cFS integration receipt: the generated relative-nav app, compiled
against the emitted kernels, must reproduce the Python EKF reference on
a recorded measurement stream. Same init as the app; same %a transport
as the golden vectors, so the comparison is exact up to the documented
matmul-order class already covered in test_cemit."""

import shutil
import subprocess

import numpy as np
import pytest

from podium.core import cw
from podium.emit import cemit, cfsapp
from podium.emit.kernels import FLIGHT_KERNELS
from podium.nav import ekf

GCC = shutil.which("gcc")

# must match PODIUM_NAV_Init in the generated app
N = 0.0011313666536110223
DT = 1.0
Q_ACCEL = 2.0e-6
R_VAR = 0.01  # EVA-proven envelope floor ((10 cm)^2); matches the app


def _python_reference(meas):
    """The app's EKF loop in Python: update_sequential then predict."""
    x = np.zeros(6)
    P = np.diag([100.0**2] * 3 + [1.0] * 3).astype(float)
    phi = cw.stm(N, DT)
    Q = ekf.process_noise_wna(DT, Q_ACCEL)
    out = []
    for z in meas:
        x, P = ekf.update_sequential(x, P, np.asarray(z), R_VAR)
        x, P = ekf.predict(x, P, phi, Q)
        out.append(x.copy())
    return np.array(out)


@pytest.fixture(scope="module")
def app_exe(tmp_path_factory):
    if GCC is None:
        pytest.skip("gcc not available")
    d = tmp_path_factory.mktemp("cfsapp")
    cfsapp.generate(str(d), cemit.emit_module(FLIGHT_KERNELS))
    exe = d / "nav_app"
    subprocess.run(
        [GCC, "-std=c99", "-O2", "-ffp-contract=off", "-Wall", "-Werror",
         "-I", str(d), str(d / "podium_nav_app.c"), str(d / "cfe_shim.c"),
         str(d / "podium_kernels.c"), "-lm", "-o", str(exe)],
        check=True, capture_output=True, text=True)
    return exe


def test_cfs_app_matches_python_ekf(app_exe):
    rng = np.random.default_rng(36)
    # a closing V-bar approach with docking-grade measurement noise
    truth = np.linspace([0.0, -50.0, 0.0], [0.0, -2.0, 0.0], 60)
    meas = truth + rng.normal(0.0, 0.05, truth.shape)
    lines = "\n".join(" ".join(float(v).hex() for v in row)
                      for row in meas) + "\n"
    r = subprocess.run([str(app_exe)], input=lines, capture_output=True,
                       text=True, check=True)
    got = np.array([[float.fromhex(t) for t in ln.split()]
                    for ln in r.stdout.strip().split("\n")])
    ref = _python_reference(meas)
    assert got.shape == ref.shape == (60, 6)
    # bit-exact but for the emitted predict's BLAS-vs-naive matmul order
    # (the documented tolerance class from test_cemit); scale by state
    scale = np.maximum(1.0, np.max(np.abs(ref), axis=1, keepdims=True))
    assert np.max(np.abs(got - ref) / scale) < 1e-11
    # and the filter actually converged toward the closing target
    assert abs(got[-1, 1] - (-2.0)) < 1.0


def test_cfs_app_emits_evs_and_runs(app_exe):
    """The app announces init on EVS (stderr) and processes the stream
    to completion — the cFS run loop and SB plumbing work end to end."""
    r = subprocess.run([str(app_exe)],
                       input="0x0p+0 -0x1p+5 0x0p+0\n",
                       capture_output=True, text=True, check=True)
    assert "PODIUM_NAV initialized" in r.stderr
    assert len(r.stdout.strip().split("\n")) == 1
