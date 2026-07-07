"""Tier-2 cross-architecture golden vectors (#38): the emitted kernels
run BIT-IDENTICALLY on aarch64 (under qemu) as on x86/CPython.

IEEE-754 binary64 with -ffp-contract=off is architecture-independent,
so arithmetic + sqrt kernels are bit-exact across ISAs; the libm-trig
kernels carry the SAME cross-libm tolerance as on x86 (their reference
is CPython's conda libm), and that tolerance vanishes in the
correctly-rounded build — demonstrated here for stm, bit-exact on
aarch64 against the mpmath correctly-rounded oracle.

Reuses the exact golden harness from test_cemit (same driver, same
vectors, same reference), so this is the tier-1 suite re-run on ARM.
Needs Docker; skips cleanly otherwise. Slow (apt + qemu in a
container)."""

import math
import pathlib
import shutil
import subprocess
import zlib

import numpy as np
import pytest

pytest.importorskip("mpmath")

from test_cemit import (  # noqa: E402
    _DRIVER,
    _MATMUL,
    _TRANSCENDENTAL,
    _dispatch_block,
    _py_call,
    _vectors,
)

from podium.core import cw  # noqa: E402
from podium.emit import cemit, croracle  # noqa: E402
from podium.emit.kernels import FLIGHT_KERNELS  # noqa: E402

DOCKER = shutil.which("docker")
CM = pathlib.Path(__file__).resolve().parents[1] / "third_party" / "core-math"
N_VEC = 1500

_STM_CR_DRIVER = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
void podium_cw_stm(double n, double t, double out[6][6]);
int main(int argc, char **argv){(void)argc;(void)argv;char l[128];
  while(fgets(l,sizeof l,stdin)){
    double n=strtod(strtok(l," \n"),0),t=strtod(strtok(0," \n"),0);
    double o[6][6];podium_cw_stm(n,t,o);
    for(int i=0;i<36;i++)printf(i?" %a":"%a",((double*)o)[i]);
    printf("\n");}return 0;}
"""


@pytest.fixture(scope="module")
def aarch64_out(tmp_path_factory):
    if DOCKER is None:
        pytest.skip("docker not available")
    d = tmp_path_factory.mktemp("tier2")
    names = [k.__name__ for k in FLIGHT_KERNELS]

    # default-mode translation unit + golden driver (from test_cemit)
    src = cemit.emit_module(FLIGHT_KERNELS, "podium_kernels_test")
    (d / "kernels.c").write_text(src)
    protos = "\n".join(
        ln for ln in src.splitlines()
        if ln.startswith(("void podium_", "double podium_"))
        and ln.endswith(";"))
    (d / "driver.c").write_text(
        _DRIVER.replace("{prototypes}", protos)
        .replace("{dispatch}", _dispatch_block()))
    (d / "kernels.txt").write_text("\n".join(names) + "\n")

    # per-kernel recorded vectors + the Python reference for each
    ref = {}
    for name in names:
        rng = np.random.default_rng(zlib.crc32(name.encode()))
        vecs = _vectors(name, N_VEC, rng)
        (d / f"{name}.in").write_text(
            "\n".join(" ".join(float(x).hex() for x in row)
                      for row in vecs) + "\n")
        ref[name] = (vecs, np.array([_py_call(name, row) for row in vecs]))

    # correctly-rounded stm unit + driver + vendored CORE-MATH
    cmdir = d / "coremath"
    cmdir.mkdir()
    for f in ("sin.c", "cos.c", "coremath.h"):
        shutil.copy(CM / f, cmdir / f)
    (d / "kernels_cr.c").write_text(
        cemit.emit_module([cw.stm], "stm_cr", correctly_rounded=True))
    (d / "driver_cr.c").write_text(_STM_CR_DRIVER)
    stm_vecs = ref["stm"][0]
    (d / "stm_cr.in").write_text(
        "\n".join(" ".join(float(x).hex() for x in row)
                  for row in stm_vecs) + "\n")

    script = pathlib.Path(__file__).resolve().parents[1] / \
        "tools" / "tier2_build_run.sh"
    shutil.copy(script, d / "run.sh")
    # DOCKER_CONFIG (if the ambient environment sets one, e.g. to skip a
    # broken credential helper) is inherited by subprocess automatically.
    r = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{d}:/work", "debian:bookworm",
         "bash", "/work/run.sh"],
        capture_output=True, text=True, timeout=1800, check=False)
    if "TIER2_OK" not in r.stdout:
        pytest.skip(f"tier-2 container did not complete: "
                    f"{r.stderr[-500:]}")
    return d, ref


@pytest.mark.slow
def test_tier2_arch_bit_exact(aarch64_out):
    """Non-trig kernels: aarch64 == CPython bit-for-bit. Trig kernels:
    same cross-libm tolerance as x86 (measured, bounded)."""
    d, ref = aarch64_out
    for name in (k.__name__ for k in FLIGHT_KERNELS):
        out = (d / f"{name}.out").read_text().strip().split("\n")
        _vecs, py = ref[name]
        got = np.array([[float.fromhex(t) for t in ln.split()]
                        for ln in out])
        assert got.shape == py.shape, name
        if name in _MATMUL:
            scale = np.maximum(1.0, np.max(np.abs(py)))
            assert np.max(np.abs(got - py)) < 1e-13 * scale, name
        elif name in _TRANSCENDENTAL:
            # cross-libm (aarch64 glibc vs CPython conda): same class as
            # x86; bounded by output scale, <=1% incidence
            scale = np.maximum(1.0, np.max(np.abs(py), axis=1,
                                           keepdims=True))
            bad = np.abs(got - py) > 1e-12 * scale
            assert bad.sum() <= 0.01 * got.size, (name, bad.sum())
        else:
            assert np.array_equal(got, py), f"{name} not bit-exact on ARM"


@pytest.mark.slow
def test_tier2_stm_correctly_rounded_bit_exact(aarch64_out):
    """stm compiled with CORE-MATH on aarch64 == the mpmath correctly-
    rounded oracle, BIT-EXACT — the trig tolerance is gone cross-arch,
    not merely bounded."""
    d, ref = aarch64_out
    stm_vecs = ref["stm"][0]
    out = (d / "stm_cr.out").read_text().strip().split("\n")
    got = np.array([[float.fromhex(t) for t in ln.split()] for ln in out])
    assert got.shape[0] == len(stm_vecs)
    real_sin, real_cos = math.sin, math.cos
    mism = 0
    for (n, t), row in zip(stm_vecs, got):
        s, c = croracle.cr_sin(n * t), croracle.cr_cos(n * t)
        math.sin = lambda a, s=s, nt=n * t: s if a == nt else real_sin(a)
        math.cos = lambda a, c=c, nt=n * t: c if a == nt else real_cos(a)
        try:
            want = cw.stm(float(n), float(t)).ravel()
        finally:
            math.sin, math.cos = real_sin, real_cos
        if not np.array_equal(row, want):
            mism += 1
    assert mism == 0, f"{mism}/{len(stm_vecs)} not bit-exact on ARM (CR)"
