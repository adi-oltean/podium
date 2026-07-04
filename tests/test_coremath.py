"""CORE-MATH correctly-rounded receipts (#37): the last tier-1
tolerance retired.

1. The vendored cr_sin/cr_cos, compiled, equal the mpmath correctly-
   rounded oracle BIT-EXACT over the kernels' argument ranges — proving
   both sides are correctly rounded and agree.
2. The stm kernel emitted in correctly-rounded mode is BIT-EXACT
   against a stm reference that uses the oracle for sin/cos — the
   measured 21/72000 ULP exception (documented in test_cemit) is gone,
   not merely bounded.
"""

import math
import pathlib
import shutil
import subprocess

import numpy as np
import pytest

pytest.importorskip("mpmath")

from podium.core import cw  # noqa: E402
from podium.emit import cemit, croracle  # noqa: E402

GCC = shutil.which("gcc")
CM = pathlib.Path(__file__).resolve().parents[1] / "third_party" / "core-math"

_FN_DRIVER = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "coremath.h"
int main(void) {
    char line[128];
    while (fgets(line, sizeof line, stdin)) {
        double x = strtod(line, NULL);
        printf("%a %a\n", cr_sin(x), cr_cos(x));
    }
    return 0;
}
"""


@pytest.fixture(scope="module")
def cr_fn_exe(tmp_path_factory):
    if GCC is None:
        pytest.skip("gcc not available")
    d = tmp_path_factory.mktemp("crfn")
    (d / "drv.c").write_text(_FN_DRIVER)
    exe = d / "crfn"
    subprocess.run(
        [GCC, "-std=c99", "-O2", "-ffp-contract=off", "-I", str(CM),
         str(d / "drv.c"), str(CM / "sin.c"), str(CM / "cos.c"),
         "-lm", "-o", str(exe)],
        check=True, capture_output=True, text=True)
    return exe


@pytest.mark.slow
def test_coremath_functions_match_oracle(cr_fn_exe):
    """cr_sin/cr_cos (compiled) == mpmath correctly-rounded oracle, bit
    for bit, over the trig arguments the kernels actually produce
    (n*t up to ~n*20000 for stm; plus a dense band around the seam)."""
    rng = np.random.default_rng(37)
    xs = np.concatenate([
        rng.uniform(-1.0, 1.0, 2000),
        rng.uniform(-3.2, 3.2, 2000),          # anomaly range
        rng.uniform(0.0, 0.0011313 * 20000, 2000),  # n*t for stm
        rng.uniform(-10.0, 10.0, 2000),
    ])
    lines = "\n".join(float(x).hex() for x in xs) + "\n"
    r = subprocess.run([str(cr_fn_exe)], input=lines, capture_output=True,
                       text=True, check=True)
    out = [ln.split() for ln in r.stdout.strip().split("\n")]
    assert len(out) == len(xs)
    for x, (cs, cc) in zip(xs, out):
        assert float.fromhex(cs) == croracle.cr_sin(float(x)), x
        assert float.fromhex(cc) == croracle.cr_cos(float(x)), x


def _stm_oracle(n, t):
    """cw.stm with the correctly-rounded sin/cos oracle substituted for
    libm's — all other arithmetic is plain float64, identical to the
    emitted kernel, so any difference from the C output would be a real
    translation bug, not rounding."""
    s = croracle.cr_sin(n * t)
    c = croracle.cr_cos(n * t)
    real_sin, real_cos = math.sin, math.cos
    math.sin = lambda a: s if a == n * t else real_sin(a)
    math.cos = lambda a: c if a == n * t else real_cos(a)
    try:
        return cw.stm(n, t)
    finally:
        math.sin, math.cos = real_sin, real_cos


@pytest.fixture(scope="module")
def stm_cr_exe(tmp_path_factory):
    if GCC is None:
        pytest.skip("gcc not available")
    d = tmp_path_factory.mktemp("stmcr")
    src = cemit.emit_module([cw.stm], "stm_cr", correctly_rounded=True)
    assert "cr_sin" in src and '#include "coremath.h"' in src
    (d / "k.c").write_text(src)
    drv = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
void podium_cw_stm(double n, double t, double out[6][6]);
int main(void){char l[128];while(fgets(l,sizeof l,stdin)){
  double n=strtod(strtok(l," \n"),0),t=strtod(strtok(0," \n"),0);
  double o[6][6];podium_cw_stm(n,t,o);
  for(int i=0;i<36;i++)printf(i?" %a":"%a",((double*)o)[i]);
  printf("\n");}return 0;}
"""
    (d / "d.c").write_text(drv)
    exe = d / "stmcr"
    subprocess.run(
        [GCC, "-std=c99", "-O2", "-ffp-contract=off", "-I", str(CM),
         str(d / "k.c"), str(d / "d.c"), str(CM / "sin.c"),
         str(CM / "cos.c"), "-lm", "-o", str(exe)],
        check=True, capture_output=True, text=True)
    return exe


@pytest.mark.slow
def test_stm_correctly_rounded_is_bit_exact(stm_cr_exe):
    """stm in CR mode == the CR-oracle stm reference, BIT-EXACT — no
    tolerance, no incidence bound. This retires the 21/72000 exception
    that libm sin/cos forced in test_cemit."""
    rng = np.random.default_rng(1)
    vecs = np.column_stack([rng.uniform(1e-4, 2e-3, 3000),
                            rng.uniform(0.0, 20_000.0, 3000)])
    lines = "\n".join(f"{float(n).hex()} {float(t).hex()}"
                      for n, t in vecs) + "\n"
    r = subprocess.run([str(stm_cr_exe)], input=lines, capture_output=True,
                       text=True, check=True)
    rows = [ln.split() for ln in r.stdout.strip().split("\n")]
    assert len(rows) == len(vecs)
    mism = 0
    for (n, t), row in zip(vecs, rows):
        got = np.array([float.fromhex(x) for x in row])
        ref = _stm_oracle(float(n), float(t)).ravel()
        if not np.array_equal(got, ref):
            mism += 1
    assert mism == 0, f"{mism}/{len(vecs)} not bit-exact"
