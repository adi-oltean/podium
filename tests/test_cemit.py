"""C-emitter receipts: tier-1 golden vectors (BIT-exact Python <-> C on
the host under pinned FP semantics), subset rejection, ACSL rendering,
and deterministic emission."""

import shutil
import subprocess

import numpy as np
import pytest

from podium.core import cw, quat
from podium.emit import cemit
from podium.verify import Interval, contract

GCC = shutil.which("gcc")

KERNELS = [quat.normalize, quat.multiply, quat.conjugate, quat.rotate,
           quat.deriv, quat.error, cw.mean_motion, cw.cw_deriv, cw.stm]

# per-function input specs: (list of (kind, size), sampler ranges)
CASES = {
    "normalize": ([("a", 4)], (-2.0, 2.0)),
    "multiply": ([("a", 4), ("a", 4)], (-1.0, 1.0)),
    "conjugate": ([("a", 4)], (-1.0, 1.0)),
    "rotate": ([("a", 4), ("a", 3)], (-1.0, 1.0)),
    "deriv": ([("a", 4), ("a", 3)], (-1.0, 1.0)),
    "error": ([("a", 4), ("a", 4)], (-1.0, 1.0)),
    "mean_motion": ([("s", 1), ("s", 1)], (1.0e6, 4.0e14)),
    "cw_deriv": ([("a", 6), ("s", 1), ("a", 3)], (-100.0, 100.0)),
    "stm": ([("s", 1), ("s", 1)], (1e-4, 2e-3)),
}
OUT_LEN = {"normalize": 4, "multiply": 4, "conjugate": 4, "rotate": 3,
           "deriv": 4, "error": 3, "mean_motion": 1, "cw_deriv": 6,
           "stm": 36}

_DRIVER = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
{prototypes}
int main(int argc, char **argv) {
    double in[16], out[64];
    char line[4096];
    while (fgets(line, sizeof line, stdin)) {
        int n = 0;
        char *tok = strtok(line, " \n");
        while (tok) { in[n++] = strtod(tok, NULL); tok = strtok(NULL, " \n"); }
        int m = 0;
        {dispatch}
        for (int i = 0; i < m; i++) printf(i ? " %a" : "%a", out[i]);
        printf("\n");
    }
    return 0;
}
"""


def _dispatch_block() -> str:
    rows = []
    for f in KERNELS:
        name = f.__name__
        spec, _ = CASES[name][0], CASES[name][1]
        c_name = f"podium_{f.__module__.rsplit('.', 1)[-1]}_{name}"
        args = []
        off = 0
        for kind, size in spec:
            if kind == "s":
                args.append(f"in[{off}]")
                off += 1
            else:
                args.append(f"&in[{off}]")
                off += size
        if OUT_LEN[name] == 1 and name == "mean_motion":
            call = f"out[0] = {c_name}({', '.join(args)}); m = 1;"
        else:
            call = (f"{c_name}({', '.join(args + ['(void *)out'])});"
                    f" m = {OUT_LEN[name]};")
        rows.append(f'if (!strcmp(argv[1], "{name}")) {{ {call} }}')
    return "\n        ".join(rows)


def _py_call(name, vec):
    f = dict((k.__name__, k) for k in KERNELS)[name]
    spec = CASES[name][0]
    args = []
    off = 0
    for kind, size in spec:
        if kind == "s":
            args.append(float(vec[off]))
            off += 1
        else:
            args.append(np.asarray(vec[off:off + size], dtype=np.float64))
            off += size
    r = f(*args)
    return np.atleast_1d(np.asarray(r, dtype=np.float64)).ravel()


def _vectors(name, n_vec, rng):
    spec, (lo, hi) = CASES[name]
    total = sum(s for _, s in spec)
    if name == "mean_motion":
        v = np.column_stack([rng.uniform(3.0e14, 4.5e14, n_vec),
                             rng.uniform(6.6e6, 8.0e6, n_vec)])
    elif name == "stm":
        v = np.column_stack([rng.uniform(1e-4, 2e-3, n_vec),
                             rng.uniform(0.0, 20_000.0, n_vec)])
    else:
        v = rng.uniform(lo, hi, (n_vec, total))
    if name == "normalize":
        v[-1, :] = 1e-200  # exercise the zero-guard branch
    if name == "error":
        v[: n_vec // 2, 0] = -np.abs(v[: n_vec // 2, 0])  # sign branch
    return v


@pytest.fixture(scope="module")
def compiled(tmp_path_factory):
    if GCC is None:
        pytest.skip("gcc not available")
    d = tmp_path_factory.mktemp("cemit")
    src = cemit.emit_module(KERNELS, "podium_kernels_test")
    (d / "kernels.c").write_text(src)
    protos = "\n".join(
        line for line in src.splitlines()
        if line.startswith(("void podium_", "double podium_")) and
        line.endswith(";"))
    driver = _DRIVER.replace("{prototypes}", protos) \
                    .replace("{dispatch}", _dispatch_block())
    (d / "driver.c").write_text(driver)
    exe = d / "kern"
    subprocess.run(
        [GCC, "-std=c99", "-O2", "-ffp-contract=off", "-Wall", "-Werror",
         "-o", str(exe), str(d / "kernels.c"), str(d / "driver.c"), "-lm"],
        check=True, capture_output=True, text=True)
    return exe


# Kernels whose only libm calls are IEEE-exact (sqrt is correctly
# rounded by the standard) are BIT-exact across toolchains. sin/cos are
# not correctly rounded in either libm, and this interpreter's (conda)
# libm differs from the system glibc that gcc links — measured: 21 of
# 72,000 stm values differ by exactly 1 ulp. Full cross-toolchain
# bit-exactness for transcendentals is what the roadmap's CORE-MATH
# item exists to buy; until then the honest tier-1 claim is:
# bit-exact for arithmetic+sqrt kernels, <=1 ulp at <0.1% incidence
# for sin/cos-bearing ones.
_TRANSCENDENTAL = {"stm"}


@pytest.mark.slow
@pytest.mark.parametrize("name", list(CASES))
def test_tier1_bit_exact(compiled, name):
    """2000 seeded vectors per kernel: C output equals Python output bit
    for bit (hex-float round trip) — except documented sin/cos libm
    divergence, bounded to 1 ulp at <0.1% incidence."""
    rng = np.random.default_rng(hash(name) % 2**32)
    vecs = _vectors(name, 2000, rng)
    # hex floats: exact in both directions (C99 strtod parses %a form)
    lines = "\n".join(" ".join(float(x).hex() for x in row)
                      for row in vecs) + "\n"
    r = subprocess.run([str(compiled), name], input=lines,
                       capture_output=True, text=True, check=True)
    out_lines = r.stdout.strip().split("\n")
    assert len(out_lines) == len(vecs)
    mismatches = 0
    total = 0
    for row, line in zip(vecs, out_lines):
        py = _py_call(name, row)
        c_vals = [float.fromhex(t) for t in line.split()]
        assert len(c_vals) == len(py)
        for a, b in zip(py, c_vals):
            total += 1
            if a != b and not (np.isnan(a) and np.isnan(b)):
                mismatches += 1
                # 1-ulp libm divergence amplified by the downstream
                # arithmetic: bounded by a few ulps, never more
                assert abs(a - b) <= 4 * np.spacing(max(abs(a), abs(b))), \
                    (a, b)
    if name in _TRANSCENDENTAL:
        assert mismatches <= 0.001 * total, f"{name}: {mismatches}/{total}"
    else:
        assert mismatches == 0, f"{name}: {mismatches} bit mismatches"


def test_rejects_outside_subset():
    def looped(x):  # noqa: ANN001
        s = 0.0
        for i in range(3):
            s = s + x[i]
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([looped])

    def heapy(x):  # noqa: ANN001
        return [x, x]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([heapy])


def test_acsl_rendering():
    @contract(n=Interval(1e-4, 1e-2), x=Interval(-10.0, 10.0))
    def scaled(x, n):  # noqa: ANN001
        out = np.empty(2)
        out[0] = x[0] * n
        out[1] = x[1] * n
        return out

    src = cemit.emit_module([scaled], "acsl_test")
    assert "requires 0.0001 <= n <= 0.01;" in src
    assert "\\forall integer i; 0 <= i < 2" in src
    assert "[in, range(-10.0,10.0)] x;" in src
    assert "/* [spec]" in src


def test_emission_deterministic():
    a = cemit.emit_module(KERNELS)
    b = cemit.emit_module(KERNELS)
    assert a == b