"""Embedded solver generation for Layer-0 problems (CVXPYgen/ECOS).

Import explicitly — cvxpygen is NOT a dev dependency (its import pulls
a Julia sidecar through a pdaqp transitive dependency; generation is a
local/offline step, not a CI step). The product is a self-contained C
tree: `build_and_run` compiles it with plain gcc — no cmake, no Python —
and runs the generated example, returning the reported objective. The
receipt is cross-solver: the embedded ECOS binary must reproduce the
cvxpy/Clarabel optimum for the same instance.

Build recipe knowledge (paid for once, encoded here):
  - -std=gnu99 (ECOS timers need POSIX timespec)
  - -fcommon   (pre-C11 tentative definitions in generated headers)
  - exclude runecos*.c (bring their own main) and test/demo sources
  - SuiteSparse_config include path must be explicit
"""

from __future__ import annotations

import math
import pathlib
import re
import subprocess

import numpy as np

from podium import constants as const
from podium.core import cw


def make_problem(k: int, dt: float, n: float) -> "tuple":
    """Fixed-grid impulsive rendezvous with x0/xf as parameters (the
    embedded-deployment shape: grid and dynamics baked, endpoints live)."""
    import cvxpy as cp

    x = cp.Variable((6, k + 1), name="x")
    v = cp.Variable((3, k + 1), name="v")
    x0 = cp.Parameter(6, name="x0")
    xf = cp.Parameter(6, name="xf")
    b = np.zeros((6, 3))
    b[3:6] = np.eye(3)
    phi = cw.stm(n, dt)
    cons = [x[:, 0] == x0, x[:, k] + b @ v[:, k] == xf]
    for i in range(k):
        cons.append(x[:, i + 1] == phi @ x[:, i] + (phi @ b) @ v[:, i])
    obj = cp.Minimize(cp.sum([cp.norm(v[:, i], 2) for i in range(k + 1)]))
    return cp.Problem(obj, cons), x0, xf


def generate(out_dir: str, k: int = 3, dt: float = 300.0,
             a_ref: float = 6_778_137.0,
             x0_val: np.ndarray | None = None,
             xf_val: np.ndarray | None = None) -> float:
    """Generate the C tree; returns the Python-side (Clarabel) optimum
    for the baked-in example instance (the cross-solver reference)."""
    import cvxpy as cp
    from cvxpygen import cpg

    n = math.sqrt(const.MU_EARTH / a_ref**3)
    prob, x0, xf = make_problem(k, dt, n)
    x0.value = (np.array([0.0, -500.0, 0.0, 0.0, 0.0, 0.0])
                if x0_val is None else x0_val)
    xf.value = np.zeros(6) if xf_val is None else xf_val
    prob.solve(solver=cp.CLARABEL)
    ref = float(prob.value)
    cpg.generate_code(prob, code_dir=out_dir, solver="ECOS", wrapper=False)
    return ref


def build_and_run(out_dir: str) -> float:
    """Compile the generated tree with gcc and run the example binary;
    returns the objective it reports."""
    c = pathlib.Path(out_dir) / "c"
    srcs = []
    for pat in ("src", "solver_code/src", "solver_code/external/ldl/src",
                "solver_code/external/amd/src",
                "solver_code/external/SuiteSparse_config"):
        for f in sorted((c / pat).glob("*.c")):
            if re.search(r"runtest|demo|_test|runecos", f.name):
                continue
            srcs.append(str(f))
    incs = ["-I" + str(c / p) for p in
            ("include", "solver_code/include",
             "solver_code/external/ldl/include",
             "solver_code/external/amd/include",
             "solver_code/external/SuiteSparse_config")]
    exe = c / "cpg_example_bin"
    subprocess.run(
        ["gcc", "-O2", "-std=gnu99", "-fcommon", "-o", str(exe)]
        + srcs + incs + ["-lm"],
        check=True, capture_output=True, text=True)
    r = subprocess.run([str(exe)], capture_output=True, text=True,
                       check=True)
    m = re.search(r"obj = ([-0-9.eE+]+)", r.stdout)
    if not m:
        raise RuntimeError(f"no objective in solver output: {r.stdout[:200]}")
    return float(m.group(1))
