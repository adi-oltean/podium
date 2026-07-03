"""Embedded-solver receipts: the CVXPYgen/ECOS binary, compiled with
plain gcc and run with zero Python, reproduces the cvxpy/Clarabel
optimum. Skipped when cvxpygen/gcc are absent (cvxpygen is deliberately
not a dev dependency — see podium.emit.solvergen)."""

import shutil

import pytest

pytest.importorskip("cvxpygen")

from podium.emit import solvergen  # noqa: E402


@pytest.mark.slow
def test_generated_solver_matches_clarabel(tmp_path):
    if shutil.which("gcc") is None:
        pytest.skip("gcc not available")
    out = tmp_path / "gen"
    ref = solvergen.generate(str(out))
    assert ref > 0.0
    obj_c = solvergen.build_and_run(str(out))
    # different solver (ECOS vs Clarabel), same optimum
    assert abs(obj_c - ref) <= 1e-5 * max(1.0, abs(ref)), (obj_c, ref)