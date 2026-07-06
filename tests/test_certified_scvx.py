"""Certified successive convexification (P3 one-cut experiment).

A 2-D min-energy transfer must round a spherical keep-out. SCvx
linearizes the nonconvex ||r|| >= R into a half-space cut at each node;
applied HARD with a rational-unit normal, each cut is a SOUND convex
inner-approximation, machine-certified exact (scvx_cut). The payoff we
assert: EVERY node of EVERY iterate -- not just the converged one --
satisfies the true nonconvex keep-out, and each cut's soundness is an
exact certificate that can be re-synthesized from an untrusted float
SDP and validated.
"""

import math
from fractions import Fraction as F

import numpy as np
import pytest

from podium.verify import scvx_cut, sos

R = 5.0                                   # keep-out radius (origin)
N = 24                                    # nodes
DT = 0.5
START = np.array([-8.0, -1.0])
GOAL = np.array([8.0, -1.0])


def test_snap_returns_exact_rational_unit():
    """The snapped normal is a rational unit vector exactly (so the
    certificate arithmetic stays in Q)."""
    for vx, vy in [(1.0, 0.3), (-0.7, 0.7), (0.1, -0.99), (-1.0, -0.2)]:
        nx, ny = scvx_cut.snap_rational_unit(vx, vy)
        assert nx * nx + ny * ny == 1
        assert float(nx) * vx + float(ny) * vy > 0     # same half-plane


def test_halfspace_cut_is_sound_certificate():
    """A rational-unit keep-out cut certifies as a sound inner-approx:
    the SOS block and the S-procedure identity both check exactly."""
    n = scvx_cut.snap_rational_unit(0.6, 0.8)
    rep = scvx_cut.certify_halfspace_koz(n, F(5), (F(2), F(-3)))
    assert rep.sound(), rep.problems
    assert rep.sos_certified and rep.identity_ok


def test_non_unit_normal_is_rejected():
    """A non-unit normal is not accepted -- the certificate requires an
    exact rational unit vector."""
    rep = scvx_cut.certify_halfspace_koz((F(1), F(1)), F(5), (F(0), F(0)))
    assert not rep.sound()
    assert any("unit" in p for p in rep.problems)


def _initial_reference():
    """A feasible detour: a semicircle of radius R+1.5 above the KOZ,
    with the straight run-in/run-out to start and goal."""
    ref = np.zeros((N, 2))
    for k in range(N):
        s = k / (N - 1)
        ang = math.pi * (1.0 - s)          # pi -> 0, sweeping over the top
        arc = (R + 1.5) * np.array([math.cos(ang), math.sin(ang)])
        ref[k] = (1 - s) * START + s * GOAL + arc * math.sin(math.pi * s)
    ref[0], ref[-1] = START, GOAL
    return ref


def _subproblem(ref):
    """Min sum ||a||^2 double-integrator transfer with HARD rational-unit
    keep-out cuts at every node. Returns (nodes, normals) or None."""
    cp = pytest.importorskip("cvxpy")
    r = cp.Variable((N, 2))
    v = cp.Variable((N, 2))
    a = cp.Variable((N - 1, 2))
    cons = [r[0] == START, v[0] == np.zeros(2),
            r[N - 1] == GOAL, v[N - 1] == np.zeros(2)]
    for k in range(N - 1):
        cons += [r[k + 1] == r[k] + DT * v[k] + 0.5 * DT**2 * a[k],
                 v[k + 1] == v[k] + DT * a[k]]
    normals = []
    for k in range(N):
        d = ref[k] / max(np.linalg.norm(ref[k]), 1e-9)
        nx, ny = scvx_cut.snap_rational_unit(float(d[0]), float(d[1]))
        normals.append((nx, ny))
        cons.append(float(nx) * r[k, 0] + float(ny) * r[k, 1] >= R)
    prob = cp.Problem(cp.Minimize(cp.sum_squares(a)), cons)
    prob.solve(solver=cp.CLARABEL)
    if r.value is None:
        return None
    return r.value, normals


@pytest.mark.slow
def test_every_iterate_is_certified_feasible():
    """The headline: with hard certified cuts, every node of every SCvx
    iterate satisfies the true nonconvex keep-out ||r|| >= R, and every
    cut is an exact sound inner-approximation."""
    ref = _initial_reference()
    # the initial reference is itself node-feasible
    assert all(np.linalg.norm(ref[k]) >= R - 1e-9 for k in range(N))

    costs = []
    for _ in range(6):
        out = _subproblem(ref)
        assert out is not None, "hard-cut subproblem stayed feasible"
        nodes, normals = out

        # (1) per-iterate feasibility for the TRUE nonconvex constraint
        for k in range(N):
            assert np.linalg.norm(nodes[k]) >= R - 1e-6, (k, nodes[k])

        # (2) every cut is a machine-checked sound inner-approximation
        for nx, ny in normals:
            rep = scvx_cut.certify_halfspace_koz((nx, ny), F(5),
                                                 (F(0), F(0)))
            assert rep.sound(), rep.problems

        costs.append(float(np.sum(np.diff(nodes, axis=0) ** 2)))
        ref = nodes
    # SCvx improves (or holds) the path length as it tightens
    assert costs[-1] <= costs[0] + 1e-6


@pytest.mark.slow
def test_cut_certificate_validated_from_untrusted_sdp():
    """One cut's SOS soundness Gram, synthesized by an untrusted float
    SDP, is validated to an exact rational certificate (the P3
    float-synthesis -> exact-certificate mechanism, on the keep-out
    cut)."""
    cp = pytest.importorskip("cvxpy")
    n = scvx_cut.snap_rational_unit(0.6, 0.8)
    basis, lhs = scvx_cut.gram_constraints(n, F(5), (F(0), F(0)))

    # float SDP: find a PSD Gram matching the lhs coefficients
    g = cp.Variable((3, 3), symmetric=True)
    cons = [g >> 0]
    prods: dict[tuple, list] = {}
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            m = (bi[0] + bj[0], bi[1] + bj[1])
            prods.setdefault(m, []).append(g[i, j])
    for m, terms in prods.items():
        coeff = float(lhs.get(m, F(0)))
        cons.append(sum(terms) == coeff)
    cp.Problem(cp.Minimize(0), cons).solve(solver=cp.CLARABEL)

    g_exact = sos.validate_gram(lhs, basis, g.value.tolist())
    assert g_exact is not None
    ok, prob = sos.is_sos(lhs, basis, g_exact)
    assert ok, prob
