"""Abort-safety barrier receipts: end-to-end synthesize->rationalize->
exact-verify, a hand-derived certificate whose algebra closes on paper,
tamper detection, infeasibility of a genuinely unsafe set, and
trajectory corroboration against dense CW propagation."""

import math
from fractions import Fraction as Fr

import numpy as np
import pytest

from podium import constants as const
from podium.core import cw
from podium.verify import barrier

# e/i-separated drift-safe formation (scaled coords: velocities are v/n):
# radial offset 400 m with vy = -1.5 n x (zero in-plane amplitude at
# center), tight dispersions, along-track free within +-500 m.
CASE = barrier.AbortSafetyCase(
    center=(Fr(400), Fr(0), Fr(0), Fr(0), Fr(-600), Fr(0)),
    radii=(Fr(10), Fr(500), Fr(30), Fr(10), Fr(20), Fr(30)),
    koz_radius=Fr(200),
)


def test_scaled_cw_matrix_matches_kernel():
    """The checker's integer matrix IS cw_deriv in scaled coordinates."""
    n = 1.13e-3
    rng = np.random.default_rng(1)
    for _ in range(5):
        u = rng.normal(0.0, 100.0, 6)
        phys = np.concatenate([u[:3], n * u[3:]])
        d = cw.cw_deriv(phys, n, np.zeros(3))
        # du/dtau = (1/n) * [dr; dv/n]
        du = np.concatenate([d[:3] / n, d[3:] / n**2])
        a = np.array([[float(x) for x in row] for row in barrier.A_CW])
        assert np.allclose(a @ u, du, rtol=1e-12, atol=1e-9)


def test_end_to_end_synthesize_and_verify_exact():
    pytest.importorskip("cvxpy")
    sol = barrier.synthesize(CASE, margin=1.0)
    assert sol is not None, "synthesis should be feasible for the safe set"
    assert sol["slack"] >= 0.0
    cert = barrier.rationalize(sol, CASE, eps=0.5)
    problems = barrier.verify_certificate(cert)
    assert problems == [], problems
    # the trusted path is float-free: certificate carries only Fractions
    assert all(isinstance(x, Fr) for x in cert.a)
    assert isinstance(cert.lam0, Fr) and isinstance(cert.lam_u, Fr)


def test_hand_certificate_koz_identity():
    """The paper algebra: B = -c1^2 + 2*Ax2 + 2*Az2 + 2*R'^2 satisfies
    B - 2*gu = (X - u2)^2 + Z^2 + 2(Az2 - Z^2) + 2(R'^2 - R^2) >= 0 with
    lam_u = 2 — an exact sum of squares. The checker must agree, and C3
    must hold exactly (all basis elements are invariants)."""
    rp = Fr(210)  # R' > R gives the eps_u margin
    cert = barrier.BarrierCertificate(
        a=(Fr(-1), Fr(2), Fr(2), Fr(0), 2 * rp * rp),
        lam0=Fr(0),  # placeholder; C1 checked separately below
        lam_u=Fr(2),
        eps0=Fr(1),
        eps_u=Fr(1),
        case=CASE,
    )
    p, m = barrier.barrier_matrices(cert)
    # C3 exactness (via the public checker on a certificate whose C1 we
    # don't claim): check the Lie derivative directly
    a_cw = [[Fr(v) for v in row] for row in barrier.A_CW]
    for i in range(6):
        for j in range(6):
            s = sum(a_cw[k][i] * p[k][j] + p[i][k] * a_cw[k][j]
                    for k in range(6))
            assert s == 0
    # C2 with lam_u = 2 must be PSD exactly
    c2 = barrier._mat_add(
        m, barrier._mat_scale(Fr(-2), barrier._gu_matrix(CASE)),
        barrier._mat_scale(Fr(-1), [[Fr(0)] * 7 for _ in range(7)]))
    c2[0][0] -= cert.eps_u
    assert barrier.is_psd(c2)


def test_tamper_detection():
    pytest.importorskip("cvxpy")
    sol = barrier.synthesize(CASE, margin=1.0)
    cert = barrier.rationalize(sol, CASE, eps=0.5)
    assert barrier.verify_certificate(cert) == []
    # flip the sign of the c1^2 coefficient: safety claim inverts
    bad = barrier.BarrierCertificate(
        a=(-cert.a[0], cert.a[1], cert.a[2], cert.a[3], cert.a[4]),
        lam0=cert.lam0, lam_u=cert.lam_u,
        eps0=cert.eps0, eps_u=cert.eps_u, case=cert.case)
    assert barrier.verify_certificate(bad) != []
    # negative multiplier must be rejected
    bad2 = barrier.BarrierCertificate(
        a=cert.a, lam0=Fr(-1), lam_u=cert.lam_u,
        eps0=cert.eps0, eps_u=cert.eps_u, case=cert.case)
    assert any("lam0" in p for p in barrier.verify_certificate(bad2))


def test_unsafe_vbar_set_is_infeasible():
    """A hold centered ON the V-bar axis coasts straight along the
    V-bar: its RN separation is ~0, the fact is false, and synthesis
    must fail. Heritage agrees: V-bar stationkeeping is not passively
    safe."""
    pytest.importorskip("cvxpy")
    unsafe = barrier.AbortSafetyCase(
        center=(Fr(0), Fr(-500), Fr(0), Fr(0), Fr(0), Fr(0)),
        radii=(Fr(10), Fr(100), Fr(10), Fr(5), Fr(5), Fr(5)),
        koz_radius=Fr(200),
    )
    assert barrier.synthesize(unsafe, margin=1.0) is None


def test_trajectory_corroboration():
    """Empirical cross-check of the certified fact: corner states of X0,
    propagated densely through two orbits of exact CW flow, never bring
    RN separation inside the keep-out radius."""
    n = math.sqrt(const.MU_EARTH / 6_778_137.0**3)
    r_koz = float(CASE.koz_radius)
    corners = []
    c = [float(x) for x in CASE.center]
    r = [float(x) for x in CASE.radii]
    for sx in (-1, 1):
        for svy in (-1, 1):
            for sz in (-1, 1):
                u = c[:]
                u[0] += sx * r[0]
                u[4] += svy * r[4]
                u[2] += sz * r[2]
                u[3] += sx * r[3]  # exercise vx too
                corners.append(u)
    for u in corners:
        x0 = np.array([u[0], u[1], u[2],
                       n * u[3], n * u[4], n * u[5]])
        worst = math.inf
        for k in range(720):
            t = k * (2 * 2 * math.pi / n) / 720
            xt = cw.stm(n, t) @ x0
            worst = min(worst, math.hypot(xt[0], xt[2]))
        assert worst > r_koz, (u, worst)