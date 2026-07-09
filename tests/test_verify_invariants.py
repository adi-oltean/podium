"""Invariants that were previously re-checked by unreachable runtime guards,
now asserted where the property is actually established (test it, don't guard
it with dead code)."""
from fractions import Fraction as F

from podium.verify import barrier, scvx_cut


def test_barrier_matrices_is_flow_invariant():
    # P is a combination of the conserved integer basis, so A'P + PA == 0
    # EXACTLY for ANY coefficients -- the C3 invariant verify_certificate used
    # to re-check per cert (it could never fail, since it depends on the basis,
    # not on cert.a).
    case = barrier.AbortSafetyCase(
        center=(F(400), F(0), F(0), F(0), F(-600), F(0)),
        radii=(F(10), F(500), F(30), F(10), F(20), F(30)), koz_radius=F(200))
    cert = barrier.BarrierCertificate(
        a=(F(3), F(-2), F(5), F(1), F(-4)),
        lam0=F(1), lam_u=F(1), eps0=F(1), eps_u=F(1), case=case)
    p, _m = barrier.barrier_matrices(cert)
    n = barrier._N
    a_cw = [[F(v) for v in row] for row in barrier.A_CW]
    lie = [[sum((a_cw[k][i] * p[k][j] + p[i][k] * a_cw[k][j] for k in range(n)),
                F(0)) for j in range(n)] for i in range(n)]
    assert all(lie[i][j] == 0 for i in range(n) for j in range(n))


def test_scvx_koz_cut_s_procedure_identity_holds():
    # lhs + 2R h == ||u||^2 - R^2 holds by construction of _koz_polys; this is
    # what the removed `if not identity_ok` branch used to assert.
    rep = scvx_cut.certify_halfspace_koz((F(3, 5), F(4, 5)), F(2), (F(1), F(1)))
    assert rep.identity_ok


def test_scvx_koz_cut_refuses_negative_radius():
    rep = scvx_cut.certify_halfspace_koz((F(3, 5), F(4, 5)), F(-1), (F(0), F(0)))
    assert any("radius" in p for p in rep.problems)
