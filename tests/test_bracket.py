"""Exact optimality-gap bracket for a nonconvex QCQP (opt-C prototype).

    min ||x||^2  s.t.  ||x - c||^2 >= R^2 ,  c=(3,0), R=5.

The origin sits inside the keep-out ball, so the problem is nonconvex;
the analytic optimum is the nearest boundary point x*=(-2,0), J*=4. We
bracket J* by two exact-rational certificates: an S-procedure LMI dual
lower bound and a sound-cut + KKT upper bound. When they meet, the
bracket closes into an exact certificate of the global optimum of a
nonconvex program.
"""

from fractions import Fraction as F

from podium.verify import bracket, kkt, scvx_cut

C = (F(3), F(0))
R = F(5)
J_STAR = F(4)                      # analytic global optimum


def _qcqp():
    return bracket.keepout_qcqp(C, R)


def _upper_bound(n):
    """Certified upper bound from a sound half-space cut n.(x-c) >= R and
    the KKT-certified min-norm point on that half-space."""
    sound = scvx_cut.certify_halfspace_koz(n, R, C).sound()
    ncx = n[0] * C[0] + n[1] * C[1]
    t = R + ncx                                     # x = t n on the cut
    x = [t * n[0], t * n[1]]
    rep = kkt.verify_qp(
        p=[[F(2), F(0)], [F(0), F(2)]], q=[F(0), F(0)],
        g=[[-n[0], -n[1]]], h=[-(R + ncx)], a=[], b=[],
        x=x, mu=[2 * t], nu=[])
    return x[0]**2 + x[1]**2, sound, rep.certified()


def test_lower_bound_certificate_is_exact_and_sound():
    """The S-procedure LMI certifies valid lower bounds (t <= J*) and
    rejects invalid ones (t > J*)."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    # dual optimum lambda* = 2/5 gives the tight bound t = 4 = J*
    assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, F(2, 5), F(4))
    # a smaller t is also certified (weak duality)...
    assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, F(2, 5),
                                       F(39, 10))
    # ...but t > J* is NOT a valid lower bound and must be rejected
    assert not bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, F(2, 5),
                                           F(41, 10))
    # a negative multiplier is not dual-feasible
    assert not bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, F(-1, 5),
                                           F(4))


def test_bracket_closes_to_exact_global_optimum():
    """Optimal cut n = unit(x*-c) = (-1,0): the certified upper bound
    equals the certified lower bound, so the bracket CLOSES -- an exact
    certificate of the nonconvex global optimum."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    j_lb = F(4)
    assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, F(2, 5), j_lb)

    j_ub, sound, kkt_ok = _upper_bound((F(-1), F(0)))
    assert sound and kkt_ok
    assert j_lb <= J_STAR <= j_ub          # valid bracket
    assert j_lb == j_ub == J_STAR          # closed -> exact global optimum


def test_bracket_is_a_real_exact_gap_for_a_suboptimal_cut():
    """A suboptimal (but still sound) rational-unit cut gives J_ub > J*,
    so the bracket is open -- a real optimality gap, both ends exact."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    j_lb = F(4)
    assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, F(2, 5), j_lb)

    n = scvx_cut.snap_rational_unit(-4.0, -3.0)
    j_ub, sound, kkt_ok = _upper_bound(n)
    assert sound and kkt_ok
    assert j_lb <= J_STAR < j_ub           # strict gap
    assert j_ub == F(169, 25)              # exact rational upper bound
