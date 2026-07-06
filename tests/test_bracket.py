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
    the min-norm point on that half-space. The bound is gated on EXACT
    feasibility (certify_upper_bound), not on the KKT tolerance."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    sound = scvx_cut.certify_halfspace_koz(n, R, C).sound()
    ncx = n[0] * C[0] + n[1] * C[1]
    t = R + ncx                                     # x = t n on the cut
    x = [t * n[0], t * n[1]]
    kkt_ok = kkt.verify_qp(
        p=[[F(2), F(0)], [F(0), F(2)]], q=[F(0), F(0)],
        g=[[-n[0], -n[1]]], h=[-(R + ncx)], a=[], b=[],
        x=x, mu=[2 * t], nu=[]).certified()
    j_ub = bracket.certify_upper_bound(p0, q0, r0, p1, q1, r1, x)
    return j_ub, sound, kkt_ok


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


def test_dual_value_gives_certified_lower_bounds():
    """g(lam) = dual_value(...) is an exact certified lower bound at every
    rational lam in the PD region, and equals J* at the optimal lam*=2/5."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    for lam in [F(1, 3), F(3, 7), F(41, 100), F(2, 5)]:
        g = bracket.dual_value(p0, q0, r0, p1, q1, r1, lam)
        assert g is not None
        assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, lam, g)
        assert g <= J_STAR
    assert bracket.dual_value(p0, q0, r0, p1, q1, r1, F(2, 5)) == J_STAR


def test_recovery_converges_to_optimum_as_lambda_approaches_lam_star():
    """As lam -> lam*=2/5 the certified rational lower bound g(lam) -> J*,
    with the gap shrinking quadratically (~100x per 10x in |lam-lam*|)."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    gaps = []
    for lam in [F(1, 3), F(39, 100), F(399, 1000)]:   # |lam-2/5| = 1/15, 1/100, 1/1000
        g = bracket.dual_value(p0, q0, r0, p1, q1, r1, lam)
        gaps.append(J_STAR - g)
    assert gaps[0] > gaps[1] > gaps[2] > 0            # monotone to 0
    # quadratic: a 10x closer lam cuts the gap by ~100x (allow slack)
    assert gaps[1] / gaps[2] > 50


def test_recover_exact_optimum_from_float_dual():
    """Rounding a float SDP dual (~0.4003) to low denominator recovers the
    exact rational optimum lam*=2/5, giving t = J* exactly -- a closed
    exact bracket recovered from an untrusted floating-point solve."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    rec = bracket.recover_lower_bound(p0, q0, r0, p1, q1, r1, 0.4003,
                                      max_den=100)
    assert rec is not None
    lam, t = rec
    assert lam == F(2, 5) and t == J_STAR             # exact recovery
    # and it composes into a closed exact bracket with the upper leg
    j_ub, _sound, _kkt = _upper_bound((F(-1), F(0)))
    assert bracket.closes(t, j_ub) and t == J_STAR


def test_certified_optimum_binds_provenance_across_problems():
    """certified_optimum binds both legs to ONE problem's data, so a lower
    certificate valid for problem A cannot form a false closed bracket with
    an upper point from a different problem B (the provenance hole the
    round-3 ChatGPT audit flagged in the bare `closes` combiner)."""
    # Problem A: min x^2 s.t. 1 >= 0.  J*_A = 0; t=0 IS a valid lower bound.
    a = ([[F(1)]], [F(0)], F(0), [[F(0)]], [F(0)], F(1))
    t_lb, j_ub, closed = bracket.certified_optimum(*a, F(0), F(0), [F(0)])
    assert closed and t_lb == F(0) == j_ub          # correctly certifies J*=0

    # Problem B: min x^2 - 1 s.t. 1 >= 0.  J*_B = -1; t=0 is NOT <= J*_B.
    b = ([[F(1)]], [F(0)], F(-1), [[F(0)]], [F(0)], F(1))
    # the naive combiner would falsely close: closes(0, f0_B(1)=0) is True
    assert bracket.closes(F(0), F(0))
    # but certified_optimum re-checks the lower leg against B's data and
    # rejects t=0 (M not PSD), so no false closure forms.
    t_lb_b, j_ub_b, closed_b = bracket.certified_optimum(
        *b, F(0), F(0), [F(1)])
    assert t_lb_b is None                            # t=0 not certified for B
    assert not closed_b


def test_multi_constraint_bracket():
    """Theorem 3 (multiple constraints): the exact-rational bracket
    generalizes to m constraints. Two keep-outs
    min ||x||^2 s.t. ||x-(2,0)|| >= 3 and ||x+(2,0)|| >= 3 have J* = 5 at
    (0, +/-sqrt5) with BOTH constraints active. certify_lower_bound_multi
    certifies t = 5 (rejects 6); the upper bound gates on feasibility for
    ALL constraints; the exact bracket contains J*."""
    p0 = [[F(1), F(0)], [F(0), F(1)]]
    q0 = [F(0), F(0)]
    r0 = F(0)

    def keepout(cx, cy, rad):
        return ([[F(1), F(0)], [F(0), F(1)]], [-2 * cx, -2 * cy],
                cx * cx + cy * cy - rad * rad)

    cons = [keepout(F(2), F(0), F(3)), keepout(F(-2), F(0), F(3))]
    half = [F(1, 2), F(1, 2)]
    assert bracket.certify_lower_bound_multi(p0, q0, r0, cons, half, F(5))
    assert not bracket.certify_lower_bound_multi(p0, q0, r0, cons, half, F(6))
    # a loose multiplier gives a valid but weaker certified lower bound
    assert bracket.certify_lower_bound_multi(p0, q0, r0, cons,
                                             [F(1, 4), F(1, 4)], F(5, 2))
    # upper bound: exactly feasible for BOTH keep-outs -> valid; inside one
    # ball -> rejected (None)
    j_ub = bracket.certify_upper_bound_multi(p0, q0, r0, cons, [F(0), F(3)])
    assert j_ub == F(9)                              # 5 <= J* <= 9
    assert bracket.certify_upper_bound_multi(
        p0, q0, r0, cons, [F(0), F(0)]) is None      # origin is inside both


def test_hard_case_singular_optimum():
    """Theorem 2 (the singular 'hard case'): when A(lam*)=P0-lam*P1 is
    rank-deficient at the dual optimum (the trust-region hard case),
    dual_value returns None (it needs A > 0) -- yet the EXACT certificate
    at (lam*, J*) still verifies, so soundness is unaffected; only the
    smooth recovery is. The interior approximation from the PD side
    converges to J* only LINEARLY (not quadratically).

    Instance: min -x1^2 + 2 x2^2 + 2 x2 s.t. 1 - ||x||^2 >= 0; A(lam) =
    diag(-1+lam, 2+lam) is singular at lam*=1; analytic J* = -4/3.
    """
    p0 = [[F(-1), F(0)], [F(0), F(2)]]
    q0 = [F(0), F(2)]
    r0 = F(0)
    p1 = [[F(-1), F(0)], [F(0), F(-1)]]
    q1 = [F(0), F(0)]
    r1 = F(1)
    j = F(-4, 3)
    # A(1) = diag(0, 3) is singular -> dual_value cannot use it
    assert bracket.dual_value(p0, q0, r0, p1, q1, r1, F(1)) is None
    # but the exact rank-deficient certificate at (lam*, J*) verifies
    assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, F(1), j)
    # interior approximation (lam -> 1+) -> J*, LINEARLY: gap shrinks ~10x
    # per 10x closer to lam*=1 (quadratic recovery would be ~100x)
    ga = j - bracket.dual_value(p0, q0, r0, p1, q1, r1, F(101, 100))
    gb = j - bracket.dual_value(p0, q0, r0, p1, q1, r1, F(1001, 1000))
    assert ga > gb > 0
    assert 8 < ga / gb < 12


def test_upper_bound_rejects_tolerance_feasible_point():
    """A soundness hole the round-2 audit found: a point that PASSES
    kkt.verify_qp at the default tolerance can be EXACTLY infeasible for
    the true keep-out, with objective BELOW J* -- which would close the
    bracket beneath the global optimum. certify_upper_bound gates on exact
    feasibility, not the solver tolerance, so it rejects such a point and
    the bracket cannot form a false closure."""
    p0, q0, r0, p1, q1, r1 = _qcqp()
    eps = F(1, 10**10)
    x = [F(-2) + eps, F(0)]                # just inside the keep-out
    # it passes KKT at the default tolerance (violation 1e-10 <= 1e-9)...
    rep = kkt.verify_qp(
        p=[[F(2), F(0)], [F(0), F(2)]], q=[F(0), F(0)],
        g=[[F(1), F(0)]], h=[F(-2)], a=[], b=[],
        x=x, mu=[F(4) - 2 * eps], nu=[])
    assert rep.certified()
    # ...but it is EXACTLY infeasible, so it is not a valid upper bound
    j_ub = bracket.certify_upper_bound(p0, q0, r0, p1, q1, r1, x)
    assert j_ub is None
    assert not bracket.closes(F(4), j_ub)   # no false closed bracket
