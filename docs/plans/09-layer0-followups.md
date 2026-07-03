# 09 — Layer-0 follow-ups: LCvx, passive-safety constraints, QP tracking

GitHub issue: https://github.com/adi-oltean/podium/issues/9

## Problem

Layer-0 (#8) planned impulsive burns only, treated passive safety as an
a-posteriori check, and had no finite-burn (min-throttle) capability —
the constraint class real thrusters impose. LCvx relaxations are lossless
in continuous time under controllability/normality; in discrete time the
guarantee weakens (violations bounded by the state dimension,
arXiv:2410.09748), so the relaxation must ship with a validity audit, not
a hand wave.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `src/podium/guidance/convex.py` | PARTIAL | FiniteBurnPlanner + LCvx audit; PassiveSafetySpec; SafeSetSpec; qp_tracking |
| `tests/test_convex.py` | PARTIAL | new receipts |
| `docs/roadmap.md` | PARTIAL | check off LCvx/scenario items |

## Fix

1. **FiniteBurnPlanner (LCvx)**: ZOH CW discrete dynamics (exact
   `lqr.cw_discrete`), continuous accel u_k with thrust annulus
   rho_min <= ||u_k|| <= rho_max via the classical slack relaxation
   (||u_k|| <= G_k, rho bounds on G_k, min sum G_k dt). Shipped validity
   checks: (a) controllability precondition (Kalman rank of (Ad, Bd));
   (b) post-solve losslessness audit — per-node gaps G_k - ||u_k||,
   the discrete-time theory bound (#inactive nodes <= n_x), and the
   max gap reported on the plan. CW-only v0 (YA ZOH input matrix is a
   follow-up).
2. **Breger-How passive-safety scenarios** in RendezvousPlanner:
   `PassiveSafetySpec(radius, horizon, n_samples, failure_nodes)` — for
   each failure node j, the free-drift trajectory from the pre-burn
   state X_j is linear in the decision variables, so KOZ avoidance at
   sampled drift times is one linear constraint per (node, sample):
   c @ X_j >= R with c = normal . S Phi(tau) folded numerically into a
   single cvxpy Parameter (keeps DPP). Normals from the same bounded
   two-pass reference refinement as the KOZ.
3. **Convex e/i safe-set terminal** in RoePlanner:
   `SafeSetSpec(direction, e_min, i_min, cone_angle, da_tol)` replaces
   the terminal equality with alignment cones + minimum magnitudes
   (convex sufficient condition for e/i separation); the exact
   `safety.rn_margin` scan verifies the achieved geometry post-solve.
4. **QP tracking objective** (`objective="qp_tracking"`): quadratic
   state-tracking to a Parameter reference + control effort, for
   MPC-style re-solves on the compiled problem.

## Acceptance Criteria

- [x] All new receipts green (suite, ruff, mypy)
- [x] LCvx: annulus respected; losslessness audit clean on the normal
      near-minimum-time scenario (max gap ~3e-9, 0 non-tight nodes,
      rho_min riding at 12/16 nodes); rho_min genuinely bites; replay
      exact; controllability check. ALSO: the audit provably catches the
      degenerate excess-capacity case (15 non-tight nodes > n_x bound,
      gap ~1.6e-3) — the discrete analogue of coast arcs, where the
      relaxation is invalid and must not be flown
- [x] Passive safety: free plan's worst failure drift 149 m inside the
      200 m sphere; guarded plan keeps every scenario at >= 200 m
      (hyperplane implies distance); cost +5%
- [x] Safe set: constraints met, rn_margin positive, alignment cone held
      (required internal ROE scaling x1e5 — raw magnitudes sat at
      interior-point tolerances and returned optimal_inaccurate)
- [x] QP tracking: DPP re-solve against different references on the
      compiled problem; tracking beats ignoring the reference
- [x] Roadmap updated; residuals filed (#12)

## Push/merge instructions

Single commit on main: `09 — Layer-0 follow-ups: LCvx, passive safety,
QP tracking (#9)`; push; close #9; file #12 for residuals.

## Verification steps

Full suite; inspect the LCvx audit numbers on the test problem.
