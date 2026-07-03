# 22 — STL mission constraints in SCP via smooth robustness

GitHub issue: https://github.com/adi-oltean/podium/issues/22

## Problem

Mission timelines carry temporal constraints ("visit the inspection box
between 200 and 400 s") that neither Layer-0 nor plain PTR expressed.

## Fix (landed) — `EventuallyBoxSpec` in `podium.guidance.scp`

Timed-reach STL fragment inside the PTR loop. The encoding that is
actually SOUND (and the unsound one we caught):

- node box margins m_k = min_i(half_i - |r_i - c_i|) are CONCAVE ->
  they enter the subproblem EXACTLY as hypograph variables
  (mk <= each of 6 affine face expressions — convex constraints);
- the smooth max LSE/tau over margins is CONVEX in the margins -> it is
  tangent-linearized at the reference; a convex function dominates its
  tangent, so tangent >= eps + ln(K)/tau  =>  LSE >= target  =>  true
  max-margin >= eps. Conservatism ln(K)/tau is consumed explicitly
  (tau = 0.5 default: 3.2 m for a 5-node window).
- FAILURE MODE CAUGHT: first implementation linearized the concave
  margins with a single-face subgradient. Concave functions lie BELOW
  their linearizations, so the constraint was a relaxation — the
  optimizer exited the box through faces the subgradient didn't pick,
  producing a stable period-2 oscillation (rho alternating +4.2 /
  -86 with zero slack). Root-caused via the numeric probe, fixed by the
  hypograph split. Also: trust-region expansion is now floored
  (w_tr -> 0 destabilizes linearized constraints, observed).

STL shortfall (eps - true robustness) is folded into the PTR loop's
true-violation metric, so convergence requires the EXACT semantics.

## Receipts (all green)

- bite: without the spec the optimal trajectory has negative true
  robustness; with it, converged plans meet TRUE non-smooth robustness
  >= eps while the keep-out sphere stays clean and fuel rises;
- engine flight judged by the spec registry over a custom box-distance
  channel: eventually_below margin positive — guidance (smooth proxy),
  truth (nonlinear flight), and monitoring (STL registry) agree on the
  same temporal property.

## Deferred

MIP reference validation (no MIP dependency by policy), richer STL
fragments (until, nesting), robustness maximization, D-GMSR-style
smooth semantics for full formulas.

## Push/merge instructions

Single commit on main: `22 — STL in SCP (#22)`; push; close.
