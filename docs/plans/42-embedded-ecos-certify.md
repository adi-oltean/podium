# 42 — Embedded ECOS solve of a Layer-0 SOCP with exact KKT certification

GitHub issue: https://github.com/adi-oltean/podium/issues/42

## Fix (landed) — `kkt.certify_ecos`

Closes the embedded-solver + verified-KKT loop on a REAL Layer-0
guidance problem. `certify_ecos(problem)`:

1. `problem.get_problem_data(cp.ECOS)` — cvxpy compiles the guidance
   SOCP to ECOS standard form (c, A, b, G, h, cone dims).
2. Convert G/A to csc, read `dims.nonneg` / `dims.soc`.
3. `ecos.solve(...)` — the EMBEDDED ECOS solver (the branchless,
   library-free interior-point C built for flight targets, invoked
   through its Python binding) returns native standard-form x/y/z/s.
4. `verify_socp(...)` (#41) re-verifies that solution in exact
   rational arithmetic.

Returns `(ecos_solution, SOCPReport)`. cvxpy/ecos/scipy are optional,
imported lazily so `podium.verify` stays dependency-free.

## Receipt (tests/test_kkt.py)

A min-fuel rendezvous SOCP — per-step thrust cones ||u_k|| <= t_k,
minimize sum t_k, subject to the CW-STM reach condition over 4 impulse
nodes — is solved by the embedded ECOS solver and certified exactly:
stationarity/conic residual/comp-slack all tiny, both cones satisfied,
certified within 1e-4; the exact-certified objective matches the
solver's pcost; and a burn is active (thrust cone tight). A tampered
primal is rejected at the checker level (the certificate is not a
rubber stamp).

## Why this is the online analogue of the offline certificates

Podium already re-checks offline artifacts exactly — barrier
certificates (#20), golden vectors (#26/#34/#37/#38), the mission
audit bundle. `certify_ecos` extends that discipline to the ONLINE
solver: the flight guidance solver's answer is trusted only after an
independent exact re-verification of its KKT optimality, not on the
solver's own float convergence flag.

## Deferred

Feeding a cvxpygen/QOCOGEN GENERATED-C solver's KKT dump to the
checker (vs the ECOS binding) — blocked only by cvxpygen's Julia
sidecar, and semantically identical since ECOS is the embedded C;
exponential/PSD cones; certifying inside the guidance planners
themselves (return a certified plan).

## Push/merge instructions

Single commit on main: `42 — Embedded ECOS solve + exact KKT
certification (#42)`; push; close.
