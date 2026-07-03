# 26 — CVXPYgen embedded solver generation

GitHub issue: https://github.com/adi-oltean/podium/issues/26

## Fix (landed) — `podium.emit.solvergen`

Fixed-grid Layer-0 rendezvous (grid + CW dynamics baked, x0/xf live
parameters — the embedded-deployment shape) generated to a
self-contained C tree via CVXPYgen with the ECOS backend;
`build_and_run` compiles it with plain gcc — no cmake, no Python — and
runs the generated example.

Build-recipe knowledge paid for and encoded in the module: -std=gnu99
(ECOS timers need POSIX timespec), -fcommon (pre-C11 tentative
definitions in the generated headers duplicate under gcc>=10),
runecos*/demo sources excluded (own mains), SuiteSparse_config include
path explicit, and the LDL/AMD long-int defines must NOT be set (cpg
uses int indices).

## Receipt (green)

The embedded ECOS binary reproduces the cvxpy/Clarabel optimum for the
same instance to 1e-5 relative (measured: 1.201844 vs 1.2018436700) —
cross-SOLVER agreement, standalone, zero Python at runtime.

## Dependency decision

cvxpygen is NOT a dev dependency: importing it pulls a Julia sidecar
through a pdaqp transitive dependency (measured: full Julia install at
first import). Generation is a local/offline step; the test
importorskips, so CI skips it. Revisit if cvxpygen drops the sidecar.

## Deferred (v0.6, already authored)

QOCOGEN alternate backend; the verified-KKT-checker pattern
(certificate checked by exact/interval arithmetic) wrapping whichever
backend ships.

## Push/merge instructions

Single commit on main: `26 — CVXPYgen embedded solvers (#26)`; push;
close.
