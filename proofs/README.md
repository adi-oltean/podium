# Lean formalization — exact-rational trajectory-optimality certificates

Machine-checked (Lean 4 + Mathlib, no `sorry`) formalization of the mathematics behind the
podium exact-ℚ optimality certifier for block-banded trajectory QCQPs. **~60 theorems across 18
files, all `axioms`-clean** (only `propext`/`Classical.choice`/`Quot.sound`; several use none).
The `IsPSD` predicate decided here provably coincides with Mathlib's standard `Matrix.PosSemidef`
(`PosSemidefBridge`), so the checker decides the genuine library notion of PSD.

- **Check everything:** `bash check-all.sh` → `ALL PASS` (prints per-file result + `#print axioms`).
- **The mathematics, on paper:** `PROOFS.md` — the human-readable argument for every theorem,
  reviewable without reading Lean.
- **Per-theorem status + Mathlib support map:** `STATUS.md`.

## What is proved

**The trusted checker is correct — both directions.** The certifier accepts a symmetric matrix
iff it is PSD, decided by exact-ℚ `LDLᵀ`:
- *soundness* (`LDLTSound.ldlt_isPSD`): a nonnegative-pivot factorization ⇒ PSD;
- *completeness* (`SchurStep.ldlt_complete`): PSD ⇒ such a factorization exists;
- combined as one iff (`SchurStep.isPSD_iff_ldlt`) — the decision-procedure spec.
- the vector-state (`d>1`) block Riccati step (`BlockSchur.block_schur_complement`).

**The certificate implies global optimality — both legs.** From an accepted certificate:
- *lower bound* (S-procedure weak duality, scalar/CDT/general-`n`): `t ≤ J*`;
- *attainment* (`ClosureOptimum`): a KKT witness attains it, so the bracket closes to `f0(x̄)=J*`;
- assembled into one **headline theorem** (`Certificate.certified_optimum_pipeline`) and
  demonstrated on a **concrete keep-out QCQP over ℚ** (`Certificate.keepout_certified`).

**Supporting structure.** The PSD convex cone closed under aggregation + congruence
(`StructuralPSD`) — why the method is sound; general-`n` diagonal-dominance ⇒ PSD
(`GershgorinPSD`); horizon-independence via the DARE fixed point (`DareFixedPoint`, `MatrixDare`);
telescoping continuants and a constructed family with linear bit-growth (`TelescopingProduct`,
`ContinuantClosedForm`); the public technical note's theorems (`note/PodiumNote`).

## Honest scope

Every theorem here machine-checks a **classical or elementary** fact (weak duality, the
Schur-complement / `LDLᵀ` PSD criterion, Gershgorin, continuants). The value is a tolerance-free,
`sorry`-free **trusted path**, not new mathematics — state results as "machine-checked," not novel.
The genuine research delta is the fusion (exact-ℚ, machine-checked certifier specialized to RPOD
trajectory QCQPs), documented in the paper.

## The one remaining gap

Not mathematics: the **algorithmic refinement** — proving the concrete band-elimination code in
`podium.verify.riccati` computes exactly the factorization `ldlt_complete` shows exists. This is
implementation correctness (a natural formal-methods collaboration target), not a new theorem.
Also not mechanized: the analytic recovery *rates* (Taylor / one-sided derivatives) and a general
bit-complexity bound (open) — see `STATUS.md`.

## Build

The corpus builds against a pinned Mathlib, whose toolchain, revision and full
dependency graph are fixed in `ci/lean-toolchain`, `ci/lakefile.toml` and
`ci/lake-manifest.json`. Two scripts do everything:

```
bash setup.sh        # installs elan if absent, builds the pinned project into .leanwork/
bash check-all.sh    # type-checks every file, prints PASS/FAIL and the axiom audit
```

`setup.sh` fetches Mathlib's prebuilt cache rather than compiling it (minutes on a cold
machine, seconds afterwards). To reuse a Mathlib project you already have, point
`LEANWORK` at it; to use a `lake` that is not on `PATH`, set `LAKE`. A single file:

```
env -C .leanwork lake env lean <File>.lean
```

All 18 top-level files (19 with `note/PodiumNote.lean`) compile with exit 0 and no
`sorryAx`.
