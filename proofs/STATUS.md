# Lean formalization — feasibility assessment + status

Goal: mechanize the theorems of the exact-arithmetic certificate note
(`docs/exact-arithmetic-certificates/note.tex`) in Lean 4 + Mathlib. Toolchain and Mathlib
revision are pinned in `ci/` and built by `setup.sh`.

## Honest scope of these proofs
Every lemma here machine-checks a **classical or elementary** fact: weak duality (S-lemma
soundness), the continuant telescoping product (Bareiss/continuant theory), the Möbius/DARE
fixed point (periodic continued fractions), Gershgorin (diagonal-dominance PSD), and the
matrix lifting. **The value is the machine-checking (a tolerance-free, sorry-free trusted
path), not new mathematics.** State them in the paper as "machine-checked", not as novel
results. The genuine research delta is the fusion (see docs/pending-paper-revisions.md §E),
not these lemmas individually.

## Verdict
**Feasible.** Mathlib has the needed foundations — `Matrix.PosSemidef`, `ConvexOn` /
`IsMinOn`, quadratic forms, Schur-complement lemmas, `Nat` bit-size. A *full*
formalization is a multi-week effort (the PSD/matrix/convexity plumbing is the bulk), but
the load-bearing **soundness** results are tractable and are where mechanization is most
valuable (they certify the trusted checker). First proof already lands:
**`DareFixedPoint.lean` — the scalar DARE horizon-independence lemma, no `sorry`**
(`#print axioms` = the three standard axioms only).

## Per-theorem map (target → Mathlib support → difficulty)
| Theorem | Mathlib support | Difficulty | Priority |
|---|---|---|---|
| **DARE fixed point** (paper Lemma dare, scalar) | field arith | **DONE** | — |
| **S-procedure soundness** `λ≥0 ∧ minorant≥t ⇒ t≤J*` — single + multi (chain/CDT) + bracket (note thm:sound/thm:multi; paper Thm 1 lower leg / Prop) | ordered field, `Finset.sum_nonneg`, `linarith` | **DONE** (minorant level) | — |
| **Closure** (paper Thm 1): Lagrangian-min witness + compl. slackness ⇒ global min | `Finset.sum_nonneg`, `linarith` | **DONE** (`ClosureOptimum`) | — |
| **Telescoping** `∏ d_j = K_k` (the round-2 fraction-free correction) | field arith, induction | **DONE** | — |
| **Bit law** `bit(d_k)=4k+1` (paper Lemma bitlaw) | `Nat`, continuant recurrence, `Nat.size` | medium (bit-size fiddly) | 3 |
| **Schur lemma** (note lem:schur) | `Matrix.schur_complement*` | medium | 3 |
| **Generic bracket** rate (Prop) | dual value, first-order optimality | medium | 4 |
| **LDLᵀ soundness leg** `M=LᵀDL ∧ pivots≥0 ⇒ PSD` (paper Lemma band, accept ⇒ PSD) | `mulVec_mulVec`, `dotProduct_mulVec`, `Finset.sum_nonneg` | **DONE** (`LDLTSound`) | — |
| **LDLᵀ completeness — elimination STEP** (`IsPSD` ⇒ pivot ≥ 0 ∧ Schur complement PSD) | Fin.cons block algebra, minimize over pivot | **DONE** (`SchurStep`) | — |
| **LDLᵀ completeness — degenerate `a=0` key fact** (`IsPSD` ∧ pivot `0` ⇒ border `b=0`) | linear-form unboundedness | **DONE** (`SchurStep`) | — |
| **LDLᵀ completeness — full theorem** (`IsPSD` ⇒ `M = ∑ⱼ dⱼvⱼvⱼᵀ`, `dⱼ ≥ 0`) | induction peeling one pivot/rank-1 per step | **DONE** (`ldlt_complete`) | — |
| **Diagonal dominance ⇒ PSD, general `n`** (paper Lemma dd, Gershgorin) | quadratic form, complete-the-square, symmetry reindex | **DONE** (`GershgorinPSD`) | — |
| **Schur-complement preservation of dominance** (paper Lemma dd, second half) | Carlson–Markham (absent) | hard | 5 |
| **Recovery** nonsingular/singular (note thm:recovery/hard) | Taylor, Schur, hard case | hard | 5 |

## Highest-value target: the trusted-checker soundness
The single most valuable thing to mechanize is **S-procedure soundness** — that
`λ≥0` and `M(λ,t)⪰0` imply `t≤J*` — because it is exactly the guarantee the exact-ℚ
checker relies on. It is weak duality: for feasible `x` (so `f_k(x)≥0`) and `λ≥0`,
`f_0(x) ≥ f_0(x) − Σλ_k f_k(x) = [x;1]ᵀ M(λ,t) [x;1] + t ≥ t`. In Mathlib this is a
`PosSemidef.inner_nonneg`-style step over the lifted quadratic form. This mechanizes the
soundness of *both* the paper's bracket and the note's thm:sound in one lemma — do it next.

## Done
- `DareFixedPoint.lean` — `dare_fixed_point`, `dare_orbit_constant` (seeded Riccati orbit
  is constant = m on the Pell locus). Compiles clean, no `sorry`.
- `SProcedureSoundness.lean` — `sprocedure_lower_bound` (single keep-out / S-lemma),
  `sprocedure_lower_bound_multi` (chain / CDT, any `m`), `sprocedure_bracket_sound`. Over a
  general `LinearOrderedField` (so `ℚ` and `ℝ`). Compiles clean, no `sorry`. This is the
  weak-duality core of the trusted-checker guarantee; the remaining plumbing is the matrix
  identity `M(λ,t)⪰0 ⇔ minorant ≥ 0` (the PSD definition of the lifted quadratic form),
  which is where Mathlib's `Matrix.PosSemidef` would connect the checker to these lemmas.
- `RecoveryRate.lean` — `recovery_rate`: the `O(1/D²)` recovery bound (note `thm:recovery`),
  mechanized as an EXPLICIT bound (not asymptotic). For a concave-quadratic dual
  `g(λ) = g* − k(λ−λ*)²` (`k ≥ 0`, the nonsingular interior-maximizer Taylor model), rationalizing
  the multiplier to precision `1/D` gives `J* − g(λ) ≤ k/D²` — the certified gap shrinks
  QUADRATICALLY in the budget `D`. Also `recovery_rate_hard`: the singular / trust-region hard
  case (note `thm:hard`), a dual with a kink `g(λ) = g* − k|λ−λ*|` ⇒ `J* − g(λ) ≤ k/D` (LINEAR
  shrinkage). Both over ℚ, no `sorry`. (Only the fully general/asymptotic statement remains
  analytic; both concrete rates — nonsingular `O(1/D²)` and hard-case `O(1/D)` — are mechanized.)

- `ContinuantClosedForm.lean` — `continuant_closed_form`: the constructed continuant
  `K_{k+2}=5K_{k+1}−4K_k`, `K_0=1`, `K_1=4` has closed form `K_k = 4^k` (roots 4, 1), so
  `bit(K_k)=2k+1` grows linearly in the horizon — a concrete, machine-checked instance of the
  paper's family-specific `Θ(N)` per-entry bit-size (the general bound is open). No `sorry`.

- `TelescopingProduct.lean` — `pivotProd_telescopes`: the product of the tridiagonal LDLᵀ
  pivots `d_j = K_{j+1}/K_j` telescopes to `K_k / K_0`. Compiles clean, no `sorry`.
  NB (honest scope): this mechanizes only the algebraic telescoping IDENTITY `∏ d_j = K_k/K_0`
  under `K_j ≠ 0`. The substance of the round-2 correction — the BIT-SIZE claim
  `bit(K_k) = Θ(N)` (so fraction-free is Θ(N), not Θ(N²)) — is analytic and NOT mechanized here.

- `SProcedureSoundness.lean` (matrix bridge) — `IsPSD M := ∀v, 0 ≤ v ⬝ᵥ M *ᵥ v` (the exact
  form the band-`LDLᵀ` checker establishes) and `sprocedure_matrix_sound`: `IsPSD M` + the
  lifting identity `(lift x)ᵀ M (lift x) = f0 x − Σλᵢ fᵢ x − t` ⇒ `t ≤ J*`. Compiles clean,
  no `sorry`. This closes the `M⪰0 ⇒ soundness` step at the matrix level.

- `LiftingIdentity.lean` — `lifting2`: the concrete bordered 2×2 `M(λ,t)`'s quadratic form
  at `(x,1)` equals the minorant `f0 x − λ f1 x − t` (proved by `simp`+`ring`); and
  `scalar_sprocedure_sound`: `IsPSD M ⇒ t ≤ f0` on the feasible set, **no lifting hypothesis**.
  This is END-TO-END soundness for the scalar single keep-out (the shipped `bracket` primitive)
  at the matrix level. Compiles clean, no `sorry`.

- `CdtLifting.lean` — `cdt_lifting` + `cdt_sprocedure_sound`: the same, end-to-end, for the
  scalar TWO-keep-out (Celis–Dennis–Tapia / adversarial) case: PSD of the aggregate bordered
  matrix ⇒ `t` lower-bounds `f0` on the doubly-feasible set. Compiles clean, no `sorry`.

- `GeneralLifting.lean` — `bordered_lifting`: the GENERAL-`n` lifting. For the bordered
  `(n+1)×(n+1)` matrix from Hessian `P`, half-linear `b`, corner `c`, the quadratic form at
  `snoc x 1` equals `xᵀ P x + 2 bᵀ x + c` (any `n`). Removes the scalar restriction. Now also
  `general_sprocedure_sound`: composing the general-`n` lifting with weak duality gives the
  END-TO-END S-procedure lower bound for the general block-banded QCQP — `IsPSD(bordered P b c)`
  + `λ≥0` + the form-identity `⇒ t ≤ f0` on the feasible set (any `n`, any `m`). No `sorry`.

- `PosSemidefBridge.lean` — `posSemidef_iff_isPSD`: the `IsPSD` predicate decided throughout this
  development coincides with Mathlib's STANDARD `Matrix.PosSemidef` for symmetric matrices over ℚ
  (`M.PosSemidef ↔ Mᵀ = M ∧ IsPSD M`). So the checker decides the genuine, library-standard notion
  of positive semidefiniteness — not a bespoke predicate. No `sorry`.

- `BlockSchur.lean` — `block_schur_complement` (+ `block_form`, `block_pivot_psd`): the BLOCK
  (`d>1` vector-state) Riccati elimination step — the real RPOD case where each stage pivot is a
  `d×d` block, not a scalar. For `[[A,B],[Bᵀ,D]]` (via `Matrix.fromBlocks` over a `⊕` index) with
  `A` invertible (two-sided inverse given as data, avoiding the nonsingular-inverse module),
  `IsPSD` ⇒ block Schur complement `D − Bᵀ Ai B` is PSD, by minimizing the block form over the
  eliminated block coordinate `x = −Ai(Bw)`. (Correctness of the scalar sweep already covers this
  via `ldlt_complete`; the block version validates the genuine vector-state Riccati step directly.)
  No `sorry`.

- `SchurStep.lean` — `schur_complement_psd` (+ `headBordered_form`, `headBordered_pivot_nonneg`):
  the inductive ENGINE of the checker COMPLETENESS leg (the hard direction, converse of
  `ldlt_isPSD`). One step of symmetric banded/Schur elimination on `[[a,bᵀ],[b,C]]`: `IsPSD` ⇒
  pivot `a ≥ 0`, and (with `a>0`) the Schur complement `C − a⁻¹bbᵀ` is PSD — proved by
  minimizing the quadratic form over the eliminated coordinate `x = −bᵀw/a`. Iterating this over
  the head coordinate IS the exact band-`LDLᵀ` sweep, so this is the mathematically load-bearing
  content of `IsPSD ⇒ LᵀDL with pivots ≥ 0`. Also `headBordered_border_zero`: the degenerate
  `a=0` case's key fact (`IsPSD` ∧ zero pivot ⇒ zero border `b`), so the induction passes through
  rank-deficient steps. And the decomposition is made EXPLICIT: `rank1_isPSD` (`vvᵀ` PSD),
  `embed_isPSD` (`[[0,0],[0,S]]` PSD iff `S` is), and `headBordered_decomp` (the equality
  `[[a,bᵀ],[b,C]] = a⁻¹(a;b)(a;b)ᵀ + [[0,0],[0,C−a⁻¹bbᵀ]]`). Iterating that equality peels one
  nonnegative-weighted rank-1 per step, so `M = ∑ⱼ dⱼvⱼvⱼᵀ` — the `LDLᵀ`. **This is now fully
  assembled**: `ldlt_complete` proves `IsPSD M ⇒ ∃ dⱼ ≥ 0, vⱼ, M = ∑ⱼ dⱼ vⱼvⱼᵀ` by induction on
  the size (peel `a⁻¹(a;b)(a;b)ᵀ`, recurse on the Schur complement; degenerate `a=0` decouples).
  With `ldlt_isPSD` this makes **both directions of the checker's PSD criterion machine-checked**.
  No `sorry`.

- `Certificate.lean` (concrete instances) — `keepout_certified` (1-D) and `disk_certified` (2-D):
  fully worked, machine-checked exact-ℚ certificates. **1-D:** `min (x−½)² s.t. |x|≥1`, cert
  (`λ=½`, `t=¼`, border `[[½,−½],[−½,½]]`, `x̄=1`) ⇒ `x̄=1` optimal, `J*=¼`. **2-D (vector-state,
  RPOD-flavored):** `min (x−½)²+y² s.t. x²+y²≥1` (disk keep-out), cert (`λ=½`, `t=¼`, bordered
  `3×3` PSD matrix, `(1,0)`) ⇒ `(1,0)` optimal, `J*=¼`. Both exercise the whole pipeline on real
  rationals. No `sorry`.

- `Certificate.lean` — `certified_optimum_pipeline` (+ `accept_isPSD`): the HEADLINE capstone,
  the whole pipeline as ONE self-contained theorem. From a single accepted exact certificate
  (`M = LᵀDL`, pivots ≥ 0) + the S-procedure form-identity + a KKT witness `x̄`, it returns all
  three: `t` is a certified lower bound on the feasible set (`t ≤ J*`), `x̄` attains the global
  minimum, and `t ≤ f0 x̄` — i.e. `f0 x̄ = J*` is machine-checked. Chains accept⇒PSD, weak
  duality, and attainment inline (self-contained). Over any `LinearOrderedField`, no `sorry`.

- `GershgorinPSD.lean` — `isPSD_of_wdd`: general `n×n` weak diagonal dominance ⇒ PSD (paper
  Lemma dd / Gershgorin, at FULL generality — generalizes the 2×2 `DiagDominance.lean`). For
  symmetric `M` with every row `Mᵢᵢ ≥ ∑_{j≠i}|Mᵢⱼ|`, the quadratic form is bounded below by
  `∑ᵢ vᵢ²(Mᵢᵢ − ∑_{j≠i}|Mᵢⱼ|) ≥ 0` via a per-pair complete-the-square and a symmetry reindex
  `∑ᵢ∑ⱼ|Mᵢⱼ|vᵢ² = ∑ᵢ∑ⱼ|Mᵢⱼ|vⱼ²`. Over any `LinearOrderedField`, no `sorry`. (Needs the `ring`
  tactic — `lake build Mathlib.Tactic.Ring` once for the slice.)

- `StructuralPSD.lean` — the STRUCTURAL BACKBONE: the PSD matrices form a convex cone closed
  under congruence. `isPSD_add` (aggregate certificates), `isPSD_smul` (nonneg multiplier),
  `isPSD_zero`, `isPSD_nonneg_combo` (`M0 + Σλᵢ Mᵢ` PSD for `λ≥0`, `Mᵢ` PSD — the S-procedure
  aggregate), and `isPSD_congr` (`IsPSD M ⇒ IsPSD (SᵀMS)`). This is *why the whole method is
  sound*: S-procedure aggregation is a nonneg combination (stays in the cone), and every exact
  Gaussian/Schur elimination step is a congruence `M ↦ SᵀMS` (preserves PSD), so the exact band
  sweep never leaves the cone — `ldlt_isPSD` is just `isPSD_congr` at `S=L`, `M=diagonal(pivots)`.
  Over any `LinearOrderedField`, no `sorry`.

- `LDLTSound.lean` — `ldlt_isPSD` (+ `isPSD_diagonal`): the SOUNDNESS leg of the exact `LDLᵀ`
  checker (paper Lemma band, accept ⇒ PSD). If `M = Lᵀ D L` with `D = diagonal d` and every
  pivot `d i ≥ 0`, then `M` is PSD — because `vᵀMv = (Lv)ᵀ D (Lv) = Σᵢ dᵢ(Lv)ᵢ² ≥ 0` (proved
  via `mulVec_mulVec`/`dotProduct_mulVec`/`vecMul_transpose`, full generality, any `n`). This is
  exactly the guarantee the exact-ℚ checker relies on: its accept condition (all exact pivots
  nonnegative) implies the certified PSD verdict. No `sorry`. (The converse — PSD ⇒ such a
  factorization exists, and that band elimination computes it — is the harder completeness leg.)

- `ClosureOptimum.lean` — `certified_global_optimum` (+ `bracket_exact`): the ATTAINMENT /
  closure leg (paper Thm 1, second leg). A witness `x̄` that minimizes the Lagrangian
  `L = f0 − Σλᵢ fᵢ` (= stationarity when `L` is convex) and satisfies complementary slackness
  (`Σλᵢ fᵢ x̄ = 0`) is a GLOBAL optimum (`f0 x̄ ≤ f0 x` on feasibles); with the lower-bound leg
  the bracket closes to the exact value `f0 x̄ = J*`. Over any `LinearOrderedField`, no `sorry`.

- `MatrixDare.lean` — `dareOrbit_constant` (+ `dareStep_fixed`): the MATRIX DARE
  horizon-independence (paper Lemma dare, matrix version). Seeding the block-Schur/Riccati
  recursion `S ↦ diag − corr S` at a stationary `S*` with the DARE-generating running cost
  `diag = S* + corr S*` (in the matrix setting `corr S = GᵀS⁻¹G`) makes every pivot `= S*`.
  The pivot is constant by pure additive cancellation, so the lemma holds at the `AddCommGroup`
  level (instantiable at `d×d` matrices). Compiles clean — `#print axioms` = **no axioms at all**.

- `DiagDominance.lean` — `wdd_form_nonneg` + `isPSD_of_wdd_2x2`: weak diagonal dominance
  (`a,c ≥ |b|`) ⇒ the `2×2` symmetric form / matrix is PSD (paper Lemma dd's Gershgorin),
  proved from the quadratic form via SOS squares (`nlinarith`), no eigenvalue theory. The
  chain's interior row `[-1, 2+w−λ, -1]`. Compiles clean, no `sorry`.

## Remaining Lean gap to a fully end-to-end mechanized checker
Both DIRECTIONS of the checker's PSD criterion are now FULLY mechanized: SOUNDNESS
(`LDLTSound.ldlt_isPSD`: `LᵀDL` with pivots ≥ 0 ⇒ `IsPSD`) and COMPLETENESS
(`SchurStep.ldlt_complete`: `IsPSD` ⇒ `M = ∑ⱼ dⱼvⱼvⱼᵀ` with `dⱼ ≥ 0`). So the abstract checker
criterion `accept ⇔ PSD` is machine-checked in both directions, no `sorry`. The ONLY remaining
gap to a *fully* end-to-end verified checker is the ALGORITHMIC refinement: proving that the
concrete band elimination in `podium.verify.riccati` computes exactly such a factorization (an
implementation-correctness / refinement task, not new mathematics). Everything with mathematical
content — lifting (scalar/CDT/general-`n`), `IsPSD ⇒ minorant ⇒ t≤J*`, closure ⇒ global optimum,
`pivots≥0 ⇒ PSD`, and `IsPSD ⇒ LDLᵀ` (both checker directions) — is mechanized, no `sorry`.

**Second gap:** the paragraph above is true of
`LDLTSound`/`SchurStep`'s OWN `IsPSD`, but that `IsPSD` is not itself machine-linked to
Mathlib's `Matrix.PosSemidef`. `PosSemidefBridge.lean`'s `posSemidef_iff_isPSD` proves
`M.PosSemidef ↔ IsPSD M` for a THIRD, separately-declared `IsPSD` — `LDLTSound.lean`,
`SchurStep.lean`, and `PosSemidefBridge.lean` each declare a top-level `def IsPSD` with no
`namespace` wrapper, so no single file can `import` more than one of them (Lean rejects the
redeclaration). The three formulas are definitionally identical after specializing to `ℚ`,
confirmed by inspection, but nothing compiles that chain end-to-end: `PosSemidefBridge.lean`
is an isolated leaf, never imported by any other corpus file. Fixing this for real needs a
namespace/rename refactor of at least one of the three declarations — which touches files
already cited by exact identifier in `papers/riccati-cert/main.tex` (`ldlt_isPSD`,
`isPSD_iff_ldlt`), so treat it as its own reviewed change, not a drive-by patch. Tracked as a
Lean formalization worklist item; `main.tex`'s corresponding sentence was narrowed
to disclose the gap honestly rather than overstate it.

## How to build
All 18 top-level files (19 with `note/PodiumNote.lean`) type-check clean against the
pinned Mathlib project that `setup.sh` builds into `.leanwork/`, from the pins in
`ci/{lakefile.toml,lake-manifest.json,lean-toolchain}`.

```
bash setup.sh        # once: toolchain + pinned Mathlib (via lake exe cache get)
bash check-all.sh    # regression: PASS/FAIL per file, then ALL PASS
```

`check-all.sh` also surfaces `#print axioms` for every theorem — expect only
propext/Classical.choice/Quot.sound, never `sorryAx`. A single file:

```
env -C .leanwork lake env lean <File>.lean
```

Set `LEANWORK` to reuse a Mathlib project built elsewhere, `LAKE` for a `lake` binary
that is not on `PATH`.
