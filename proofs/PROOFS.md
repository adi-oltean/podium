# Proofs — human-readable arguments behind the Lean formalization

The mathematics behind every mechanized lemma in this directory, written so the ideas can be
checked on paper independently of Lean. Each entry names the Lean theorem and file. The
arguments are elementary; the value of the mechanization is a tolerance-free, `sorry`-free
trusted path, not new mathematics (see `STATUS.md`).

Notation. `M ⪰ 0` (PSD) means `∀ v, vᵀMv ≥ 0` (`IsPSD`). The S-procedure / bordered "lifting"
sends a state `x` to the homogenized vector `[x;1]`. `J* = inf{ f0(x) : f_k(x) ≥ 0 ∀k }` is the
optimal value of the (nonconvex) QCQP.

---

## 1. Soundness — a passing certificate lower-bounds the optimum

### 1.1 S-procedure weak duality (`note_thm_sound`, `note_thm_multi`; `sprocedure_*`)
**Claim.** If `λ_k ≥ 0` and the minorant `f0(x) − Σ_k λ_k f_k(x) ≥ t` holds for all `x`, then
`t ≤ f0(x)` for every feasible `x` (so `t ≤ J*`).
**Proof.** For feasible `x`, each `f_k(x) ≥ 0` and `λ_k ≥ 0`, so `Σ_k λ_k f_k(x) ≥ 0`. Hence
`f0(x) ≥ f0(x) − Σ_k λ_k f_k(x) ≥ t`. ∎

### 1.2 Lifting identity (`bordered_lifting`, general `n`; `lifting2`, `cdt_lifting`)
**Claim.** For the bordered matrix `M = [[P, b],[bᵀ, c]]`, `[x;1]ᵀ M [x;1] = xᵀPx + 2bᵀx + c`.
**Proof.** Expand the block product: `[x;1]ᵀ M [x;1] = xᵀPx + xᵀb + bᵀx + c = xᵀPx + 2bᵀx + c`
(using `bᵀx = xᵀb`). ∎
With `P = P0 − Σλ_k P_k`, `b = ½(q0 − Σλ_k q_k)`, `c = r0 − Σλ_k r_k − t`, the RHS is exactly
the minorant `f0(x) − Σλ_k f_k(x) − t`; so `M ⪰ 0` gives the hypothesis of §1.1.

### 1.3 Matrix soundness, end-to-end (`sprocedure_matrix_sound`, `general_sprocedure_sound`,
`certified_optimum_pipeline`)
**Claim.** `M ⪰ 0` + `λ ≥ 0` + the lifting/form identity ⇒ `t ≤ f0` on the feasible set.
**Proof.** `M ⪰ 0` applied to `[x;1]` gives `[x;1]ᵀM[x;1] ≥ 0`; by §1.2 this is the minorant
`≥ 0`; apply §1.1. ∎

### 1.4 Attainment / closure — the bracket closes (`certified_global_optimum`, `bracket_exact`)
**Claim.** If additionally the witness `x̄` minimizes the Lagrangian `L(x) = f0(x) − Σλ_k f_k(x)`
(= stationarity when `L` is convex) and complementary slackness `Σλ_k f_k(x̄) = 0` holds, then
`x̄` is a global optimum: `f0(x̄) ≤ f0(x)` for all feasible `x`. Combined with §1.3, `f0(x̄) = J*`.
**Proof.** For feasible `x`: `f0(x̄) = f0(x̄) − Σλ_k f_k(x̄)` (compl. slackness) `= L(x̄) ≤ L(x)`
(x̄ minimizes L) `= f0(x) − Σλ_k f_k(x) ≤ f0(x)` (feasibility + `λ≥0`). ∎

---

## 2. The checker's PSD criterion — both directions

The trusted checker decides `M ⪰ 0` by exact-ℚ `LDLᵀ`: accept iff all pivots `d_i ≥ 0`.

### 2.0 The criterion as one iff (`isPSD_iff_ldlt`) — the decision-procedure spec
**Claim.** Symmetric `M ⪰ 0` **iff** `M = Σⱼ dⱼvⱼvⱼᵀ` for some `dⱼ ≥ 0` — i.e. the exact checker's
`accept` (a nonnegative-pivot `LDLᵀ` is found) is equivalent to PSD. Forward = §2.5
(`ldlt_complete`); backward = §2.1 idea (each `dⱼvⱼvⱼᵀ ⪰ 0`, sum of PSD is PSD).

### 2.1 Soundness: accept ⇒ PSD (`ldlt_isPSD`, `isPSD_diagonal`; `accept_isPSD`)
**Claim.** If `M = LᵀDL` with `D = diag(d)`, `d_i ≥ 0`, then `M ⪰ 0`.
**Proof.** `vᵀMv = vᵀLᵀDLv = (Lv)ᵀD(Lv) = Σ_i d_i (Lv)_i² ≥ 0`. ∎

### 2.2 Completeness — the elimination step (`schur_complement_psd`, `headBordered_pivot_nonneg`)
Write `M = [[a, bᵀ],[b, C]]` (head pivot `a`). Quadratic form `q(x,w) = a x² + 2x(bᵀw) + wᵀCw`.
**Pivot ≥ 0.** `M ⪰ 0` at `[1;0]` gives `q(1,0) = a ≥ 0`.
**Schur complement PSD (a > 0).** For any tail `w`, `q(·,w)` is a quadratic in `x` with positive
leading coeff, minimized at `x* = −bᵀw/a` with value `q(x*,w) = wᵀCw − (bᵀw)²/a = wᵀ(C − a⁻¹bbᵀ)w`.
Since `M ⪰ 0`, `q(x*,w) ≥ 0`, i.e. `C − a⁻¹bbᵀ ⪰ 0`. ∎

### 2.2b Block Schur complement — the `d>1` Riccati step (`block_schur_complement`, `block_pivot_psd`)
**Claim.** For `[[A, B],[Bᵀ, D]]` with `A` invertible (inverse `Ai`), `IsPSD ⇒ D − Bᵀ Ai B ⪰ 0`,
and the pivot block `A ⪰ 0`.
**Proof.** Same minimize-out-the-pivot idea, block-valued. Block form
`q(x,w) = xᵀAx + 2xᵀBw + wᵀDw`; minimize over the block coordinate `x` at `x* = −Ai(Bw)`, giving
value `wᵀDw − (Bw)ᵀAi(Bw) = wᵀ(D − BᵀAiB)w ≥ 0`. Pivot: evaluate at `[x;0]` → `xᵀAx ≥ 0`.
(Scalar `ldlt_complete` already covers correctness; this validates the genuine vector-state step.)

### 2.3 Completeness — the degenerate pivot (`headBordered_border_zero`)
**Claim.** `[[0, bᵀ],[b, C]] ⪰ 0 ⇒ b = 0`.
**Proof.** With `a = 0`, `q(x,w) = 2x(bᵀw) + wᵀCw` is *linear* in `x`. If `bᵀw = c ≠ 0` for some
`w`, set `x = −(wᵀCw + 1)/(2c)` to get `q = −1 < 0`, contradicting `M ⪰ 0`. So `bᵀw = 0` for all
`w`, hence `b = 0`. ∎

### 2.4 Explicit rank-1 + Schur decomposition (`headBordered_decomp`, `rank1_isPSD`, `embed_isPSD`)
**Claim (a ≠ 0).** `[[a,bᵀ],[b,C]] = a⁻¹ (a;b)(a;b)ᵀ + [[0,0],[0, C − a⁻¹bbᵀ]]`.
**Proof.** `(a;b)(a;b)ᵀ = [[a², a bᵀ],[a b, bbᵀ]]`; scaling by `a⁻¹` gives `[[a, bᵀ],[b, a⁻¹bbᵀ]]`;
adding `[[0,0],[0, C − a⁻¹bbᵀ]]` restores the `(2,2)` block to `C`. ∎
Supporting: `vvᵀ ⪰ 0` since `wᵀ(vvᵀ)w = (vᵀw)² ≥ 0` (`rank1_isPSD`); `[[0,0],[0,S]] ⪰ 0 ⇔ S ⪰ 0`
since its form is `wᵀSw` regardless of the head coordinate (`embed_isPSD`).

### 2.5 Completeness — the assembly (`ldlt_complete`; MECHANIZED, no `sorry`)
**Claim.** `M ⪰ 0` (symmetric, `n×n`) ⇒ `M = Σ_{j} d_j v_j v_jᵀ` with `d_j ≥ 0` — i.e. `M = VᵀDV`,
the `LDLᵀ` with nonnegative pivots (the converse of §2.1).
**Proof (induction on `n`).** `n = 0`: empty sum. `n+1`: pivot `a = M₀₀ ≥ 0` by §2.2.
- If `a > 0`: by §2.4, `M = a⁻¹ v vᵀ + embed(S)` with `v = (a;b)`, `S = C − a⁻¹bbᵀ`. By §2.2,
  `S ⪰ 0`, so by induction `S = Σ_j d_j' w_j w_jᵀ`, `d_j' ≥ 0`. Embedding (a linear map on the
  tail) gives `embed(S) = Σ_j d_j' (0;w_j)(0;w_j)ᵀ`. Then
  `M = a⁻¹ v vᵀ + Σ_j d_j' (0;w_j)(0;w_j)ᵀ` — a nonneg combination of rank-1s, with `a⁻¹ ≥ 0`.
- If `a = 0`: by §2.3, `b = 0`, so `M = embed(C)` with `C ⪰ 0` (restrict `M ⪰ 0` to the tail).
  By induction `C = Σ_j d_j' w_j w_jᵀ`; embed. The eliminated coordinate contributes pivot `0`. ∎
**Status.** MECHANIZED as `ldlt_complete` (`SchurStep.lean`), no `sorry` — the induction above
transcribed directly, carrying the family `{(d_j, v_j)}` as `Fin.cons`-extended tuples and the
`Fin n ↪ Fin(n+1)` embedding as `Fin.cons 0`. With §2.1 (`ldlt_isPSD`) this makes **both
directions** of the checker's PSD criterion `accept ⇔ PSD` machine-checked. The only remaining
gap to a fully end-to-end verified checker is the algorithmic refinement (that the concrete band
elimination computes such a factorization) — implementation correctness, not new mathematics.

---

## 3. Structural cone — why aggregation and elimination are sound (`StructuralPSD`)
**Claims / proofs (all one line on the quadratic form).**
- `isPSD_add`: `vᵀ(M+N)v = vᵀMv + vᵀNv ≥ 0`.
- `isPSD_smul` (`c ≥ 0`): `vᵀ(cM)v = c(vᵀMv) ≥ 0`.
- `isPSD_congr`: `vᵀ(SᵀMS)v = (Sv)ᵀM(Sv) ≥ 0`.
- `isPSD_nonneg_combo`: `M0 + Σλ_i M_i ⪰ 0` for `λ_i ≥ 0`, `M_i ⪰ 0`, by add + smul.
These say the PSD matrices form a convex cone closed under congruence: S-procedure aggregation is a
nonneg combination (stays in the cone), and each exact elimination step is a congruence
`M ↦ SᵀMS` (preserves PSD) — so the exact band sweep never leaves the cone.

---

## 4. Diagonal dominance ⇒ PSD, general `n` (`isPSD_of_wdd`, Gershgorin)
**Claim.** Symmetric `M` with `M_ii ≥ Σ_{j≠i} |M_ij|` for every row is PSD.
**Proof.** Bound the form below by completing the square on each off-diagonal pair:
`vᵀMv = Σ_i M_ii v_i² + Σ_{i≠j} M_ij v_i v_j`. For each pair, `M_ij v_i v_j ≥ −½|M_ij|(v_i²+v_j²)`
(since `|M_ij|(v_i ± v_j)² ≥ 0`). Summing and using the symmetry reindex
`Σ_i Σ_{j≠i}|M_ij|v_i² = Σ_i Σ_{j≠i}|M_ij|v_j²` (from `|M_ij| = |M_ji|`), the off-diagonal penalty
is `Σ_i v_i² Σ_{j≠i}|M_ij|`. Hence `vᵀMv ≥ Σ_i v_i²(M_ii − Σ_{j≠i}|M_ij|) ≥ 0` by dominance. ∎

---

## 5. Horizon-independence and bit-size

### 5.1 Scalar DARE fixed point (`dare_fixed_point`, `dare_orbit_constant`)
On the Pell locus `c = m + 1/m`, the seeded Riccati/continued-fraction recursion `s ↦ c − 1/s`
has fixed point `m` (`m = c − 1/m` ⇔ `m² − cm + 1 = 0`), so seeding at `m` gives a constant orbit —
horizon-independent pivots.

### 5.2 Matrix DARE (`dareOrbit_constant`)
The block-Schur step `S ↦ diag − corr(S)` (with `corr(S) = GᵀS⁻¹G`) seeded at a stationary `S*`
with running cost `diag = S* + corr(S*)` is constant: `step(S*) = (S* + corr(S*)) − corr(S*) = S*`
— pure additive cancellation, so it holds at the `AddCommGroup` level.

### 5.3 Telescoping continuant (`pivotProd_telescopes`)
The tridiagonal `LDLᵀ` pivots `d_j = K_{j+1}/K_j` (ratios of continuants) telescope:
`∏_{j<k} d_j = K_k / K_0`. (This is the algebraic identity; the *bit-size* claim `bit(K_k) = Θ(N)`
is analytic and not mechanized — see `STATUS.md`.)

---

## 6. Not mechanized (analytic — stated, not proved here)
- Recovery *rates* (`thm:recovery` `O(1/D²)`, `thm:hard` linear): Taylor / one-sided derivative
  of the dual value; need analysis. The *soundness* content is §1 (mechanized); only the rates are
  analytic.
- General bit-complexity beyond constructed families: open (Bienstock-hard), see the paper.
