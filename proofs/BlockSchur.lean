import Mathlib.Data.Matrix.Block
import Mathlib.Data.Matrix.Mul
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-
Block Schur complement (the `d>1` vector-state Riccati elimination step — the real RPOD case,
where each stage pivot is a `d×d` block, not a scalar).  For the block matrix `[[A, B],[Bᵀ, D]]`
with `A ≻ 0` (given an inverse `Ai`), `IsPSD ⇒ D − Bᵀ Ai B ⪰ 0` — obtained by minimizing the
block quadratic form over the eliminated block coordinate `x = −Ai (B w)`.  We take the inverse
as data (a two-sided inverse of `A`) to avoid the nonsingular-inverse machinery.  Over any
`LinearOrderedField`, no `sorry`.
-/

open Matrix
variable {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m]
  {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- PSD predicate over an arbitrary finite index. -/
def IsPSD {ι : Type*} [Fintype ι] (M : Matrix ι ι K) : Prop := ∀ v : ι → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- Block quadratic form: `[x;w]ᵀ [[A,B],[C,D]] [x;w] = xᵀAx + xᵀBw + wᵀCx + wᵀDw`. -/
theorem block_form (A : Matrix m m K) (B : Matrix m n K) (C : Matrix n m K) (D : Matrix n n K)
    (x : m → K) (w : n → K) :
    (Sum.elim x w) ⬝ᵥ (fromBlocks A B C D) *ᵥ (Sum.elim x w)
      = x ⬝ᵥ A *ᵥ x + x ⬝ᵥ B *ᵥ w + w ⬝ᵥ C *ᵥ x + w ⬝ᵥ D *ᵥ w := by
  rw [fromBlocks_mulVec]
  simp only [dotProduct, Fintype.sum_sum_type, Sum.elim_inl, Sum.elim_inr,
    Sum.elim_comp_inl, Sum.elim_comp_inr, Pi.add_apply, mulVec, mul_add]
  simp only [Finset.sum_add_distrib]
  ring

/-- Block pivot PSD: `IsPSD [[A,B],[Bᵀ,D]]` ⇒ the pivot block `A` is PSD (test `[x;0]`). -/
theorem block_pivot_psd (A : Matrix m m K) (B : Matrix m n K) (D : Matrix n n K)
    (hA : IsPSD (fromBlocks A B Bᵀ D)) : IsPSD A := by
  intro x
  have hval := hA (Sum.elim x 0)
  rw [block_form] at hval
  simpa using hval

/-- **Block Schur complement preserves PSD (the `d>1` Riccati step).** For the symmetric block
matrix `[[A, B],[Bᵀ, D]]` with `A` invertible (given a two-sided inverse `Ai`), `IsPSD` ⇒ the
block Schur complement `D − Bᵀ Ai B` is PSD.  Proof: minimize the block quadratic form over the
eliminated block coordinate `x = −Ai(Bw)`; the minimum value is `wᵀ(D − Bᵀ Ai B)w`. -/
theorem block_schur_complement
    (A : Matrix m m K) (B : Matrix m n K) (D : Matrix n n K) (Ai : Matrix m m K)
    (hA : IsPSD (fromBlocks A B Bᵀ D)) (hAi : A * Ai = 1) :
    IsPSD (D - Bᵀ * Ai * B) := by
  intro w
  set u := B *ᵥ w with hu
  set xstar : m → K := -(Ai *ᵥ u) with hx
  have adj : ∀ z : m → K, w ⬝ᵥ Bᵀ *ᵥ z = u ⬝ᵥ z := by
    intro z; rw [dotProduct_mulVec, vecMul_transpose, ← hu]
  have hAx : A *ᵥ xstar = -u := by
    rw [hx, mulVec_neg, mulVec_mulVec, hAi, one_mulVec]
  have hT1 : xstar ⬝ᵥ A *ᵥ xstar = u ⬝ᵥ (Ai *ᵥ u) := by
    rw [hAx, hx, neg_dotProduct, dotProduct_neg, neg_neg, dotProduct_comm]
  have hT2 : xstar ⬝ᵥ B *ᵥ w = -(u ⬝ᵥ (Ai *ᵥ u)) := by
    rw [← hu, hx, neg_dotProduct, dotProduct_comm]
  have hT3 : w ⬝ᵥ Bᵀ *ᵥ xstar = -(u ⬝ᵥ (Ai *ᵥ u)) := by
    rw [adj, hx, dotProduct_neg]
  have hRHS : w ⬝ᵥ (D - Bᵀ * Ai * B) *ᵥ w = w ⬝ᵥ D *ᵥ w - u ⬝ᵥ (Ai *ᵥ u) := by
    rw [sub_mulVec, dotProduct_sub]
    congr 1
    rw [← mulVec_mulVec, ← mulVec_mulVec, adj, ← hu]
  have hval : (Sum.elim xstar w) ⬝ᵥ (fromBlocks A B Bᵀ D) *ᵥ (Sum.elim xstar w)
      = w ⬝ᵥ (D - Bᵀ * Ai * B) *ᵥ w := by
    rw [block_form, hT1, hT2, hT3, hRHS]; ring
  rw [← hval]
  exact hA (Sum.elim xstar w)

/-- Block elimination as an additive decomposition: `[[A,B],[Bᵀ,D]]` splits into the pivot part
`[[A,B],[Bᵀ,BᵀAiB]]` plus the embedded Schur complement `[[0,0],[0, D−BᵀAiB]]`. -/
theorem block_decomp (A : Matrix m m K) (B : Matrix m n K) (D : Matrix n n K) (Ai : Matrix m m K) :
    fromBlocks A B Bᵀ D
      = fromBlocks A B Bᵀ (Bᵀ * Ai * B) + fromBlocks 0 0 0 (D - Bᵀ * Ai * B) := by
  ext i j
  cases i <;> cases j <;>
    simp [Matrix.fromBlocks, Matrix.add_apply, Matrix.sub_apply]

#print axioms block_form
#print axioms block_pivot_psd
#print axioms block_schur_complement
#print axioms block_decomp
