import Mathlib.Data.Matrix.Mul
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.Linarith

/-
Checker-SOUNDNESS direction of exact `LDLᵀ` (paper Lemma band, the easy/soundness leg).
The exact-ℚ trusted checker ACCEPTS a symmetric matrix exactly when its `LDLᵀ` factorization
has all pivots `d i ≥ 0`.  For that acceptance to be trustworthy, pivots-nonneg must imply
PSD.  This file mechanizes exactly that, at full generality: if `M = Lᵀ D L` with `D` the
diagonal of the pivots and every pivot `≥ 0`, then `M` is PSD — because
`vᵀ M v = (L v)ᵀ D (L v) = Σᵢ dᵢ (Lv)ᵢ² ≥ 0`.  Over any `LinearOrderedField`, no `sorry`.
(The converse — PSD ⇒ such a factorization with pivots ≥ 0 exists — is the harder
completeness leg and is not needed for checker soundness.)
-/

open Matrix
variable {n : ℕ} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- PSD predicate the checker decides. -/
def IsPSD (M : Matrix (Fin n) (Fin n) K) : Prop := ∀ v : Fin n → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- A diagonal matrix with nonnegative entries is PSD (the base pivots case). -/
theorem isPSD_diagonal (d : Fin n → K) (hd : ∀ i, 0 ≤ d i) : IsPSD (diagonal d) := by
  intro v
  have hform : v ⬝ᵥ (diagonal d) *ᵥ v = ∑ i, d i * (v i * v i) := by
    simp only [dotProduct, mulVec_diagonal]
    exact Finset.sum_congr rfl (fun i _ => mul_left_comm (v i) (d i) (v i))
  rw [hform]
  exact Finset.sum_nonneg (fun i _ => mul_nonneg (hd i) (mul_self_nonneg _))

/-- **Checker soundness (LDLᵀ pivots ≥ 0 ⇒ PSD).** If `M = Lᵀ D L` with `D = diagonal d` and
every pivot `d i ≥ 0`, then `M` is PSD. This is the guarantee the exact-ℚ checker relies on:
its accept condition (all exact pivots nonnegative) implies the certified PSD verdict. -/
theorem ldlt_isPSD (L : Matrix (Fin n) (Fin n) K) (d : Fin n → K) (hd : ∀ i, 0 ≤ d i)
    (M : Matrix (Fin n) (Fin n) K) (hM : M = Lᵀ * diagonal d * L) : IsPSD M := by
  intro v
  have key : v ⬝ᵥ M *ᵥ v = (L *ᵥ v) ⬝ᵥ (diagonal d) *ᵥ (L *ᵥ v) := by
    rw [hM, Matrix.mul_assoc, ← Matrix.mulVec_mulVec, dotProduct_mulVec, vecMul_transpose,
      ← Matrix.mulVec_mulVec]
  rw [key]
  exact isPSD_diagonal d hd (L *ᵥ v)

#print axioms isPSD_diagonal
#print axioms ldlt_isPSD
