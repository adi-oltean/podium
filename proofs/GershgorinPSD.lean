import Mathlib.Data.Matrix.Mul
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-
General `n×n` weak diagonal dominance ⇒ PSD (paper Lemma dd / Gershgorin, at full generality,
generalizing the 2×2 `DiagDominance.lean`).  For a SYMMETRIC matrix `M` with every row weakly
diagonally dominant (`M i i ≥ ∑_{j≠i} |M i j|`), `M` is PSD.  Proof: bound the quadratic form
below by completing the square on each off-diagonal pair,
  `vᵀMv ≥ ∑ᵢ vᵢ²(Mᵢᵢ − ∑_{j≠i}|Mᵢⱼ|) ≥ 0`,
the row-excess being nonnegative by dominance.  The cross terms are handled by
`|Mᵢⱼ|(vᵢ²+vⱼ²) + 2 Mᵢⱼ vᵢvⱼ = |Mᵢⱼ|(vᵢ + sgn·vⱼ)² ≥ 0` and a symmetry reindex
`∑ᵢ∑ⱼ|Mᵢⱼ|vᵢ² = ∑ᵢ∑ⱼ|Mᵢⱼ|vⱼ²`.  Over any `LinearOrderedField`, no `sorry`.
-/

open Matrix Finset
variable {n : ℕ} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

def IsPSD (M : Matrix (Fin n) (Fin n) K) : Prop := ∀ v : Fin n → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- **Weak diagonal dominance ⇒ PSD (general `n`).** Symmetric `M` with `Mᵢᵢ ≥ ∑_{j≠i}|Mᵢⱼ|`
for every row is positive semidefinite. -/
theorem isPSD_of_wdd (M : Matrix (Fin n) (Fin n) K)
    (hsymm : ∀ i j, M i j = M j i)
    (hdd : ∀ i, ∑ j ∈ univ.erase i, |M i j| ≤ M i i) :
    IsPSD M := by
  intro v
  -- quadratic form as a double sum
  have hform : v ⬝ᵥ M *ᵥ v = ∑ i, ∑ j, v i * M i j * v j := by
    simp only [dotProduct, mulVec, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ =>
      (mul_assoc (v i) (M i j) (v j)).symm
  rw [hform]
  -- erase-sum = full-sum minus the removed term
  have herase : ∀ (F : Fin n → K) (i : Fin n),
      ∑ j ∈ univ.erase i, F j = (∑ j, F j) - F i := by
    intro F i
    have h := Finset.add_sum_erase univ F (Finset.mem_univ i)
    linarith [h]
  -- symmetry reindex on the full square: ∑ᵢ∑ⱼ |Mᵢⱼ| vᵢ² = ∑ᵢ∑ⱼ |Mᵢⱼ| vⱼ²
  have hfull : (∑ i, ∑ j, |M i j| * v i ^ 2) = (∑ i, ∑ j, |M i j| * v j ^ 2) := by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => by rw [hsymm y x]
  -- per-row lower bound: gᵢᵢ + Σ_{j≠i} gᵢⱼ  ≥  vᵢ Mᵢᵢ vᵢ − Σ_{j≠i} |Mᵢⱼ|(vᵢ²+vⱼ²)/2
  have hlow : ∀ i, v i * M i i * v i + ∑ j ∈ univ.erase i, (-(|M i j| * (v i ^ 2 + v j ^ 2) / 2))
      ≤ ∑ j, v i * M i j * v j := by
    intro i
    rw [← Finset.add_sum_erase univ (fun j => v i * M i j * v j) (Finset.mem_univ i)]
    apply add_le_add_right
    apply Finset.sum_le_sum
    intro j _
    have h1 : (0 : K) ≤ (|M i j| - M i j) * (v i - v j) ^ 2 :=
      mul_nonneg (by linarith [le_abs_self (M i j)]) (sq_nonneg _)
    have h2 : (0 : K) ≤ (|M i j| + M i j) * (v i + v j) ^ 2 :=
      mul_nonneg (by linarith [neg_abs_le (M i j)]) (sq_nonneg _)
    nlinarith [h1, h2]
  -- the summed lower bound is nonnegative
  have hnn : (0 : K) ≤ ∑ i, (v i * M i i * v i
      + ∑ j ∈ univ.erase i, (-(|M i j| * (v i ^ 2 + v j ^ 2) / 2))) := by
    -- rewrite each row's penalty sum, then use the reindex to collapse to the row excess
    have hrow : ∀ i, v i * M i i * v i
        + ∑ j ∈ univ.erase i, (-(|M i j| * (v i ^ 2 + v j ^ 2) / 2))
        = v i * M i i * v i
          - ((∑ j ∈ univ.erase i, |M i j| * v i ^ 2) / 2
             + (∑ j ∈ univ.erase i, |M i j| * v j ^ 2) / 2) := by
      intro i
      have hp : ∑ j ∈ univ.erase i, (-(|M i j| * (v i ^ 2 + v j ^ 2) / 2))
          = -((∑ j ∈ univ.erase i, |M i j| * v i ^ 2) / 2
               + (∑ j ∈ univ.erase i, |M i j| * v j ^ 2) / 2) := by
        rw [Finset.sum_div, Finset.sum_div, ← Finset.sum_add_distrib, ← Finset.sum_neg_distrib]
        exact Finset.sum_congr rfl fun j _ => by ring
      linarith [hp]
    simp_rw [hrow]
    -- Aerase = Berase by the reindex + herase
    have hAB : (∑ i, ∑ j ∈ univ.erase i, |M i j| * v j ^ 2)
             = (∑ i, ∑ j ∈ univ.erase i, |M i j| * v i ^ 2) := by
      have e1 : (∑ i, ∑ j ∈ univ.erase i, |M i j| * v j ^ 2)
          = (∑ i, ∑ j, |M i j| * v j ^ 2) - ∑ i, |M i i| * v i ^ 2 := by
        rw [← Finset.sum_sub_distrib]
        exact Finset.sum_congr rfl fun i _ => herase (fun j => |M i j| * v j ^ 2) i
      have e2 : (∑ i, ∑ j ∈ univ.erase i, |M i j| * v i ^ 2)
          = (∑ i, ∑ j, |M i j| * v i ^ 2) - ∑ i, |M i i| * v i ^ 2 := by
        rw [← Finset.sum_sub_distrib]
        exact Finset.sum_congr rfl fun i _ => herase (fun j => |M i j| * v i ^ 2) i
      rw [e1, e2, hfull]
    rw [Finset.sum_sub_distrib, Finset.sum_add_distrib]
    -- the two halves are equal, so the subtracted quantity = Aerase; bound each row by dominance
    have hcollapse :
        (∑ i, (∑ j ∈ univ.erase i, |M i j| * v i ^ 2) / 2)
        + (∑ i, (∑ j ∈ univ.erase i, |M i j| * v j ^ 2) / 2)
        = ∑ i, (∑ j ∈ univ.erase i, |M i j|) * v i ^ 2 := by
      rw [← Finset.sum_div, ← Finset.sum_div, hAB]
      have hx : (∑ i, ∑ j ∈ univ.erase i, |M i j| * v i ^ 2)
          = ∑ i, (∑ j ∈ univ.erase i, |M i j|) * v i ^ 2 :=
        Finset.sum_congr rfl fun i _ => by rw [Finset.sum_mul]
      rw [hx]; ring
    rw [hcollapse, ← Finset.sum_sub_distrib]
    apply Finset.sum_nonneg
    intro i _
    have : (∑ j ∈ univ.erase i, |M i j|) * v i ^ 2 ≤ M i i * v i ^ 2 :=
      mul_le_mul_of_nonneg_right (hdd i) (sq_nonneg _)
    nlinarith [this]
  calc (0 : K) ≤ ∑ i, (v i * M i i * v i
          + ∑ j ∈ univ.erase i, (-(|M i j| * (v i ^ 2 + v j ^ 2) / 2))) := hnn
    _ ≤ ∑ i, ∑ j, v i * M i j * v j := Finset.sum_le_sum fun i _ => hlow i

#print axioms isPSD_of_wdd
