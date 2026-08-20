import Mathlib.LinearAlgebra.Matrix.PosDef

/-
Bridge: the `IsPSD` predicate used throughout this development (the exact form the band-`LDLᵀ`
checker decides) coincides with Mathlib's standard `Matrix.PosSemidef` for symmetric matrices
over ℚ.  So the checker decides the genuine, library-standard notion of positive semidefiniteness.
No `sorry`.
-/

open Matrix
variable {n : ℕ}

/-- `IsPSD` predicate as used in this development. -/
def IsPSD (M : Matrix (Fin n) (Fin n) ℚ) : Prop := ∀ v : Fin n → ℚ, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- Over ℚ, Mathlib's `Matrix.PosSemidef` is exactly "symmetric and `IsPSD`". -/
theorem posSemidef_iff_isPSD (M : Matrix (Fin n) (Fin n) ℚ) :
    M.PosSemidef ↔ Mᵀ = M ∧ IsPSD M := by
  rw [posSemidef_iff_dotProduct_mulVec]
  constructor
  · intro h
    refine ⟨?_, fun v => ?_⟩
    · have := h.1
      rwa [Matrix.IsHermitian, conjTranspose_eq_transpose_of_trivial] at this
    · have := h.2 v
      simpa using this
  · rintro ⟨hsym, hp⟩
    refine ⟨?_, fun v => ?_⟩
    · rw [Matrix.IsHermitian, conjTranspose_eq_transpose_of_trivial]; exact hsym
    · simpa using hp v

#print axioms posSemidef_iff_isPSD
