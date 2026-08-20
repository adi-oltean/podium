import Mathlib.Tactic.Linarith
import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Notation

/-
Diagonal dominance ⇒ PSD (paper Lemma dd / the benign-regime design rule), at the concrete
2×2 level — the chain's interior structure `[a, b; b, c]` (e.g. the row `[-1, 2+w-λ, -1]`).
Proved from the quadratic form via AM-GM squares (`nlinarith`), no eigenvalue theory.  No `sorry`.
-/

open Matrix
variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- Quadratic-form nonnegativity for a symmetric 2×2 weakly-diagonally-dominant matrix:
`a,c ≥ 0` and `a,c ≥ |b|` (as `±b ≤ a`, `±b ≤ c`) ⇒ the form is `≥ 0`. -/
theorem wdd_form_nonneg (a b c v0 v1 : K)
    (ha : 0 ≤ a) (hc : 0 ≤ c)
    (hab : b ≤ a) (hab' : -b ≤ a) (hcb : b ≤ c) (hcb' : -b ≤ c) :
    0 ≤ a * v0 ^ 2 + 2 * b * v0 * v1 + c * v1 ^ 2 := by
  rcases le_or_gt 0 b with hb | hb
  · -- b ≥ 0 :  form = b(v0+v1)² + (a−b)v0² + (c−b)v1²
    nlinarith [mul_nonneg hb (sq_nonneg (v0 + v1)),
      mul_nonneg (by linarith : (0:K) ≤ a - b) (sq_nonneg v0),
      mul_nonneg (by linarith : (0:K) ≤ c - b) (sq_nonneg v1)]
  · -- b < 0 :  form = (−b)(v0−v1)² + (a+b)v0² + (c+b)v1²
    nlinarith [mul_nonneg (by linarith : (0:K) ≤ -b) (sq_nonneg (v0 - v1)),
      mul_nonneg (by linarith : (0:K) ≤ a + b) (sq_nonneg v0),
      mul_nonneg (by linarith : (0:K) ≤ c + b) (sq_nonneg v1)]

/-- The same as PSD of the 2×2 matrix `!![a,b; b,c]` (`IsPSD := ∀ v, 0 ≤ vᵀMv`). -/
theorem isPSD_of_wdd_2x2 (a b c : K)
    (ha : 0 ≤ a) (hc : 0 ≤ c)
    (hab : b ≤ a) (hab' : -b ≤ a) (hcb : b ≤ c) (hcb' : -b ≤ c) :
    ∀ v : Fin 2 → K, 0 ≤ v ⬝ᵥ (!![a, b; b, c] : Matrix (Fin 2) (Fin 2) K) *ᵥ v := by
  intro v
  have h := wdd_form_nonneg a b c (v 0) (v 1) ha hc hab hab' hcb hcb'
  simp [dotProduct, mulVec, Fin.sum_univ_two]
  nlinarith [h]

#print axioms wdd_form_nonneg
#print axioms isPSD_of_wdd_2x2
