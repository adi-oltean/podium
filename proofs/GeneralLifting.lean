import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Data.Matrix.Mul
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-
General-n lifting identity: the bordered (n+1)×(n+1) matrix built from Hessian `P`, linear
half-vector `b` (= q/2), and corner `c` has quadratic form at the lift `snoc x 1` equal to
`xᵀ P x + 2 bᵀ x + c` (= f0 x − Σλ fₖ x − t with `P = P0−Σλ Pₖ`, `q = q0−Σλ qₖ`, `c = r0−Σλ rₖ−t`).
This removes the scalar restriction of LiftingIdentity.lean.  No `sorry`.
-/

open Matrix Finset

variable {K : Type*} [CommRing K] {n : ℕ}

/-- Bordered matrix: `P` in the top-left block, `b` on the last row/col, `c` in the corner. -/
def bordered (P : Matrix (Fin n) (Fin n) K) (b : Fin n → K) (c : K) :
    Matrix (Fin (n + 1)) (Fin (n + 1)) K :=
  Matrix.of fun i j =>
    Fin.lastCases (Fin.lastCases c (fun j' => b j') j)
      (fun i' => Fin.lastCases (b i') (fun j' => P i' j') j) i

theorem bordered_lifting (P : Matrix (Fin n) (Fin n) K) (b : Fin n → K) (c : K)
    (x : Fin n → K) :
    (Fin.snoc x 1) ⬝ᵥ (bordered P b c) *ᵥ (Fin.snoc x 1)
      = x ⬝ᵥ P *ᵥ x + 2 * (b ⬝ᵥ x) + c := by
  simp only [dotProduct, mulVec, bordered, Matrix.of_apply, Fin.sum_univ_castSucc,
    Fin.snoc_castSucc, Fin.snoc_last, Fin.lastCases_last, Fin.lastCases_castSucc,
    mul_one, one_mul, mul_add, Finset.sum_add_distrib]
  have h : (∑ i : Fin n, x i * b i) = (∑ i : Fin n, b i * x i) := by
    apply Finset.sum_congr rfl; intro i _; ring
  rw [h]; ring

section Soundness
variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K] {n : ℕ}

/-- PSD predicate (what the exact band-`LDLᵀ` checker decides). -/
def IsPSD {ι : Type*} [Fintype ι] (M : Matrix ι ι K) : Prop := ∀ v : ι → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- IsPSD of the bordered matrix ⇒ the quadratic form `xᵀPx + 2bᵀx + c ≥ 0` for all `x`. -/
theorem bordered_form_nonneg (P : Matrix (Fin n) (Fin n) K) (b : Fin n → K) (c : K)
    (h : IsPSD (bordered P b c)) (x : Fin n → K) :
    0 ≤ x ⬝ᵥ P *ᵥ x + 2 * (b ⬝ᵥ x) + c := by
  have := h (Fin.snoc x 1)
  rwa [bordered_lifting] at this

/-- **General-`n` end-to-end soundness.** IsPSD of the bordered certificate matrix, `λ ≥ 0`,
and the QCQP form-identity (`xᵀPx + 2bᵀx + c = f0 x − Σλᵢ fᵢ x − t`) ⇒ `t ≤ f0` on the feasible
set — the S-procedure lower bound for the general block-banded QCQP. -/
theorem general_sprocedure_sound {X : Type*}
    (P : Matrix (Fin n) (Fin n) K) (b : Fin n → K) (c : K) (h : IsPSD (bordered P b c))
    (lift : X → (Fin n → K)) (f0 : X → K) {m : ℕ} (f : Fin m → X → K) (lam : Fin m → K) (t : K)
    (hlam : ∀ i, 0 ≤ lam i)
    (hform : ∀ y, (lift y) ⬝ᵥ P *ᵥ (lift y) + 2 * (b ⬝ᵥ (lift y)) + c
                    = f0 y - ∑ i, lam i * f i y - t) :
    ∀ y, (∀ i, 0 ≤ f i y) → t ≤ f0 y := by
  intro y hy
  have hnn := bordered_form_nonneg P b c h (lift y)
  rw [hform y] at hnn
  have hsum : 0 ≤ ∑ i, lam i * f i y :=
    Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hy i))
  linarith
end Soundness

#print axioms bordered_lifting
#print axioms general_sprocedure_sound
