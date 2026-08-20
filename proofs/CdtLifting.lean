import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Notation

/-
Concrete lifting + soundness for the scalar TWO-keep-out (Celis–Dennis–Tapia) case — the
adversarial regime.  The two-quadratic S-procedure is not lossless in general, but the
lower-bound SOUNDNESS holds for any λ₁,λ₂ ≥ 0: if the bordered aggregate matrix is PSD then
t lower-bounds f0 on the doubly-feasible set {x : f1 x ≥ 0 ∧ f2 x ≥ 0}.  No `sorry`.
-/

open Matrix
variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- Bordered aggregate matrix `P0 − λ₁P1 − λ₂P2` for scalar quadratics `fᵢ = pᵢ x² + qᵢ x + rᵢ`. -/
def Mcdt (p0 q0 r0 p1 q1 r1 p2 q2 r2 lam1 lam2 t : K) : Matrix (Fin 2) (Fin 2) K :=
  !![p0 - lam1 * p1 - lam2 * p2, (q0 - lam1 * q1 - lam2 * q2) / 2;
     (q0 - lam1 * q1 - lam2 * q2) / 2, r0 - lam1 * r1 - lam2 * r2 - t]

def lift2 (x : K) : Fin 2 → K := ![x, 1]

theorem cdt_lifting (p0 q0 r0 p1 q1 r1 p2 q2 r2 lam1 lam2 t x : K) :
    (lift2 x) ⬝ᵥ (Mcdt p0 q0 r0 p1 q1 r1 p2 q2 r2 lam1 lam2 t) *ᵥ (lift2 x)
      = (p0 * x ^ 2 + q0 * x + r0) - lam1 * (p1 * x ^ 2 + q1 * x + r1)
        - lam2 * (p2 * x ^ 2 + q2 * x + r2) - t := by
  simp [Mcdt, lift2, dotProduct, mulVec, Fin.sum_univ_two]
  ring

/-- **CDT soundness.** PSD of the aggregate bordered matrix + `λ₁,λ₂ ≥ 0` ⇒ `t` lower-bounds
`f0` on the doubly-feasible set — no lifting hypothesis. -/
theorem cdt_sprocedure_sound
    (p0 q0 r0 p1 q1 r1 p2 q2 r2 lam1 lam2 t : K) (h1 : 0 ≤ lam1) (h2 : 0 ≤ lam2)
    (hPSD : ∀ v : Fin 2 → K, 0 ≤ v ⬝ᵥ (Mcdt p0 q0 r0 p1 q1 r1 p2 q2 r2 lam1 lam2 t) *ᵥ v) :
    ∀ x, 0 ≤ p1 * x ^ 2 + q1 * x + r1 → 0 ≤ p2 * x ^ 2 + q2 * x + r2 →
      t ≤ p0 * x ^ 2 + q0 * x + r0 := by
  intro x hx1 hx2
  have h := hPSD (lift2 x)
  rw [cdt_lifting] at h
  have hm1 : 0 ≤ lam1 * (p1 * x ^ 2 + q1 * x + r1) := mul_nonneg h1 hx1
  have hm2 : 0 ≤ lam2 * (p2 * x ^ 2 + q2 * x + r2) := mul_nonneg h2 hx2
  linarith

#print axioms cdt_lifting
#print axioms cdt_sprocedure_sound
