import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Notation

/-
Concrete lifting identity for the scalar single-keep-out QCQP (the base S-lemma case, the
primitive shipped in `podium.verify.bracket`).  This removes the abstract `lifting`
hypothesis of `sprocedure_matrix_sound` for this case: the bordered 2×2 matrix M(λ,t)'s
quadratic form at the lift (x,1) IS the Lagrangian minorant f0 x − λ f1 x − t.  Combined
with IsPSD (what the exact LDLᵀ checker decides), this is end-to-end soundness for the
scalar single keep-out.  No `sorry`.
-/

open Matrix
variable {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- Bordered S-procedure matrix for `f0 x = p0 x² + q0 x + r0`, `f1 x = p1 x² + q1 x + r1`. -/
def M2 (p0 q0 r0 p1 q1 r1 lam t : K) : Matrix (Fin 2) (Fin 2) K :=
  !![p0 - lam * p1, (q0 - lam * q1) / 2;
     (q0 - lam * q1) / 2, r0 - lam * r1 - t]

/-- The lift `x ↦ (x, 1)`. -/
def lift2 (x : K) : Fin 2 → K := ![x, 1]

/-- **Lifting identity.** `(x,1)ᵀ M(λ,t) (x,1) = f0 x − λ f1 x − t`. -/
theorem lifting2 (p0 q0 r0 p1 q1 r1 lam t x : K) :
    (lift2 x) ⬝ᵥ (M2 p0 q0 r0 p1 q1 r1 lam t) *ᵥ (lift2 x)
      = (p0 * x ^ 2 + q0 * x + r0) - lam * (p1 * x ^ 2 + q1 * x + r1) - t := by
  simp [M2, lift2, dotProduct, mulVec, Fin.sum_univ_two]
  ring

/-- PSD predicate (what the exact band-LDLᵀ checker decides). -/
def IsPSD2 (M : Matrix (Fin 2) (Fin 2) K) : Prop := ∀ v : Fin 2 → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- **End-to-end soundness, scalar single keep-out.** If the bordered matrix is PSD, then
`t` lower-bounds `f0` on the feasible set `{x : f1 x ≥ 0}` — no lifting hypothesis. -/
theorem scalar_sprocedure_sound
    (p0 q0 r0 p1 q1 r1 lam t : K) (hlam : 0 ≤ lam)
    (hPSD : IsPSD2 (M2 p0 q0 r0 p1 q1 r1 lam t)) :
    ∀ x, 0 ≤ p1 * x ^ 2 + q1 * x + r1 → t ≤ p0 * x ^ 2 + q0 * x + r0 := by
  intro x hx
  have h := hPSD (lift2 x)
  rw [lifting2] at h
  have hmul : 0 ≤ lam * (p1 * x ^ 2 + q1 * x + r1) := mul_nonneg hlam hx
  linarith

/-- **Capstone (scalar single keep-out): PSD certificate ⇒ certified global optimum.**
Composing the lifting, soundness, and closure: if the bordered matrix is PSD (`IsPSD2`, what
the exact `LDLᵀ` checker decides), `λ ≥ 0`, and a feasible witness `xbar` attains
`f0 xbar = t`, then `xbar` is a global minimiser of `f0` on the feasible set. Fully
mechanized for the shipped `podium.verify.bracket` primitive, no `sorry`. -/
theorem scalar_certified_optimum
    (p0 q0 r0 p1 q1 r1 lam t : K) (hlam : 0 ≤ lam)
    (hPSD : IsPSD2 (M2 p0 q0 r0 p1 q1 r1 lam t))
    (xbar : K) (hfeas : 0 ≤ p1 * xbar ^ 2 + q1 * xbar + r1)
    (hclose : p0 * xbar ^ 2 + q0 * xbar + r0 = t) :
    ∀ x, 0 ≤ p1 * x ^ 2 + q1 * x + r1 →
      p0 * xbar ^ 2 + q0 * xbar + r0 ≤ p0 * x ^ 2 + q0 * x + r0 := by
  intro x hx
  rw [hclose]
  exact scalar_sprocedure_sound p0 q0 r0 p1 q1 r1 lam t hlam hPSD x hx

#print axioms lifting2
#print axioms scalar_sprocedure_sound
#print axioms scalar_certified_optimum
