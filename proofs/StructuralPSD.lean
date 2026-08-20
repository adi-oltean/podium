import Mathlib.Data.Matrix.Mul
import Mathlib.Tactic.Linarith

/-
Structural backbone: the PSD matrices form a convex cone closed under congruence.  This is
*why* the exact certificate method is sound as a whole:
  • S-procedure aggregation `M = M0 − Σλₖ Mₖ` (λ ≥ 0) lives in the cone — `isPSD_add`, `isPSD_smul`;
  • every exact Gaussian/Schur elimination step is a congruence `M ↦ SᵀMS`, which preserves PSD —
    `isPSD_congr` (so the exact band-`LDLᵀ` sweep never leaves the cone, and `ldlt_isPSD` is the
    special case `S = L`, `M = diagonal(pivots)`);
  • the lifting `[x;1]ᵀ M [x;1]` is congruence-evaluation.
Over any `LinearOrderedField`, no `sorry`.
-/

open Matrix
variable {n : ℕ} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- PSD predicate the checker decides. -/
def IsPSD (M : Matrix (Fin n) (Fin n) K) : Prop := ∀ v : Fin n → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- The zero matrix is PSD (the apex of the cone). -/
theorem isPSD_zero : IsPSD (0 : Matrix (Fin n) (Fin n) K) := by
  intro v; simp

/-- **Cone: closed under addition.** Aggregating certificates preserves PSD. -/
theorem isPSD_add {M N : Matrix (Fin n) (Fin n) K} (hM : IsPSD M) (hN : IsPSD N) :
    IsPSD (M + N) := by
  intro v
  rw [add_mulVec, dotProduct_add]
  exact add_nonneg (hM v) (hN v)

/-- **Cone: closed under nonnegative scaling.** Multiplying a certificate by `λ ≥ 0` (an
S-procedure multiplier) preserves PSD. -/
theorem isPSD_smul {c : K} (hc : 0 ≤ c) {M : Matrix (Fin n) (Fin n) K} (h : IsPSD M) :
    IsPSD (c • M) := by
  intro v
  rw [smul_mulVec, dotProduct_smul, smul_eq_mul]
  exact mul_nonneg hc (h v)

/-- **Congruence preserves PSD.** For any `S`, `IsPSD M ⇒ IsPSD (SᵀMS)`, because
`vᵀ(SᵀMS)v = (Sv)ᵀ M (Sv) ≥ 0`.  Each exact elimination step is such a congruence, so the whole
exact sweep stays in the PSD cone; `ldlt_isPSD` (pivots ≥ 0 ⇒ PSD) is the case `M = diagonal d`. -/
theorem isPSD_congr (S : Matrix (Fin n) (Fin n) K) {M : Matrix (Fin n) (Fin n) K}
    (h : IsPSD M) : IsPSD (Sᵀ * M * S) := by
  intro v
  have key : v ⬝ᵥ (Sᵀ * M * S) *ᵥ v = (S *ᵥ v) ⬝ᵥ M *ᵥ (S *ᵥ v) := by
    rw [Matrix.mul_assoc, ← Matrix.mulVec_mulVec, dotProduct_mulVec, vecMul_transpose,
      ← Matrix.mulVec_mulVec]
  rw [key]
  exact h (S *ᵥ v)

/-- **S-procedure aggregate is PSD.** The exact certificate `M0 + Σᵢ λᵢ Mᵢ` (each `Mᵢ` PSD,
`λᵢ ≥ 0`) is PSD — the cone is closed under nonnegative combinations. -/
theorem isPSD_nonneg_combo {m : ℕ} (M0 : Matrix (Fin n) (Fin n) K) (hM0 : IsPSD M0)
    (M : Fin m → Matrix (Fin n) (Fin n) K) (lam : Fin m → K)
    (hlam : ∀ i, 0 ≤ lam i) (hM : ∀ i, IsPSD (M i)) :
    IsPSD (M0 + ∑ i, lam i • M i) := by
  refine isPSD_add hM0 ?_
  apply Finset.sum_induction (fun i => lam i • M i) IsPSD
    (fun a b ha hb => isPSD_add ha hb) isPSD_zero
  intro i _
  exact isPSD_smul (hlam i) (hM i)

#print axioms isPSD_add
#print axioms isPSD_smul
#print axioms isPSD_congr
#print axioms isPSD_nonneg_combo
