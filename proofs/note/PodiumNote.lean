import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Notation

/-
Lean treatment of the previous Podium note
  "Exact-Rational Certificates in Podium" (rpod/docs/exact-arithmetic-certificates/note.tex).
Mechanizes the soundness core; the recovery-RATE theorems (thm:recovery O(1/D²),
thm:hard linear) are analytic (Taylor / one-sided derivative) and are noted, not mechanized.
Everything here is over a general `LinearOrderedField` (covers ℚ and ℝ), no `sorry`.
-/

open Matrix
variable {X : Type*} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- **Note Theorem 1 (Lower bound / thm:sound), single constraint.** `λ ≥ 0` and the
minorant `f0 − λ f1 − t ≥ 0` (which `M(λ,t) ⪰ 0` provides via the lifting) ⇒ `t ≤ f0` on
the feasible set `{x : f1 x ≥ 0}`, hence `t ≤ J*`. -/
theorem note_thm_sound (f0 f1 : X → K) (lam t : K) (hlam : 0 ≤ lam)
    (minorant : ∀ x, t ≤ f0 x - lam * f1 x) :
    ∀ x, 0 ≤ f1 x → t ≤ f0 x := by
  intro x hx
  have := minorant x
  nlinarith [mul_nonneg hlam hx]

/-- **Note Theorem 4 (Certified duality gap / thm:multi), general m.** Same, for any number
of constraints — the minorant argument does not use `m = 1`. -/
theorem note_thm_multi (f0 : X → K) {m : ℕ} (f : Fin m → X → K) (lam : Fin m → K) (t : K)
    (hlam : ∀ i, 0 ≤ lam i)
    (minorant : ∀ x, t ≤ f0 x - ∑ i, lam i * f i x) :
    ∀ x, (∀ i, 0 ≤ f i x) → t ≤ f0 x := by
  intro x hx
  have hsum : 0 ≤ ∑ i, lam i * f i x :=
    Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hx i))
  have := minorant x
  linarith

/-- **Note bracket (thm:multi closing observation).** A certified lower bound `t ≤ f0 x`
(all feasible `x`) and an exactly feasible witness `xbar` place `J*` in `[t, f0 xbar]`; the
gap `f0 xbar − t` bounds the suboptimality of `xbar`. Here: `t ≤ f0 xbar` and, for any feasible
`x`, `f0 xbar − (f0 xbar − t) = t ≤ f0 x`. -/
theorem note_bracket (f0 : X → K) (feas : X → Prop) (t : K)
    (lower : ∀ x, feas x → t ≤ f0 x) (xbar : X) (hbar : feas xbar) :
    t ≤ f0 xbar ∧ ∀ x, feas x → f0 xbar - (f0 xbar - t) ≤ f0 x := by
  refine ⟨lower xbar hbar, ?_⟩
  intro x hx; simpa using lower x hx

/-- **Note Lemma (Schur complement / lem:schur), scalar case `A = a > 0`.** For the bordered
matrix `!![a, b/2; b/2, c−t]`, PSD ⟺ `t ≤ g` where `g = c − b²/(4a)` is the dual value
`inf_x (a x² + b x + c)`. -/
theorem note_lem_schur (a b c t : K) (ha : 0 < a) :
    (∀ v : Fin 2 → K, 0 ≤ v ⬝ᵥ (!![a, b/2; b/2, c - t] : Matrix (Fin 2) (Fin 2) K) *ᵥ v)
      ↔ t ≤ c - b ^ 2 / (4 * a) := by
  constructor
  · intro h
    have hv := h ![b, -2 * a]
    simp [dotProduct, mulVec, Fin.sum_univ_two] at hv
    have key : b ^ 2 ≤ (c - t) * (4 * a) := by nlinarith [hv, ha]
    have hdiv : b ^ 2 / (4 * a) ≤ c - t := (div_le_iff₀ (by positivity)).2 key
    linarith
  · intro h v
    simp [dotProduct, mulVec, Fin.sum_univ_two]
    have hg' : b ^ 2 ≤ (c - t) * (4 * a) := by
      rw [← div_le_iff₀ (by positivity)]; linarith
    nlinarith [mul_nonneg ha.le (sq_nonneg (2 * a * v 0 + b * v 1)), ha, hg', sq_nonneg (v 1)]

#print axioms note_thm_sound
#print axioms note_thm_multi
#print axioms note_lem_schur
