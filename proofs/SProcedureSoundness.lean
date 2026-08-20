import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Matrix.Mul

/-
Mechanized SOUNDNESS of the exact-rational S-procedure lower-bound certificate
(paper Thm 1 lower leg / Prop generic; note thm:sound and thm:multi).

The trusted checker accepts a multiplier `lam ≥ 0` and scalar `t` when the lifted matrix
`M(lam,t) ⪰ 0`.  That PSD condition is exactly the statement that the Lagrangian minorant
`f0 - Σ lam_i f_i - t` is nonnegative for all x (its value is `[x;1]ᵀ M [x;1] ≥ 0`).  Given
that, soundness — `t ≤ f0 x` on the feasible set, hence `t ≤ J*` — is weak duality, proved
here with no floating point and no `sorry`.  This is the guarantee the certificate rests on.

Stated over an arbitrary `LinearOrderedField K`, so it covers both the exact rational
checker (`K = ℚ`) and the real optimum (`K = ℝ`).
-/

variable {X : Type*} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- Single keep-out (S-lemma).  If `lam ≥ 0` and the Lagrangian minorant dominates `t`
everywhere, then `t` lower-bounds the objective on the feasible set `{x : f1 x ≥ 0}`. -/
theorem sprocedure_lower_bound
    (f0 f1 : X → K) (lam t : K) (hlam : 0 ≤ lam)
    (minorant : ∀ x, t ≤ f0 x - lam * f1 x) :
    ∀ x, 0 ≤ f1 x → t ≤ f0 x := by
  intro x hx
  have hmul : 0 ≤ lam * f1 x := mul_nonneg hlam hx
  have := minorant x
  linarith

/-- Multi keep-out (chain / Celis–Dennis–Tapia, note thm:multi).  If every `lam i ≥ 0` and
the aggregated Lagrangian minorant dominates `t` everywhere, then `t` lower-bounds the
objective on the feasible set `{x : ∀ i, f i x ≥ 0}` — for any number of constraints `m`. -/
theorem sprocedure_lower_bound_multi
    (f0 : X → K) {m : ℕ} (f : Fin m → X → K) (lam : Fin m → K) (t : K)
    (hlam : ∀ i, 0 ≤ lam i)
    (minorant : ∀ x, t ≤ f0 x - ∑ i, lam i * f i x) :
    ∀ x, (∀ i, 0 ≤ f i x) → t ≤ f0 x := by
  intro x hx
  have hsum : 0 ≤ ∑ i, lam i * f i x :=
    Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hx i))
  have := minorant x
  linarith

/-- The bracket is sound: for any feasible witness `x̄`, `[t, f0 x̄]` brackets the optimum
(`t ≤ f0 x̄`).  Closure is the case `t = f0 x̄`. -/
theorem sprocedure_bracket_sound
    (f0 f1 : X → K) (lam t : K) (hlam : 0 ≤ lam)
    (minorant : ∀ x, t ≤ f0 x - lam * f1 x)
    (xbar : X) (hfeas : 0 ≤ f1 xbar) :
    t ≤ f0 xbar :=
  sprocedure_lower_bound f0 f1 lam t hlam minorant xbar hfeas

/- ## The matrix bridge: from the checker's PSD verdict to soundness -/

open Matrix

/-- Our PSD predicate: the quadratic form of `M` is nonnegative everywhere. This is exactly
what the exact-rational band-`LDLᵀ` checker establishes (all pivots ≥ 0). -/
def IsPSD {ι : Type*} [Fintype ι] (M : Matrix ι ι K) : Prop :=
  ∀ v : ι → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- **Matrix bridge.** If the certificate matrix `M` is PSD and its quadratic form at the
lift `[x;1]` equals the Lagrangian minorant `f0 x − Σ λᵢ fᵢ x − t`, then `t` lower-bounds the
objective on the feasible set. This connects the checker's PSD verdict `IsPSD M` to
soundness; the only remaining (QCQP-specific, purely algebraic) hypothesis is the lifting
identity `lifting`. -/
theorem sprocedure_matrix_sound {ι : Type*} [Fintype ι]
    (M : Matrix ι ι K) (hM : IsPSD M)
    (f0 : X → K) {m : ℕ} (f : Fin m → X → K) (lam : Fin m → K) (t : K)
    (hlam : ∀ i, 0 ≤ lam i) (lift : X → (ι → K))
    (lifting : ∀ x, (lift x) ⬝ᵥ M *ᵥ (lift x) = f0 x - ∑ i, lam i * f i x - t) :
    ∀ x, (∀ i, 0 ≤ f i x) → t ≤ f0 x := by
  refine sprocedure_lower_bound_multi f0 f lam t hlam ?_
  intro x
  have h := hM (lift x)
  rw [lifting x] at h
  linarith

/-- **Bracket closure ⇒ certified global optimum.** If `t` lower-bounds `f0` on the feasible
set (e.g. via the S-procedure) and a feasible witness `xbar` attains `f0 xbar = t`, then
`xbar` is a global minimiser over the feasible set. This is the closure step: a matching
lower bound and feasible upper bound certify the optimum. -/
theorem bracket_closes (f0 : X → K) (feas : X → Prop) (t : K)
    (lower : ∀ x, feas x → t ≤ f0 x) (xbar : X) (hfeas : feas xbar) (hclose : f0 xbar = t) :
    ∀ x, feas x → f0 xbar ≤ f0 x := by
  intro x hx
  rw [hclose]
  exact lower x hx

#print axioms sprocedure_lower_bound
#print axioms sprocedure_lower_bound_multi
#print axioms sprocedure_matrix_sound
#print axioms bracket_closes
