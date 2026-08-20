import Mathlib.Tactic.Linarith
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-
Closure / attainment leg (paper Thm 1, second leg): the S-procedure lower bound certifies
`t ≤ J*`; this file certifies the MATCHING UPPER bound — that a witness `x̄` which minimizes
the (convex) Lagrangian and satisfies complementary slackness is a *global* optimum, so the
bracket closes to equality `f0 x̄ = J*`.  Abstract (no convexity machinery needed): the only
facts used are that `x̄` minimizes the Lagrangian `L = f0 − Σλᵢ fᵢ` and that `Σλᵢ fᵢ x̄ = 0`
(aggregate complementary slackness).  Over any `LinearOrderedField`.  No `sorry`.
-/

variable {X : Type*} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- **Certified global optimum (attainment leg).** If `λ ≥ 0`, `x̄` minimizes the Lagrangian
`L(x) = f0 x − Σλᵢ fᵢ x` (which for a *convex* `L` is exactly stationarity `∇L(x̄)=0`), and
complementary slackness holds (`Σλᵢ fᵢ x̄ = 0`), then `x̄` is a global optimum: `f0 x̄ ≤ f0 x`
for every feasible `x`.  Combined with the S-procedure lower bound `t ≤ f0` this closes the
bracket — `f0 x̄` is the exact optimal value. -/
theorem certified_global_optimum (f0 : X → K) {m : ℕ} (f : Fin m → X → K) (lam : Fin m → K)
    (xbar : X) (hlam : ∀ i, 0 ≤ lam i)
    (Lmin : ∀ x, f0 xbar - ∑ i, lam i * f i xbar ≤ f0 x - ∑ i, lam i * f i x)
    (cs : ∑ i, lam i * f i xbar = 0) :
    ∀ x, (∀ i, 0 ≤ f i x) → f0 xbar ≤ f0 x := by
  intro x hx
  have hsum : 0 ≤ ∑ i, lam i * f i x :=
    Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hx i))
  have hL := Lmin x
  rw [cs] at hL
  linarith

/-- **Bracket closes to equality.** With the attainment leg above, `f0 x̄` is simultaneously a
lower bound (via the certified `t ≤ f0` on feasibles, applied at `x̄`) and the global minimum,
so any certified lower bound `t` with `t ≤ f0 x̄` and the optimality of `x̄` pins the optimal
value into the point interval `[t, f0 x̄]` collapsing at `x̄`. -/
theorem bracket_exact (f0 : X → K) {m : ℕ} (f : Fin m → X → K) (lam : Fin m → K)
    (xbar : X) (t : K) (hlam : ∀ i, 0 ≤ lam i)
    (Lmin : ∀ x, f0 xbar - ∑ i, lam i * f i xbar ≤ f0 x - ∑ i, lam i * f i x)
    (cs : ∑ i, lam i * f i xbar = 0)
    (lower : ∀ x, (∀ i, 0 ≤ f i x) → t ≤ f0 x)
    (hfeas : ∀ i, 0 ≤ f i xbar) :
    t ≤ f0 xbar ∧ ∀ x, (∀ i, 0 ≤ f i x) → f0 xbar ≤ f0 x :=
  ⟨lower xbar hfeas, certified_global_optimum f0 f lam xbar hlam Lmin cs⟩

#print axioms certified_global_optimum
#print axioms bracket_exact
