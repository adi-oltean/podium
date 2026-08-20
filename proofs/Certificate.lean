import Mathlib.Data.Matrix.Mul
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-
HEADLINE capstone — the exact-certificate pipeline as ONE self-contained theorem.
Chains the pieces proved separately elsewhere (LDLTSound, GeneralLifting, ClosureOptimum):

  exact checker ACCEPTS  (M = LᵀDL with pivots d ≥ 0)
    ⇒ IsPSD M                                              [accept ⇒ PSD]
    ⇒ t ≤ f0 on the feasible set                           [S-procedure weak duality]
  and if the witness x̄ minimizes the Lagrangian with complementary slackness
    ⇒ f0 x̄ is the global minimum, with t ≤ f0 x̄            [attainment; bracket closes]

So a single accepted exact-ℚ certificate yields a machine-checked *global optimum with a
matching lower bound* — `f0 x̄ = J*`.  Self-contained (re-proves the two small steps inline),
over any `LinearOrderedField`, no `sorry`.
-/

open Matrix
variable {n : ℕ} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- PSD predicate the exact checker decides. -/
def IsPSD (M : Matrix (Fin n) (Fin n) K) : Prop := ∀ v : Fin n → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- Checker soundness (accept ⇒ PSD): `M = LᵀDL` with pivots `d ≥ 0` ⇒ `M` PSD. -/
theorem accept_isPSD (L : Matrix (Fin n) (Fin n) K) (d : Fin n → K) (hd : ∀ i, 0 ≤ d i)
    {M : Matrix (Fin n) (Fin n) K} (hM : M = Lᵀ * diagonal d * L) : IsPSD M := by
  intro v
  have key : v ⬝ᵥ M *ᵥ v = (L *ᵥ v) ⬝ᵥ (diagonal d) *ᵥ (L *ᵥ v) := by
    rw [hM, Matrix.mul_assoc, ← Matrix.mulVec_mulVec, dotProduct_mulVec, vecMul_transpose,
      ← Matrix.mulVec_mulVec]
  rw [key]
  have hform : (L *ᵥ v) ⬝ᵥ (diagonal d) *ᵥ (L *ᵥ v) = ∑ i, d i * ((L *ᵥ v) i * (L *ᵥ v) i) := by
    simp only [dotProduct, mulVec_diagonal]
    exact Finset.sum_congr rfl (fun i _ => mul_left_comm _ _ _)
  rw [hform]
  exact Finset.sum_nonneg (fun i _ => mul_nonneg (hd i) (mul_self_nonneg _))

/-- **The exact-certificate capstone.** From a single accepted exact certificate plus a KKT
witness: `t` is a certified lower bound on the feasible set, `x̄` attains the global minimum,
and `t ≤ f0 x̄` — so `f0 x̄ = J*` is machine-checked. -/
theorem certified_optimum_pipeline {X : Type*}
    -- exact checker accepts the S-procedure matrix M:
    (L : Matrix (Fin n) (Fin n) K) (d : Fin n → K) (hd : ∀ i, 0 ≤ d i)
    (M : Matrix (Fin n) (Fin n) K) (hM : M = Lᵀ * diagonal d * L)
    -- S-procedure data: lift, objective f0, constraints f, multipliers λ ≥ 0, level t:
    (lift : X → (Fin n → K)) (f0 : X → K) {m : ℕ} (f : Fin m → X → K)
    (lam : Fin m → K) (t : K) (hlam : ∀ i, 0 ≤ lam i)
    -- the lifting / S-procedure form identity (M's quadratic form = Lagrangian − t):
    (hform : ∀ y, (lift y) ⬝ᵥ M *ᵥ (lift y) = f0 y - ∑ i, lam i * f i y - t)
    -- KKT witness x̄: feasible, minimizes the Lagrangian, complementary slackness:
    (xbar : X) (hfeas : ∀ i, 0 ≤ f i xbar)
    (Lmin : ∀ y, f0 xbar - ∑ i, lam i * f i xbar ≤ f0 y - ∑ i, lam i * f i y)
    (cs : ∑ i, lam i * f i xbar = 0) :
    (∀ y, (∀ i, 0 ≤ f i y) → t ≤ f0 y)              -- certified lower bound (t ≤ J*)
      ∧ (∀ y, (∀ i, 0 ≤ f i y) → f0 xbar ≤ f0 y)    -- x̄ is a global optimum
      ∧ t ≤ f0 xbar := by
  have hpsd : IsPSD M := accept_isPSD L d hd hM
  -- (1) certified lower bound via weak duality
  have lower : ∀ y, (∀ i, 0 ≤ f i y) → t ≤ f0 y := by
    intro y hy
    have hnn := hpsd (lift y)
    rw [hform y] at hnn
    have hsum : 0 ≤ ∑ i, lam i * f i y :=
      Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hy i))
    linarith
  -- (2) attainment: x̄ is a global optimum
  have attain : ∀ y, (∀ i, 0 ≤ f i y) → f0 xbar ≤ f0 y := by
    intro y hy
    have hsum : 0 ≤ ∑ i, lam i * f i y :=
      Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hy i))
    have hL := Lmin y
    rw [cs] at hL
    linarith
  exact ⟨lower, attain, lower xbar hfeas⟩

/-- **Verdict-level capstone.** Same conclusion as `certified_optimum_pipeline` but starting from
the checker's *verdict* `IsPSD M` directly, rather than a specific `LᵀDL` factorization. This is
the form that composes with the decision criterion `isPSD_iff_ldlt` (either accept route feeds it):
`IsPSD M` + the S-procedure form identity + a KKT witness ⇒ certified global optimum with matching
bracket, `f0 x̄ = J*`. -/
theorem certified_optimum_of_verdict {X : Type*}
    (M : Matrix (Fin n) (Fin n) K) (hpsd : IsPSD M)
    (lift : X → (Fin n → K)) (f0 : X → K) {m : ℕ} (f : Fin m → X → K)
    (lam : Fin m → K) (t : K) (hlam : ∀ i, 0 ≤ lam i)
    (hform : ∀ y, (lift y) ⬝ᵥ M *ᵥ (lift y) = f0 y - ∑ i, lam i * f i y - t)
    (xbar : X) (hfeas : ∀ i, 0 ≤ f i xbar)
    (Lmin : ∀ y, f0 xbar - ∑ i, lam i * f i xbar ≤ f0 y - ∑ i, lam i * f i y)
    (cs : ∑ i, lam i * f i xbar = 0) :
    (∀ y, (∀ i, 0 ≤ f i y) → t ≤ f0 y)
      ∧ (∀ y, (∀ i, 0 ≤ f i y) → f0 xbar ≤ f0 y)
      ∧ t ≤ f0 xbar := by
  have lower : ∀ y, (∀ i, 0 ≤ f i y) → t ≤ f0 y := by
    intro y hy
    have hnn := hpsd (lift y)
    rw [hform y] at hnn
    have hsum : 0 ≤ ∑ i, lam i * f i y :=
      Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hy i))
    linarith
  have attain : ∀ y, (∀ i, 0 ≤ f i y) → f0 xbar ≤ f0 y := by
    intro y hy
    have hsum : 0 ≤ ∑ i, lam i * f i y :=
      Finset.sum_nonneg (fun i _ => mul_nonneg (hlam i) (hy i))
    have hL := Lmin y
    rw [cs] at hL
    linarith
  exact ⟨lower, attain, lower xbar hfeas⟩

/-! ### A concrete machine-checked instance (a 1-D keep-out QCQP over ℚ)

Minimize `f0(x) = (x − ½)²` subject to the keep-out `f1(x) = x² − 1 ≥ 0` (i.e. `|x| ≥ 1`).
The unconstrained minimizer `x = ½` is infeasible; the exact certificate — multiplier `λ = ½`,
level `t = ¼`, bordered PSD matrix `M = [[½, −½],[−½, ½]]`, witness `x̄ = 1` — is machine-checked
to prove `x̄ = 1` is a *global* optimum with value `J* = ¼`.  Everything is exact rational. -/
section ConcreteInstance

/-- The bordered S-procedure matrix `M(λ,t) = [[½,−½],[−½,½]]` is PSD: `vᵀMv = ½(v₀−v₁)² ≥ 0`. -/
theorem keepout_M_isPSD :
    IsPSD (n := 2) (K := ℚ) !![1/2, -1/2; -1/2, 1/2] := by
  intro v
  have hv : v ⬝ᵥ !![1/2, -1/2; -1/2, 1/2] *ᵥ v = 1 / 2 * (v 0 - v 1) ^ 2 := by
    simp only [dotProduct, mulVec, Fin.sum_univ_two, Matrix.of_apply, Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.head_cons]
    ring
  rw [hv]
  have := sq_nonneg (v 0 - v 1)
  linarith

/-- **Certified global optimum of the keep-out QCQP.** `x̄ = 1` minimizes `(x−½)²` over `|x| ≥ 1`,
with certified value `¼`: for every feasible `y` (`y²−1 ≥ 0`), `¼ ≤ (y−½)²` and `(1−½)² ≤ (y−½)²`. -/
theorem keepout_certified :
    (∀ y : ℚ, 0 ≤ y ^ 2 - 1 → (1 / 4 : ℚ) ≤ (y - 1 / 2) ^ 2)
      ∧ (∀ y : ℚ, 0 ≤ y ^ 2 - 1 → (1 - 1 / 2 : ℚ) ^ 2 ≤ (y - 1 / 2) ^ 2)
      ∧ (1 / 4 : ℚ) ≤ (1 - 1 / 2) ^ 2 := by
  obtain ⟨lower, attain, hbar⟩ := certified_optimum_of_verdict (n := 2) (K := ℚ)
    (M := !![1/2, -1/2; -1/2, 1/2]) keepout_M_isPSD
    (lift := fun x => ![x, 1]) (f0 := fun x => (x - 1 / 2) ^ 2)
    (m := 1) (f := fun _ x => x ^ 2 - 1) (lam := fun _ => 1 / 2) (t := 1 / 4)
    (hlam := by intro i; norm_num)
    (hform := by
      intro y
      simp only [dotProduct, mulVec, Fin.sum_univ_two, Fin.sum_univ_one, Matrix.of_apply,
        Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
      ring)
    (xbar := 1) (hfeas := by intro i; norm_num)
    (Lmin := by intro y; simp only [Fin.sum_univ_one]; nlinarith [sq_nonneg (y - 1)])
    (cs := by simp only [Fin.sum_univ_one]; norm_num)
  exact ⟨fun y hy => lower y (fun _ => hy), fun y hy => attain y (fun _ => hy), hbar⟩

end ConcreteInstance

/-! ### A 2-D (vector-state) instance — a disk keep-out (RPOD-flavored)

Minimize `f0(x,y) = (x−½)² + y²` (target `(½,0)`) subject to the disk keep-out
`f1(x,y) = x²+y²−1 ≥ 0` (stay outside the unit disk).  The target is inside the disk (infeasible);
the exact certificate (`λ=½`, `t=¼`, bordered `3×3` PSD matrix, witness `(1,0)`) proves `(1,0)`
globally optimal with `J*=¼`.  This exercises the whole pipeline on a genuine vector state. -/
section DiskInstance

/-- The bordered `3×3` matrix for the disk instance is PSD: `vᵀMv = ½(v₀−v₂)² + ½v₁² ≥ 0`. -/
theorem disk_M_isPSD :
    IsPSD (n := 3) (K := ℚ) !![1/2, 0, -1/2; 0, 1/2, 0; -1/2, 0, 1/2] := by
  intro v
  have hv : v ⬝ᵥ !![1/2, 0, -1/2; 0, 1/2, 0; -1/2, 0, 1/2] *ᵥ v
      = 1 / 2 * (v 0 - v 2) ^ 2 + 1 / 2 * (v 1) ^ 2 := by
    simp only [dotProduct, mulVec, Fin.sum_univ_three, Matrix.of_apply, Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two, Matrix.tail_cons]
    ring
  rw [hv]
  nlinarith [sq_nonneg (v 0 - v 2), sq_nonneg (v 1)]

/-- **Certified global optimum of the disk keep-out.** `(1,0)` minimizes `(x−½)²+y²` over
`x²+y² ≥ 1` with certified value `¼`. -/
theorem disk_certified :
    (∀ p : Fin 2 → ℚ, 0 ≤ (p 0) ^ 2 + (p 1) ^ 2 - 1 →
        (1 / 4 : ℚ) ≤ (p 0 - 1 / 2) ^ 2 + (p 1) ^ 2)
      ∧ (∀ p : Fin 2 → ℚ, 0 ≤ (p 0) ^ 2 + (p 1) ^ 2 - 1 →
        ((1 : ℚ) - 1 / 2) ^ 2 + (0 : ℚ) ^ 2 ≤ (p 0 - 1 / 2) ^ 2 + (p 1) ^ 2) := by
  obtain ⟨lower, attain, _⟩ := certified_optimum_of_verdict (n := 3) (K := ℚ)
    (M := !![1/2, 0, -1/2; 0, 1/2, 0; -1/2, 0, 1/2]) disk_M_isPSD
    (lift := fun p => ![p 0, p 1, 1]) (f0 := fun p => (p 0 - 1 / 2) ^ 2 + (p 1) ^ 2)
    (m := 1) (f := fun _ p => (p 0) ^ 2 + (p 1) ^ 2 - 1) (lam := fun _ => 1 / 2) (t := 1 / 4)
    (hlam := by intro i; norm_num)
    (hform := by
      intro p
      simp only [dotProduct, mulVec, Fin.sum_univ_three, Fin.sum_univ_one, Matrix.of_apply,
        Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
        Matrix.tail_cons]
      ring)
    (xbar := ![1, 0]) (hfeas := by intro i; norm_num)
    (Lmin := by
      intro p
      simp only [Fin.sum_univ_one, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
      nlinarith [sq_nonneg (p 0 - 1), sq_nonneg (p 1)])
    (cs := by simp only [Fin.sum_univ_one, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.head_cons]; norm_num)
  exact ⟨fun p hp => lower p (fun _ => hp), fun p hp => attain p (fun _ => hp)⟩

end DiskInstance

#print axioms certified_optimum_pipeline
#print axioms certified_optimum_of_verdict
#print axioms keepout_certified
#print axioms disk_certified
