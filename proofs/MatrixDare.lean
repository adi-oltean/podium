import Mathlib.Algebra.Group.Basic

/-
Matrix DARE horizon-independence (paper Lemma "dare", matrix version): seeding the
block-Schur/Riccati recursion at a stationary solution `S*` with the DARE-generating
running cost `diag = S* + GᵀS*⁻¹G` makes every pivot equal `S*` (horizon-independent
verification).  The pivot is constant because the correction term `GᵀS⁻¹G` is *added* in
the running-cost block and *subtracted* in the Schur step — pure additive cancellation.
So the result holds at the `AddCommGroup` level, with the correction as an arbitrary map
`corr : M → M` (in the real setting `M = d×d` matrices and `corr S = GᵀS⁻¹G`).  No `sorry`.
-/

variable {M : Type*} [AddCommGroup M]

/-- One Schur/Riccati step with running-cost block `diag` and correction map `corr`:
`S ↦ diag − corr S`  (in the matrix setting `corr S = GᵀS⁻¹G`). -/
def dareStep (corr : M → M) (diag S : M) : M := diag - corr S

/-- The seeded orbit of the recursion. -/
def dareOrbit (corr : M → M) (diag S0 : M) : ℕ → M
  | 0 => S0
  | n + 1 => dareStep corr diag (dareOrbit corr diag S0 n)

/-- `S*` is a fixed point when the running cost is the DARE-generating datum
`diag = S* + corr S*`. -/
theorem dareStep_fixed (corr : M → M) (Sstar : M) :
    dareStep corr (Sstar + corr Sstar) Sstar = Sstar := by
  rw [dareStep, add_sub_cancel_right]

/-- **Horizon-independence.** Seeded at `S*` with the DARE running cost, every pivot is `S*`
(the certificate is horizon-independent). -/
theorem dareOrbit_constant (corr : M → M) (Sstar : M) :
    ∀ n, dareOrbit corr (Sstar + corr Sstar) Sstar n = Sstar := by
  intro n
  induction n with
  | zero => rfl
  | succ k ih =>
      rw [dareOrbit, ih, dareStep_fixed]

#print axioms dareOrbit_constant
