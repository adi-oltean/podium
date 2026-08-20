import Mathlib.Tactic.Ring.Basic

/-
Feasibility probe: mechanize the scalar DARE-fixed-point identity (paper Lemma "dare",
docs/riccati-proofs.md Lemma C').  Pure field arithmetic.

Claim: for c = m + 1/m (Pell locus) the Riccati/Schur recursion d_{k+1} = c - 1/d_k,
seeded at d_0 = m, is constant: d_k = m for all k.  (Lean's 1/0 = 0 convention makes the
identity hold even at m = 0, so no side condition is needed.)
-/

-- The recursion over the rationals.
def step (c d : Rat) : Rat := c - 1 / d

-- m is a fixed point of the step for c = m + 1/m.
theorem dare_fixed_point (m : Rat) : step (m + 1 / m) m = m := by
  unfold step; exact add_sub_cancel_right m (1 / m)

-- The seeded orbit.
def orbit (c d0 : Rat) : Nat → Rat
  | 0 => d0
  | n + 1 => step c (orbit c d0 n)

-- Hence the whole seeded orbit is constant: every pivot equals m (O(1) in the horizon).
theorem dare_orbit_constant (m : Rat) : ∀ n, orbit (m + 1 / m) m n = m := by
  intro n
  induction n with
  | zero => rfl
  | succ k ih => unfold orbit; rw [ih]; exact dare_fixed_point m

#print axioms dare_orbit_constant
