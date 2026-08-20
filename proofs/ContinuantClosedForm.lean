import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-
A constructed continuant family with an exact closed form, giving explicit linear bit-growth.
The continuant recurrence `K_{k+2} = 5 K_{k+1} − 4 K_k` with `K_0 = 1`, `K_1 = 4` has closed
form `K_k = 4^k` (characteristic roots 4 and 1). Hence `bit(K_k) = 2k+1` grows linearly in the
horizon — a concrete instance of the paper's "constructed family with Θ(N) per-entry bit-size"
(the honest, family-specific version of the bit-size discussion; the general bound is open).
Proof: two-step induction. No `sorry`.
-/

/-- `K_{k+2}=5K_{k+1}−4K_k`, `K_0=1`, `K_1=4` ⇒ `K_k = 4^k` (so bit-size `= 2k+1`, linear). -/
theorem continuant_closed_form (K : ℕ → ℤ) (h0 : K 0 = 1) (h1 : K 1 = 4)
    (hrec : ∀ k, K (k + 2) = 5 * K (k + 1) - 4 * K k) :
    ∀ k, K k = 4 ^ k := by
  have key : ∀ k, K k = 4 ^ k ∧ K (k + 1) = 4 ^ (k + 1) := by
    intro k
    induction k with
    | zero => exact ⟨by rw [h0, pow_zero], by rw [h1, pow_one]⟩
    | succ n ih =>
        refine ⟨ih.2, ?_⟩
        rw [hrec n, ih.1, ih.2, pow_succ, pow_succ]
        ring
  exact fun k => (key k).1

-- Corollary (not mechanized here — `Nat.size` module not in the slice): `bit(K_k) = bit(4^k)
-- = 2k+1`, linear in the horizon.  The closed form above is the substantive content.

#print axioms continuant_closed_form
