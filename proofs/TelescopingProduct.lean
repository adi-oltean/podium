import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-
Mechanized certificate of the round-2 review correction (docs/deeper-angles.md §2).

The tridiagonal LDLᵀ pivots are d_j = K_{j+1} / K_j, ratios of consecutive continuants
(leading minors) K_j.  An earlier draft claimed the pivot PRODUCT ∏_{j<k} d_j had bit-size
Σ bit(d_j) = Θ(k²) — the fallacy bit(∏) = Σ bit.  In fact the product TELESCOPES to the
single continuant K_k, whose bit-size is Θ(k), so fraction-free elimination is NOT
quadratically worse.  Here is the telescoping identity, machine-checked over ℚ.
-/

/-- Recursive product of the first k pivots d_j = K_{j+1}/K_j (avoids BigOperators). -/
def pivotProd (K : ℕ → ℚ) : ℕ → ℚ
  | 0 => 1
  | k + 1 => pivotProd K k * (K (k + 1) / K k)

/-- Telescoping: the product of the first k pivots is exactly the continuant ratio
`K k / K 0` — one leading minor, not a growing product of pivot bit-sizes. -/
theorem pivotProd_telescopes (K : ℕ → ℚ) (h : ∀ j, K j ≠ 0) :
    ∀ k, pivotProd K k = K k / K 0 := by
  intro k
  induction k with
  | zero => simp [pivotProd, div_self (h 0)]
  | succ n ih =>
      have hn := h n
      have h0 := h 0
      unfold pivotProd
      rw [ih]
      field_simp

#print axioms pivotProd_telescopes
