import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-
Recovery-rate bound (note `thm:recovery`, the `O(1/D²)` leg — an explicit bound, not asymptotic).
At a nonsingular interior dual maximizer the dual value is locally a concave quadratic
`g(λ) = g* − k(λ−λ*)²` (`k ≥ 0`, `g* = J*`).  Rationalizing the optimal multiplier to precision
`1/D` gives `|λ − λ*| ≤ 1/D`, so the certified suboptimality gap shrinks QUADRATICALLY in the
budget: `J* − g(λ) ≤ k / D²`.  Over ℚ, no `sorry`.
-/

/-- **Quadratic recovery rate.** For a concave-quadratic dual `g(λ) = g* − k(λ−λ*)²` with `k ≥ 0`
and a multiplier within `1/D` of the maximizer `λ*`, the gap `g* − g(λ)` is at most `k/D²`
(quadratic shrinkage in the rationalization budget `D`). -/
theorem recovery_rate (g : ℚ → ℚ) (gstar lstar k D : ℚ) (hk : 0 ≤ k) (hD : 0 < D)
    (hg : ∀ l, g l = gstar - k * (l - lstar) ^ 2)
    (l : ℚ) (hl : |l - lstar| ≤ 1 / D) :
    gstar - g l ≤ k / D ^ 2 := by
  rw [hg l]
  have hsq : (l - lstar) ^ 2 ≤ (1 / D) ^ 2 := by
    nlinarith [le_abs_self (l - lstar), neg_abs_le (l - lstar), hl]
  have hmul : k * (l - lstar) ^ 2 ≤ k * (1 / D) ^ 2 := mul_le_mul_of_nonneg_left hsq hk
  have hd2 : k * (1 / D) ^ 2 = k / D ^ 2 := by rw [div_pow, one_pow, mul_one_div]
  linarith [hmul, hd2]

/-- **Linear recovery rate (singular / trust-region hard case, note `thm:hard`).** When the dual
value has a kink at the maximizer — `g(λ) = g* − k|λ−λ*|` (a one-sided derivative, `k ≥ 0`) — the
gap shrinks only LINEARLY in the budget: `J* − g(λ) ≤ k/D`. -/
theorem recovery_rate_hard (g : ℚ → ℚ) (gstar lstar k D : ℚ) (hk : 0 ≤ k) (hD : 0 < D)
    (hg : ∀ l, g l = gstar - k * |l - lstar|)
    (l : ℚ) (hl : |l - lstar| ≤ 1 / D) :
    gstar - g l ≤ k / D := by
  rw [hg l]
  have hmul : k * |l - lstar| ≤ k * (1 / D) := mul_le_mul_of_nonneg_left hl hk
  have hd : k * (1 / D) = k / D := mul_one_div k D
  linarith [hmul, hd]

#print axioms recovery_rate
#print axioms recovery_rate_hard
