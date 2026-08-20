import Mathlib.Data.Matrix.Mul
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp

/-
Inductive engine of the checker COMPLETENESS leg (paper Lemma band, the hard direction):
one step of symmetric (banded) Gaussian / Schur elimination preserves PSD and yields a
nonnegative pivot.  Concretely, for the bordered matrix `[[a, bᵀ], [b, C]]`:
  • `IsPSD` ⇒ pivot `a ≥ 0`                                  (`headBordered_pivot_nonneg`);
  • `IsPSD` and `a > 0` ⇒ the Schur complement `C − a⁻¹ bbᵀ` is PSD  (`schur_complement_psd`),
    obtained by minimizing the quadratic form over the eliminated coordinate `x = −bᵀw/a`.
Iterating this over the head coordinate is exactly the exact band-`LDLᵀ` sweep: it proves
`IsPSD ⇒ LᵀDL with pivots ≥ 0` (completeness), the converse of `LDLTSound.ldlt_isPSD`.
Over any `LinearOrderedField`, no `sorry`.
-/

open Matrix
variable {n : ℕ} {K : Type*} [Field K] [LinearOrder K] [IsStrictOrderedRing K]

/-- PSD predicate. -/
def IsPSD {m : ℕ} (M : Matrix (Fin m) (Fin m) K) : Prop := ∀ v : Fin m → K, 0 ≤ v ⬝ᵥ M *ᵥ v

/-- Bordered matrix `[[a, bᵀ], [b, C]]` with pivot `a` at the head. -/
def headBordered (a : K) (b : Fin n → K) (C : Matrix (Fin n) (Fin n) K) :
    Matrix (Fin (n + 1)) (Fin (n + 1)) K :=
  Matrix.of (Fin.cons (Fin.cons a b) (fun i => Fin.cons (b i) (C i)))

/-- Quadratic form of the bordered matrix at `(x, w)`: `a x² + 2x (bᵀw) + wᵀ C w`. -/
theorem headBordered_form (a : K) (b : Fin n → K) (C : Matrix (Fin n) (Fin n) K)
    (x : K) (w : Fin n → K) :
    (Fin.cons x w) ⬝ᵥ (headBordered a b C) *ᵥ (Fin.cons x w)
      = a * x ^ 2 + 2 * x * (b ⬝ᵥ w) + w ⬝ᵥ C *ᵥ w := by
  simp only [dotProduct, mulVec, headBordered, Matrix.of_apply, Fin.sum_univ_succ,
    Fin.cons_zero, Fin.cons_succ]
  rw [show (∑ i, w i * (b i * x + ∑ j, C i j * w j))
        = ∑ i, (x * (b i * w i) + w i * (∑ j, C i j * w j)) from
        Finset.sum_congr rfl fun i _ => by ring,
      Finset.sum_add_distrib, ← Finset.mul_sum]
  ring

/-- Pivot nonnegativity: `IsPSD` of the bordered matrix forces `a ≥ 0` (test the head basis
vector `(1, 0)`). -/
theorem headBordered_pivot_nonneg (a : K) (b : Fin n → K) (C : Matrix (Fin n) (Fin n) K)
    (h : IsPSD (headBordered a b C)) : 0 ≤ a := by
  have := h (Fin.cons 1 0)
  rw [headBordered_form] at this
  simpa using this

/-- **Schur complement preserves PSD (the elimination step).** If the bordered matrix is PSD
and the pivot `a > 0`, then the Schur complement `C − a⁻¹ bbᵀ` is PSD.  Proof: for any tail
vector `w`, evaluate the (nonnegative) form at the minimizing head coordinate `x = −bᵀw/a`;
the minimum value is exactly `wᵀ(C − a⁻¹ bbᵀ)w`. -/
theorem schur_complement_psd (a : K) (b : Fin n → K) (C : Matrix (Fin n) (Fin n) K)
    (h : IsPSD (headBordered a b C)) (ha : 0 < a) :
    IsPSD (C - a⁻¹ • vecMulVec b b) := by
  intro w
  -- value of the Schur-complement form
  have hvmv : w ⬝ᵥ (vecMulVec b b) *ᵥ w = (b ⬝ᵥ w) * (b ⬝ᵥ w) := by
    simp only [dotProduct, mulVec, vecMulVec_apply]
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.mul_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by ring
  have hSform : w ⬝ᵥ (C - a⁻¹ • vecMulVec b b) *ᵥ w
      = w ⬝ᵥ C *ᵥ w - a⁻¹ * (b ⬝ᵥ w) ^ 2 := by
    simp only [sub_mulVec, dotProduct_sub, smul_mulVec, dotProduct_smul, smul_eq_mul]
    rw [hvmv]; ring
  -- evaluate PSD hypothesis at the minimizer x = -(bᵀw)/a
  have hx := h (Fin.cons (-(b ⬝ᵥ w) / a) w)
  rw [headBordered_form] at hx
  have ha' : a ≠ 0 := ne_of_gt ha
  -- the form at the minimizer equals the Schur-complement value
  have hmin : a * (-(b ⬝ᵥ w) / a) ^ 2 + 2 * (-(b ⬝ᵥ w) / a) * (b ⬝ᵥ w) + w ⬝ᵥ C *ᵥ w
      = w ⬝ᵥ C *ᵥ w - a⁻¹ * (b ⬝ᵥ w) ^ 2 := by
    field_simp
    ring
  rw [hmin] at hx
  rw [hSform]
  exact hx

/-- **Degenerate pivot ⇒ zero border.** If the bordered matrix is PSD and its pivot is `0`,
the border `b` vanishes.  (A zero pivot with `b ≠ 0` would make the form linear and unbounded
below in the eliminated coordinate.)  This is the fact that lets the completeness induction pass
through a rank-deficient step: at `a = 0` the exact `LDLᵀ` records a zero pivot and the row
decouples. -/
theorem headBordered_border_zero (b : Fin n → K) (C : Matrix (Fin n) (Fin n) K)
    (h : IsPSD (headBordered 0 b C)) : ∀ w, b ⬝ᵥ w = 0 := by
  intro w
  by_contra hc
  -- the form is linear in the eliminated coordinate x, hence unbounded below unless bᵀw = 0
  have key : ∀ x : K, 0 ≤ 2 * x * (b ⬝ᵥ w) + w ⬝ᵥ C *ᵥ w := by
    intro x
    have hx := h (Fin.cons x w)
    rw [headBordered_form] at hx
    simpa using hx
  have hx := key (-(w ⬝ᵥ C *ᵥ w + 1) / (2 * (b ⬝ᵥ w)))
  have hval : 2 * (-(w ⬝ᵥ C *ᵥ w + 1) / (2 * (b ⬝ᵥ w))) * (b ⬝ᵥ w) + w ⬝ᵥ C *ᵥ w = -1 := by
    field_simp
    ring
  rw [hval] at hx
  linarith

/-- Rank-1 matrices `vvᵀ` are PSD (the atoms of any `LDLᵀ`): `wᵀ(vvᵀ)w = (vᵀw)² ≥ 0`. -/
theorem rank1_isPSD {m : ℕ} (v : Fin m → K) : IsPSD (vecMulVec v v) := by
  intro w
  have hvw : w ⬝ᵥ (vecMulVec v v) *ᵥ w = (v ⬝ᵥ w) * (v ⬝ᵥ w) := by
    simp only [dotProduct, mulVec, vecMulVec_apply]
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.mul_sum, Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [hvw]
  exact mul_self_nonneg _

/-- Embedding a tail block as `[[0,0],[0,S]]` preserves PSD (the head coordinate is inert). -/
theorem embed_isPSD (S : Matrix (Fin n) (Fin n) K) (h : IsPSD S) :
    IsPSD (headBordered 0 0 S) := by
  intro v
  rw [← Fin.cons_self_tail v, headBordered_form]
  simpa using h (Fin.tail v)

/-- **Rank-1 + Schur decomposition of the bordered matrix** (`a ≠ 0`): the elimination step,
made explicit as an equality. `[[a,bᵀ],[b,C]] = a⁻¹·(a;b)(a;b)ᵀ + [[0,0],[0, C−a⁻¹bbᵀ]]`.
Iterating this — each step peeling a nonnegative-weighted rank-1 `dⱼ vⱼvⱼᵀ` (`dⱼ = pivot ≥ 0`
by `headBordered_pivot_nonneg`, Schur complement PSD by `schur_complement_psd`) — is exactly the
`LDLᵀ` factorization `M = ∑ⱼ dⱼ vⱼvⱼᵀ`, the completeness leg. -/
theorem headBordered_decomp (a : K) (b : Fin n → K) (C : Matrix (Fin n) (Fin n) K) (ha : a ≠ 0) :
    headBordered a b C
      = a⁻¹ • vecMulVec (Fin.cons a b) (Fin.cons a b)
        + headBordered 0 0 (C - a⁻¹ • vecMulVec b b) := by
  ext i j
  refine Fin.cases ?_ (fun i' => ?_) i <;> refine Fin.cases ?_ (fun j' => ?_) j <;>
    simp only [headBordered, Matrix.of_apply, Matrix.add_apply, Matrix.smul_apply,
      vecMulVec_apply, Matrix.sub_apply, Fin.cons_zero, Fin.cons_succ, smul_eq_mul,
      Pi.zero_apply, mul_zero, zero_mul, add_zero, zero_add]
  · field_simp
  · field_simp
  · field_simp
  · ring

/-! ### Assembly infrastructure: `embed = headBordered 0 0 ·` is linear and turns rank-1 into
rank-1, plus the decomposition of an arbitrary symmetric matrix into `headBordered`. -/

/-- `embed` is additive. -/
theorem embed_add (X Y : Matrix (Fin n) (Fin n) K) :
    headBordered 0 0 (X + Y) = headBordered 0 0 X + headBordered 0 0 Y := by
  ext i j
  refine Fin.cases ?_ (fun i' => ?_) i <;> refine Fin.cases ?_ (fun j' => ?_) j <;>
    simp [headBordered]

/-- `embed` commutes with scaling. -/
theorem embed_smul (c : K) (X : Matrix (Fin n) (Fin n) K) :
    headBordered 0 0 (c • X) = c • headBordered 0 0 X := by
  ext i j
  refine Fin.cases ?_ (fun i' => ?_) i <;> refine Fin.cases ?_ (fun j' => ?_) j <;>
    simp [headBordered]

/-- `embed` of a rank-1 `wwᵀ` is the rank-1 of the head-padded vector `(0;w)`. -/
theorem embed_rank1 (w : Fin n → K) :
    headBordered 0 0 (vecMulVec w w) = vecMulVec (Fin.cons 0 w) (Fin.cons 0 w) := by
  ext i j
  refine Fin.cases ?_ (fun i' => ?_) i <;> refine Fin.cases ?_ (fun j' => ?_) j <;>
    simp [headBordered, vecMulVec_apply]

/-- `embed` commutes with finite sums (from `embed_add` by induction on the index set). -/
theorem embed_sum {k : ℕ} (f : Fin k → Matrix (Fin n) (Fin n) K) :
    headBordered 0 0 (∑ j, f j) = ∑ j, headBordered 0 0 (f j) := by
  classical
  have gen : ∀ s : Finset (Fin k),
      headBordered 0 0 (∑ j ∈ s, f j) = ∑ j ∈ s, headBordered 0 0 (f j) := by
    intro s
    induction s using Finset.induction with
    | empty =>
        simp only [Finset.sum_empty]
        ext i j
        refine Fin.cases ?_ (fun i' => ?_) i <;> refine Fin.cases ?_ (fun j' => ?_) j <;>
          simp [headBordered]
    | insert a s h ih =>
        rw [Finset.sum_insert h, Finset.sum_insert h, embed_add, ih]
  exact gen Finset.univ

/-- Any symmetric `(n+1)×(n+1)` matrix is `headBordered` of its pivot, first row, and tail block. -/
theorem eq_headBordered (M : Matrix (Fin (n + 1)) (Fin (n + 1)) K)
    (hsymm : ∀ i j, M i j = M j i) :
    M = headBordered (M 0 0) (fun i => M 0 i.succ) (fun i j => M i.succ j.succ) := by
  ext i j
  refine Fin.cases ?_ (fun i' => ?_) i <;> refine Fin.cases ?_ (fun j' => ?_) j <;>
    simp only [headBordered, Matrix.of_apply, Fin.cons_zero, Fin.cons_succ]
  exact hsymm i'.succ 0

#print axioms headBordered_form
#print axioms headBordered_pivot_nonneg
#print axioms schur_complement_psd
#print axioms headBordered_border_zero
#print axioms rank1_isPSD
#print axioms embed_isPSD
#print axioms headBordered_decomp
/-- **Checker COMPLETENESS (the hard direction).** Every symmetric PSD matrix factors as a
nonnegative combination of rank-1s, `M = ∑ⱼ dⱼ vⱼvⱼᵀ` with `dⱼ ≥ 0` — i.e. `M = VᵀDV`, the
`LDLᵀ` with nonnegative pivots (the converse of `ldlt_isPSD`).  Proof by induction on the size,
peeling one pivot/rank-1 per step via `headBordered_decomp`, `schur_complement_psd`, and the
degenerate `headBordered_border_zero`.  This closes the completeness leg.  No `sorry`. -/
theorem ldlt_complete : ∀ {m : ℕ} (M : Matrix (Fin m) (Fin m) K),
    (∀ i j, M i j = M j i) → IsPSD M →
    ∃ (k : ℕ) (d : Fin k → K) (v : Fin k → Fin m → K),
      (∀ j, 0 ≤ d j) ∧ M = ∑ j, d j • vecMulVec (v j) (v j) := by
  intro m
  induction m with
  | zero =>
      intro M _ _
      exact ⟨0, Fin.elim0, Fin.elim0, fun j => j.elim0, by ext i; exact i.elim0⟩
  | succ n ih =>
      intro M hsymm hpsd
      have hM : M = headBordered (M 0 0) (fun i => M 0 i.succ) (fun i j => M i.succ j.succ) :=
        eq_headBordered M hsymm
      set a := M 0 0 with ha_def
      set b : Fin n → K := fun i => M 0 i.succ with hb_def
      set C : Matrix (Fin n) (Fin n) K := fun i j => M i.succ j.succ with hC_def
      have hpsd' : IsPSD (headBordered a b C) := hM ▸ hpsd
      have hCsymm : ∀ i j, C i j = C j i := fun i j => hsymm i.succ j.succ
      have ha : 0 ≤ a := headBordered_pivot_nonneg a b C hpsd'
      rcases eq_or_lt_of_le ha with ha0 | hapos
      · -- degenerate pivot a = 0 : the border vanishes and the row decouples
        have hpsd0 : IsPSD (headBordered 0 b C) := by rw [← ha0] at hpsd'; exact hpsd'
        have hb0 : b = 0 := by
          funext i
          have hbi := headBordered_border_zero b C hpsd0 (Pi.single i 1)
          rwa [dotProduct_single, mul_one] at hbi
        have hMe : M = headBordered 0 0 C := by rw [hM, ← ha0, hb0]
        have hCpsd : IsPSD C := by
          intro w
          have hval := (hMe ▸ hpsd) (Fin.cons 0 w)
          rw [headBordered_form] at hval
          simpa using hval
        obtain ⟨k', d', v', hd', hC⟩ := ih C hCsymm hCpsd
        refine ⟨k', d', fun j => Fin.cons 0 (v' j), hd', ?_⟩
        rw [hMe, hC, embed_sum]
        exact Finset.sum_congr rfl fun j _ => by rw [embed_smul, embed_rank1]
      · -- positive pivot a > 0 : peel a⁻¹·(a;b)(a;b)ᵀ and recurse on the Schur complement
        have ha' : a ≠ 0 := ne_of_gt hapos
        have hSsymm : ∀ i j, (C - a⁻¹ • vecMulVec b b) i j = (C - a⁻¹ • vecMulVec b b) j i := by
          intro i j
          simp only [Matrix.sub_apply, Matrix.smul_apply, vecMulVec_apply, smul_eq_mul]
          rw [hCsymm i j]; ring
        have hSpsd : IsPSD (C - a⁻¹ • vecMulVec b b) := schur_complement_psd a b C hpsd' hapos
        obtain ⟨k', d', v', hd', hS⟩ := ih (C - a⁻¹ • vecMulVec b b) hSsymm hSpsd
        refine ⟨k' + 1, Fin.cons a⁻¹ d',
          Fin.cons (Fin.cons a b) (fun j => Fin.cons 0 (v' j)), ?_, ?_⟩
        · refine Fin.cases ?_ (fun j' => ?_)
          · rw [Fin.cons_zero]; exact le_of_lt (inv_pos.2 hapos)
          · rw [Fin.cons_succ]; exact hd' j'
        · rw [hM, headBordered_decomp a b C ha', Fin.sum_univ_succ]
          simp only [Fin.cons_zero, Fin.cons_succ]
          congr 1
          rw [hS, embed_sum]
          exact Finset.sum_congr rfl fun j _ => by rw [embed_smul, embed_rank1]

/-- **Checker criterion, both directions as one iff (the decision-procedure spec).** A symmetric
matrix is PSD *iff* it admits a nonnegative-pivot `LDLᵀ` factorization `∑ⱼ dⱼvⱼvⱼᵀ`. Forward is
`ldlt_complete`; backward is `rank1_isPSD` summed with nonnegative weights. This is exactly what
the exact-ℚ checker decides: `accept` (a nonnegative-pivot factorization is found) ⇔ `PSD`. -/
theorem isPSD_iff_ldlt {m : ℕ} (M : Matrix (Fin m) (Fin m) K) (hsymm : ∀ i j, M i j = M j i) :
    IsPSD M ↔ ∃ (k : ℕ) (d : Fin k → K) (v : Fin k → Fin m → K),
      (∀ j, 0 ≤ d j) ∧ M = ∑ j, d j • vecMulVec (v j) (v j) := by
  constructor
  · intro h; exact ldlt_complete M hsymm h
  · rintro ⟨k, d, v, hd, rfl⟩
    have hadd : ∀ (X Y : Matrix (Fin m) (Fin m) K), IsPSD X → IsPSD Y → IsPSD (X + Y) := by
      intro X Y hX hY w; rw [add_mulVec, dotProduct_add]; exact add_nonneg (hX w) (hY w)
    have hz : IsPSD (0 : Matrix (Fin m) (Fin m) K) := by intro w; simp
    apply Finset.sum_induction (fun j => d j • vecMulVec (v j) (v j)) IsPSD
      (fun X Y => hadd X Y) hz
    intro j _ w
    rw [smul_mulVec, dotProduct_smul, smul_eq_mul]
    exact mul_nonneg (hd j) (rank1_isPSD (v j) w)

#print axioms embed_add
#print axioms embed_rank1
#print axioms embed_sum
#print axioms eq_headBordered
#print axioms ldlt_complete
#print axioms isPSD_iff_ldlt
