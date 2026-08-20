# Lean treatment — the exact-arithmetic certificate note

Mechanization of the technical note **"Exact-Rational Certificates in Podium"**
(`docs/exact-arithmetic-certificates/note.tex`), against the same pinned toolchain and
Mathlib as the rest of the corpus. Build:
```
bash ../setup.sh                                   # once, from proofs/
env -C ../.leanwork lake env lean PodiumNote.lean
```

## Done (`PodiumNote.lean`, all no `sorry`, over any `LinearOrderedField`)
| Note statement | Lean |
|---|---|
| **Thm 1 (Lower bound / `thm:sound`)** — `λ≥0 ∧ M⪰0 ⇒ t≤J*` | `note_thm_sound` |
| **Thm 4 (Certified duality gap / `thm:multi`)** — same for any `m` | `note_thm_multi` |
| bracket closing observation (`thm:multi`) — `J* ∈ [t, f0 x̄]`, gap bounds suboptimality | `note_bracket` |
| **Lemma (Schur complement / `lem:schur`)**, scalar `A=a>0` — `M(λ,t)⪰0 ⟺ t ≤ g = c − b²/(4a)` | `note_lem_schur` |

## Not mechanized (analytic — recovery RATES)
- `thm:recovery` (nonsingular): interior maximizer ⇒ `g'(λ*)=0` ⇒ `J*−g(λ)=O(1/D²)` by Taylor.
- `thm:hard` (singular / trust-region hard case): (a) soundness at a rank-deficient PSD `M` is
  just `thm:sound` (already mechanized, since the exact `LDLᵀ` accepts a zero pivot); (b) the
  linear `O(1/D)` rate from a one-sided derivative; (c) `J*` algebraic of degree > 1.
These are calculus (Taylor / one-sided derivatives / algebraicity), needing analysis Mathlib;
the SOUNDNESS content (a) is covered by `note_thm_sound`. The rates are stated, not mechanized.

## Honest scope
As with the main-paper lemmas: these machine-check **classical / elementary** facts (weak
duality, the Schur-complement PSD criterion). Value = a tolerance-free, `sorry`-free trusted
path, not new mathematics.
