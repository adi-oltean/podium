# Verification-modality coverage

Podium runs several independent verification lanes. Are they
complementary — each catching failure modes the others miss — or
redundant? `tools/fault_coverage.py` answers it empirically: it injects
a representative fault into each artifact class and records which lane
catches it. `tests/test_fault_coverage.py` asserts the properties
below so the answer stays true as the code evolves.

## The coverage matrix

`[x]` = the lane rejects the faulted artifact; `.` = the lane passes
(does not fire). Every lane is also run against the *good* artifacts of
the other faults; none fires (zero false alarms).

| fault (artifact perturbed) | conserv | analytic | golden | STL | barrier | KKT | Lyap | SOS |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| dynamics torque bug        | [x] | [x] |  .  |  .  |  .  |  .  |  .  |  .  |
| emitter translation        |  .  |  .  | [x] |  .  |  .  |  .  |  .  |  .  |
| spec corridor breach       |  .  |  .  |  .  | [x] |  .  |  .  |  .  |  .  |
| invalid barrier cert       |  .  |  .  |  .  |  .  | [x] |  .  |  .  |  .  |
| invalid KKT solution       |  .  |  .  |  .  |  .  |  .  | [x] |  .  |  .  |
| invalid Lyapunov P         |  .  |  .  |  .  |  .  |  .  |  .  | [x] |  .  |
| non-PSD SOS Gram           |  .  |  .  |  .  |  .  |  .  |  .  |  .  | [x] |

## What it shows

- **Six of seven faults are caught by exactly one lane.** Remove that
  lane and the fault ships silently. The stack is not redundant.
- **Certificate faults are invisible to the physics lanes.** A wrong
  *proof* (a corrupted barrier/KKT/Lyapunov/SOS certificate) does not
  violate energy, momentum, or the analytic trajectory — only the exact
  arithmetic checker for that certificate rejects it. Conservation laws
  and the Jacobi-elliptic oracle catch none of them.
- **The emitter fault is invisible to everything but golden vectors.** A
  Python↔C transcription bug is undetectable by a Python reference that
  is itself correct; only bit-for-bit cross-checking of the emitted C
  catches it.
- **A dynamics bug lights up multiple physics lanes.** An energy-
  injecting torque error breaks both conservation and the analytic
  oracle — the physics lanes overlap on dynamics faults, so they
  reinforce rather than substitute there.

The block structure is the justification for the multi-modal stack:
the lanes partition the fault space, and the exact-arithmetic
certificate checkers occupy a region — invalid proofs — that no
physics- or trajectory-level check reaches.

## Reproduce

```
python3 tools/fault_coverage.py            # prints the matrix + metrics
python3 tools/fault_coverage_figure.py OUT # renders OUT.svg + OUT.pdf
pytest tests/test_fault_coverage.py        # asserts the properties
```

`tools/fault_coverage_figure.py` renders the matrix as a publication
figure (grouped by lane kind — physics / emitter / spec / exact
certificates) directly from `build_matrix()`, so the figure can never
drift from the asserted data.
