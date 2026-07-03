# 04 — References cleanup and arXiv preprint

GitHub issue: https://github.com/adi-oltean/podium/issues/4

## Problem

docs/trajectory-optimization.md ends in a run-on citation blob — hard to
scan, no links, inconsistent formats. And the project now has enough
validated technical content (YA kernel + envelopes, truth model, pulsed
docking control with two empirically-caught failure modes, black-box sysid
case study) to warrant an arXiv preprint.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `docs/trajectory-optimization.md` | PARTIAL | replace reference blob with formatted list |
| `paper/main.tex`, `paper/references.bib` | MISS | arXiv-ready LaTeX source |
| `paper/podium-preprint.pdf` | MISS | compiled artifact (tectonic) |
| `README.md` | PARTIAL | link the preprint |

## Fix

1. Rewrite the Key References section as a numbered, linked list (authors,
   title, venue, year, DOI/arXiv), consistent style.
2. `paper/`: single-column article-class preprint (~10 pages):
   introduction & ecosystem gap; two-layer verification-oriented
   architecture (static subset, contracts); relative-motion kernels with
   the truth-model validation methodology and measured validity envelopes
   (C(e) constants, quadratic-scaling receipts); nonlinear truth model
   (exact LVLH transform incl. frame-precession term, J2/drag receipts);
   pulsed terminal docking control — deadband stall offset tau*dv/2,
   near-zone gain scheduling, sub-pulse dithering, and the frame-rate-
   dependent reincarnation of the stall; ISS-sim autopilot as a black-box
   system-identification and validation case study; planned convex
   guidance stack (clearly marked as design); verification pipeline and
   precedents. Honest v0.1 framing throughout.
3. Compile with tectonic (self-contained binary; no TeX install on box).
4. Author: Adi Oltean; acknowledgment notes AI-assisted development
   (arXiv norms). Category suggestion: eess.SY, cross-list cs.SE.

## Tests

`node tmp/ro/check_issim.mjs` unaffected; PDF builds clean (no undefined
references/citations); existing suite untouched.

## Acceptance Criteria

- [x] References section reformatted with links (21 entries, grouped)
- [x] paper/ builds to PDF without errors (tectonic 0.16.9; two cosmetic
      overfull-hbox warnings only)
- [x] README links the preprint
- [x] Committed and pushed; issue closed

## Push/merge instructions

Single commit on main: `04 — References cleanup and arXiv preprint (#4)`,
push, close #4.

## Verification steps

Open the PDF; check all citations resolve; spot-check envelope numbers
against tests/test_validity_envelopes.py and the sysid logs in issue #3.
