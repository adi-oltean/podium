# 49 — Documentation consolidation

GitHub issue: https://github.com/adi-oltean/podium/issues/49

## Fix (landed) — no code changes

The flagship docs had fallen far behind the code after the v0.3–v0.7
feature burst.

- `docs/verification.md`: rewritten from PLAN-tense to SHIPPED-reality.
  New "What is shipped" table maps the nine verification modalities
  (contracts, STL oracles, closed-loop reachability, exact barrier +
  KKT certificates, golden vectors, sound EVA static analysis,
  verified-compiler + cross-arch, independent physics/dynamics oracles)
  to their module, CI lane, and receipts. The "Python → C translation"
  and "Layered assurance" sections now describe the emitter, CORE-MATH,
  EVA, CompCert, qemu, and the KKT checker as shipped, keeping the
  design rationale and precedents. The contract-pipeline table's C-
  translation row now describes ACSL rendering discharged by EVA.
- `README.md`: status line updated (v0.1–v0.6 complete, v0.7 current;
  eight CI lanes); the "external abstract-interpretation tool" framing
  replaced by the integrated Frama-C/EVA + CompCert story; the pillar
  list, layout tree, module descriptions, design principles, and docs
  table brought current (attitude no longer "planned", SCP/6-DOF/
  allocation/environmental-torque/emitter/certificates all listed).
- `CONTRIBUTING.md`: removed the docs-upkeep row for the nonexistent
  `docs/trajectory-optimization.md`.

## Verification

Docs-only; the test suite is unaffected. Checked that no stale phrases
remain ("v0.3 in progress", "v0.4 layer", "(planned)", "external
prover/validation tool") and that every doc referenced in the README
table and the layout tree exists.

## Push/merge instructions

Single commit on main: `49 — Documentation consolidation (#49)`;
push; close.
