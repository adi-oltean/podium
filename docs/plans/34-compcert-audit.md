# 34 — CompCert audit: golden vectors through the verified compiler

GitHub issue: https://github.com/adi-oltean/podium/issues/34

## Shape of the claim

Tier-1 before this issue: "gcc with pinned FP flags reproduces CPython
bit-for-bit on the emitted kernels." After: the same golden vectors
replay through **ccomp** — CompCert, the C compiler with a
machine-checked semantics-preservation proof — so the binary's
behavior is connected to the C source's formal semantics by proof,
and to the Python semantics by the bitwise receipts. That is the
strongest form of "the flight code means what the Python meant" short
of a full source-level proof.

## Pieces

- `tests/test_cemit.py` gains PODIUM_CC: the golden fixture compiles
  with any compiler; ccomp gets `-O` (it never contracts FP and takes
  no -std/-Werror), gcc keeps the pinned flags.
- Subset tripwire test: emitted C must never contain goto / switch /
  union / long double / _Complex / malloc / alloca / setjmp / va_ /
  asm, and every declared array dimension is a decimal literal (no
  VLAs). Future emitter growth cannot silently leave the verified-
  compilable subset.
- `.github/workflows/compcert.yml`: opam-installed ccomp (setup-ocaml
  caches; first run builds CompCert ~25 min), then the full
  test_cemit suite runs with PODIUM_CC=ccomp — all 20 kernels'
  vectors through the verified compiler, on emit paths + weekly.
- Local bring-up: ocaml/opam Docker image (sandbox has no opam);
  the image's opam index is stale — `opam update` required before
  `opam install compcert` (first attempt failed on this).

## Results

Local (ocaml/opam container): **CompCert 3.17 compiles all 20 emitted
kernels cleanly** (`ccomp -c kernels.c` — COMPILE_OK). The bring-up
hit four real snags, each now encoded in build script + workflow so
nobody rediscovers them: (1) Docker Desktop's WSL credential helper
breaks anonymous pulls (clean DOCKER_CONFIG); (2) the ocaml/opam
image's default repo is a FROZEN local git clone that never fetches;
(3) `opam repo add --set-default` does not attach the repo to the
ACTIVE switch (plain add does); (4) the package is `coq-compcert` in
the coq-released archive — the old standalone `compcert` name no
longer exists upstream. The full golden-vector replay through ccomp
is compcert.yml's job (first CI run builds Coq+CompCert ~35 min,
cached thereafter).

## Push/merge instructions

Single commit on main: `34 — CompCert audit (#34)`; push; watch
compcert.yml's first run; close after green.
