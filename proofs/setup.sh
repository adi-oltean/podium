#!/usr/bin/env bash
# Build the pinned Mathlib project these proofs type-check against.
#
#   bash proofs/setup.sh          # install elan if absent, build the project
#   bash proofs/check-all.sh      # then: type-check every file, audit axioms
#
# Everything is pinned: the toolchain in ci/lean-toolchain, Mathlib's revision
# in ci/lakefile.toml, and the full dependency graph in ci/lake-manifest.json.
# The build lands in .leanwork/ next to this script (git-ignored) and is
# reused by check-all.sh; nothing is installed into the proofs directory
# itself, so the corpus stays a flat set of .lean files.
#
# Mathlib is large. The first run downloads a prebuilt cache (`lake exe cache
# get`) rather than compiling it, which takes a few minutes on a cold machine
# and seconds afterwards.
set -o errexit
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
PROJ="$HERE/.leanwork"

# elan is the Lean toolchain manager; it reads ci/lean-toolchain and fetches
# the exact compiler version pinned there.
if ! command -v lake >/dev/null 2>&1 && [ ! -x "$HOME/.elan/bin/lake" ]; then
  echo "installing elan (Lean toolchain manager)..."
  curl -sSf https://elan.lean-lang.org/elan-init.sh -o "$PROJ.elan-init.sh" \
    || curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh \
       -o "$PROJ.elan-init.sh"
  bash "$PROJ.elan-init.sh" -y --no-modify-path
  rm -f "$PROJ.elan-init.sh"
fi
if ! command -v lake >/dev/null 2>&1; then
  export PATH="$HOME/.elan/bin:$PATH"
fi

if [ ! -f "$PROJ/lake-manifest.json" ]; then
  echo "creating the pinned Mathlib project at $PROJ ..."
  mkdir -p "$PROJ/Leanwork"
  cp "$HERE/ci/lakefile.toml" "$PROJ/"
  cp "$HERE/ci/lake-manifest.json" "$PROJ/"
  cp "$HERE/ci/lean-toolchain" "$PROJ/"
  printf 'import Mathlib\n' > "$PROJ/Leanwork.lean"
  printf 'def hello := "leanwork"\n' > "$PROJ/Leanwork/Basic.lean"
fi

echo "fetching the prebuilt Mathlib cache (first run: a few minutes) ..."
( cd "$PROJ" && lake exe cache get )

echo
echo "ready. Type-check the corpus with:"
echo "  bash $HERE/check-all.sh"
