#!/usr/bin/env bash
# Regression: type-check every Lean proof file in this directory against the pinned
# Mathlib project, and print each theorem's axiom dependencies (expect only the three
# standard axioms propext/Classical.choice/Quot.sound -- or none -- and never sorryAx).
#
# Usage:  bash setup.sh      # once, builds the pinned project into .leanwork/
#         bash check-all.sh  # type-check everything, audit axioms
#
# Requires elan/lake and the pinned Mathlib project that setup.sh builds from
# ci/lakefile.toml + ci/lake-manifest.json. Set LEANWORK to point at an existing
# build elsewhere, and LAKE to a lake binary that is not on PATH.
set -u
LP="$(cd "$(dirname "$0")" && pwd)"
WD="${LEANWORK:-$LP/.leanwork}"
LAKE_BIN="${LAKE:-$(command -v lake || echo "$HOME/.elan/bin/lake")}"
LEAN="$LAKE_BIN env lean"

if [ ! -d "$WD" ]; then
  echo "no Mathlib project at $WD -- run: bash $LP/setup.sh" >&2
  exit 1
fi
TMP="$(mktemp)"
fail=0
for f in \
  "$LP/DareFixedPoint.lean" \
  "$LP/SProcedureSoundness.lean" \
  "$LP/LiftingIdentity.lean" \
  "$LP/CdtLifting.lean" \
  "$LP/GeneralLifting.lean" \
  "$LP/DiagDominance.lean" \
  "$LP/TelescopingProduct.lean" \
  "$LP/ContinuantClosedForm.lean" \
  "$LP/RecoveryRate.lean" \
  "$LP/MatrixDare.lean" \
  "$LP/ClosureOptimum.lean" \
  "$LP/LDLTSound.lean" \
  "$LP/StructuralPSD.lean" \
  "$LP/GershgorinPSD.lean" \
  "$LP/Certificate.lean" \
  "$LP/SchurStep.lean" \
  "$LP/BlockSchur.lean" \
  "$LP/PosSemidefBridge.lean" \
  "$LP/note/PodiumNote.lean"
do
  base="$(basename "$f")"
  if env -C "$WD" $LEAN "$f" >"$TMP" 2>&1
  then
    echo "PASS  $base"
    # Surface the axiom audit on the PASS path too: the #print axioms lines
    # are the attestation that only propext/Classical.choice/Quot.sound are
    # admitted and sorryAx never is. Discarding them on success would make
    # that attestation invisible on exactly the runs that matter.
    grep "depends on axioms" "$TMP" || true
  else
    echo "FAIL  $base"
    fail=1
    cat "$TMP"
  fi
done
rm -f "$TMP"
echo "----"
if [ "$fail" = 0 ]
then
  echo "ALL PASS"
else
  echo "SOME FAILED"
  exit 1
fi
