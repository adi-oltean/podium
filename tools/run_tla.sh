#!/usr/bin/env bash
# Run the complete TLA+ lane locally: every shipped TLC configuration,
# the E-NAV falsification run (which MUST fail), and the annotation
# extraction gate. Mirrors .github/workflows/tla.yml.
#
# Downloads tla2tools.jar (pinned release, checksum-verified) into the
# repo root on first use; the jar is gitignored. Only java (11+), curl,
# and python3 are required. Trace validation against real simulation
# runs is separate (tools/tla_trace_check.py; needs podium installed).
#
# Usage: tools/run_tla.sh

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
JAR="$ROOT/tla2tools.jar"
TLC_URL="https://github.com/tlaplus/tlaplus/releases/download/v1.8.0/tla2tools.jar"
TLC_SHA256="ab323b79802aedc3203b3f9af37c6aca3ed43f4e0225b36f2aa77b26de46c05f"
METADIR="$(mktemp -d "${TMPDIR:-/tmp}/tlc-meta.XXXXXX")"
LOG="$METADIR/tlc.log"
trap 'rm -rf "$METADIR" "$ROOT"/tla/*_TTrace_* 2>/dev/null || true' EXIT

sha_check() {  # portable sha256 verification (linux: sha256sum, macOS: shasum)
    if command -v sha256sum >/dev/null 2>&1; then
        echo "$TLC_SHA256  $JAR" | sha256sum -c - >/dev/null
    else
        echo "$TLC_SHA256  $JAR" | shasum -a 256 -c - >/dev/null
    fi
}

if [ ! -f "$JAR" ]; then
    echo "fetching TLC (pinned v1.8.0) -> $JAR"
    curl -fsSL -o "$JAR" "$TLC_URL"
fi
sha_check || { echo "FATAL: tla2tools.jar checksum mismatch"; exit 1; }

tlc() {  # tlc <config> <module>: -deadlock DISABLES deadlock reporting
    java -cp "$JAR" tlc2.TLC -metadir "$METADIR" -workers auto -deadlock \
        -config "$ROOT/tla/$1.cfg" "$ROOT/tla/$2.tla"
}

fail=0
echo "== model checking (all shipped configurations) =="
while read -r cfg mod; do
    if out=$(tlc "$cfg" "$mod" 2>&1) && grep -q "No error has been found" <<<"$out"; then
        states=$(grep -oE "[0-9,]+ distinct states found" <<<"$out" | tail -1)
        printf 'PASS  %-22s %s\n' "$cfg" "($states)"
    else
        printf 'FAIL  %-22s\n' "$cfg"
        tail -30 <<<"$out"
        fail=1
    fi
done <<'EOF'
ArchRendezvousSRA ArchRendezvous
ArchRendezvousSRNA ArchRendezvous
ArchRendezvousExit ArchRendezvous
Mission Mission
NavApp NavApp
NavAppBackpressure NavApp
EOF

echo "== E-NAV falsification receipt (this run MUST fail) =="
if out=$(tlc NavAppReader NavApp 2>&1); then
    echo "FAIL  NavAppReader unexpectedly verified — the atomicity obligation lost its witness"
    fail=1
elif grep -q "Invariant NoPartialPublish is violated" <<<"$out"; then
    echo "PASS  NavAppReader rejected: NoPartialPublish counterexample produced, as required"
else
    echo "FAIL  NavAppReader failed for the wrong reason:"
    tail -30 <<<"$out"
    fail=1
fi
rm -f "$ROOT"/tla/*_TTrace_*  # counterexample artifacts of the expected failure

echo "== annotation extraction gate =="
if python3 "$ROOT/tools/tla_extract.py" --strict; then
    echo "PASS  extraction"
else
    echo "FAIL  extraction"
    fail=1
fi

[ "$fail" -eq 0 ] && echo "TLA LANE PASSED" || echo "TLA LANE FAILED"
exit "$fail"
