#!/usr/bin/env python3
"""Build the release audit bundle — and GATE on it.

Produces in <outdir>:
  bundle.json     byte-deterministic reference-mission audit (fixed
                  seed): capture outcome, IDSS margins (translation +
                  rotation), STL phase margins, dv, the exact-rational
                  barrier verdict. No timestamps, no host facts — the
                  same source must produce the same bytes.
  kernels.c       the emitted flight C, exactly what EVA verifies
  eva_driver.c    the interval-input driver (contract gaps visible)
  meta.json       tag/commit/toolchain stamps — SEPARATE from
                  bundle.json so the latter stays byte-comparable
                  across rebuilds of identical code
  SHA256SUMS      sha-256 of the byte-deterministic files (kernels.c,
                  eva_driver.c, bundle.json) in `sha256sum -c` format, so
                  a third party can confirm an identical-source rebuild

Hard gate: exits nonzero unless the mission CAPTURES, every IDSS
margin is positive, and the barrier certificate verifies. A tag cannot
become a release on failed evidence. (The EVA zero-alarm gate runs as
its own release step via tools/eva_gate.py.)

Usage: python3 tools/build_audit_bundle.py <outdir> [--tag TAG] [--sha SHA]
"""

import argparse
import hashlib
import json
import pathlib
import platform
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy  # noqa: E402

from podium.emit import cemit, evagen  # noqa: E402
from podium.emit.kernels import FLIGHT_KERNELS  # noqa: E402
from podium.sim import mission  # noqa: E402

SEED = 7


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--tag", default="untagged")
    ap.add_argument("--sha", default="unknown")
    args = ap.parse_args()
    out = pathlib.Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    (out / "kernels.c").write_text(cemit.emit_module(FLIGHT_KERNELS))
    (out / "eva_driver.c").write_text(
        evagen.emit_eva_driver(FLIGHT_KERNELS))

    print(f"flying the reference mission (seed {SEED}) ...", flush=True)
    res = mission.fly(seed=SEED)
    (out / "bundle.json").write_text(mission.audit_bundle(res, SEED) + "\n")

    meta = {
        "tag": args.tag,
        "commit": args.sha,
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "kernels_emitted": len(FLIGHT_KERNELS),
    }
    (out / "meta.json").write_text(json.dumps(meta, indent=1,
                                              sort_keys=True) + "\n")

    # SHA-256 manifest over the byte-deterministic files only (identical
    # source -> identical bytes -> identical digests). meta.json is excluded
    # because it carries variable toolchain stamps by design.
    manifest = "".join(
        f"{hashlib.sha256((out / name).read_bytes()).hexdigest()}  {name}\n"
        for name in ("kernels.c", "eva_driver.c", "bundle.json"))
    (out / "SHA256SUMS").write_text(manifest)

    failures = []
    if not res.captured:
        failures.append("mission did not capture")
    for k, v in {**res.idss_translation, **res.idss_rotation}.items():
        if v <= 0.0:
            failures.append(f"IDSS margin {k} = {v}")
    for k, v in res.spec_margins.items():
        if v <= 0.0:
            failures.append(f"STL margin {k} = {v}")
    if not res.barrier_ok:
        failures.append("barrier certificate failed exact verification")
    if failures:
        print("AUDIT GATE FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"AUDIT GATE PASSED: captured at t={res.contact_time:.0f}s, "
          f"dv={res.dv_total:.2f} m/s, barrier verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
