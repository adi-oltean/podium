#!/usr/bin/env python3
"""Deploy a verified, version-pinned build of the viewer pages.

Adopts the fermi deploy discipline (fermi tools/deploy.py) for Podium's
Actions-based Pages layout. Build identity is DERIVED, never hand-typed:

  * BUILD = (highest existing viewer/builds/b<N>/) + 1 — monotonic, you
    cannot forget to bump it. First build requires --first explicitly.
  * SHA   = git rev-parse HEAD with a clean-tree guard on the shipped
    sources, so the badge always links to exactly what shipped.

The sources carry BUILD = 0 / 'HEAD' as DEV sentinels. This script
injects the real values into the LIVE pages (viewer/index.html,
viewer/3d/index.html) and writes IMMUTABLE snapshots to
viewer/builds/b<N>/{index.html, 3d.html} (relative links rewritten for
the new depth; the vendored three.js is shared via the version-pinned
viewer/vendor/three-<ver>/ directory, so upgrades ADD a pinned dir and
never mutate what old builds reference). Old builds are never touched:
every deploy only adds a directory, so all versions stay live side by
side. The builds catalog (viewer/builds/index.html) is regenerated from
the manifest. Everything written is read back and re-verified; any
mismatch aborts loudly.

Usage:
  python3 tools/deploy_viewer.py            # derive, inject, snapshot, verify
  python3 tools/deploy_viewer.py --dry-run  # verify in memory, write nothing
  python3 tools/deploy_viewer.py --first    # allow build 1 (no prior builds)

After it runs: commit everything (the deploy commit), push — the Pages
workflow ships the whole viewer/ tree.
"""

import datetime
import json
import os
import re
import subprocess
import sys

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VIEWER = os.path.join(SRC, "viewer")
BUILDS = os.path.join(VIEWER, "builds")
MANIFEST = os.path.join(BUILDS, "manifest.json")
SHIPPED = ["viewer/index.html", "viewer/3d/index.html"]
REPO_URL = "https://github.com/adi-oltean/podium"

SHA_PAT = r"(HEAD|[0-9a-f]{40})"


def die(msg: str) -> None:
    sys.exit(f"ABORT: {msg}")


def sh(*args: str) -> str:
    return subprocess.run(args, capture_output=True, text=True, check=True).stdout.strip()


def derive_sha() -> str:
    dirty = sh("git", "-C", SRC, "status", "--porcelain", *SHIPPED)
    if dirty:
        die("shipped files have uncommitted changes — commit them first so "
            f"the badge SHA matches what ships:\n{dirty}")
    return sh("git", "-C", SRC, "rev-parse", "HEAD")


def next_build(allow_first: bool) -> int:
    mx = 0
    if os.path.isdir(BUILDS):
        for x in os.listdir(BUILDS):
            m = re.fullmatch(r"b(\d+)", x)
            if m and os.path.isdir(os.path.join(BUILDS, x)):
                mx = max(mx, int(m.group(1)))
    if mx == 0 and not allow_first:
        die("no existing builds/b<N>/ found — pass --first for the first build")
    return mx + 1


def read(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        return f.read()


def inject(html: str, build: int, sha: str, label: str) -> str:
    """Replace the three identity markers; each must occur exactly once."""
    pats = [
        (r"const BUILD = \d+", f"const BUILD = {build}"),
        (r">build \d+</a>", f">build {build}</a>"),
        (rf"const SRC_COMMIT = '{SHA_PAT}'", f"const SRC_COMMIT = '{sha}'"),
        (rf"commit/{SHA_PAT}\"", f'commit/{sha}"'),
    ]
    for pat, rep in pats:
        n = len(re.findall(pat, html))
        if n != 1:
            die(f"[{label}] pattern /{pat}/ found {n}x (expected 1)")
        html = re.sub(pat, rep, html)
    return html


def verify(text: str, build: int, sha: str, label: str) -> None:
    problems = []
    if not re.search(rf"const BUILD = {build}\b", text):
        problems.append(f"const BUILD is not {build}")
    if f">build {build}</a>" not in text:
        problems.append(f"no-JS badge is not >build {build}</a>")
    if sha not in text:
        problems.append("source SHA not embedded")
    if "SRC_COMMIT = 'HEAD'" in text or 'commit/HEAD"' in text:
        problems.append("dev sentinel HEAD leaked")
    if re.search(r"const BUILD = 0\b", text) or ">build 0</a>" in text:
        problems.append("dev sentinel build 0 leaked")
    if problems:
        die(f"[{label}] " + "; ".join(problems))


def snapshot_2d(live: str) -> str:
    """Rewrite the live 2-D page's relative links for builds/bN/ depth."""
    for marker in ('href="3d/"', 'href="iss-sim/"', 'href="builds/"'):
        if live.count(marker) != 1:
            die(f"[2d snapshot] marker {marker} found {live.count(marker)}x")
    out = live.replace('href="3d/"', 'href="3d.html"')
    out = out.replace('href="iss-sim/"', 'href="../../iss-sim/"')
    out = out.replace('href="builds/"', 'href="../"')
    return out


def snapshot_3d(live: str) -> str:
    """Rewrite the live 3-D page's relative paths for builds/bN/ depth."""
    for marker in ('"../vendor/three-', 'href="../"', 'href="../builds/"'):
        if live.count(marker) != 1:
            die(f"[3d snapshot] marker {marker} found {live.count(marker)}x")
    out = live.replace('"../vendor/three-', '"../../vendor/three-')
    out = out.replace('href="../"', 'href="index.html"')
    out = out.replace('href="../builds/"', 'href="../"')
    return out


def catalog(manifest: list[dict]) -> str:
    rows = "\n".join(
        f'<tr><td><a href="b{m["build"]}/index.html">build {m["build"]}</a></td>'
        f'<td>{m["date"]}</td>'
        f'<td><a href="b{m["build"]}/index.html">2-D</a> &middot; '
        f'<a href="b{m["build"]}/3d.html">3-D</a></td>'
        f'<td><a href="{REPO_URL}/commit/{m["sha"]}"><code>{m["sha"][:10]}</code></a></td></tr>'
        for m in sorted(manifest, key=lambda m: -m["build"])
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Podium viewer — builds</title>
<style>
 body {{ background:#0a0e14; color:#c9d4e3; font:14px/1.6 "SF Mono", Consolas,
   monospace; max-width: 760px; margin: 40px auto; padding: 0 16px; }}
 a {{ color:#5ac8fa; }} h1 span {{ color:#5ac8fa; }}
 table {{ border-collapse: collapse; width: 100%; margin-top: 18px; }}
 td, th {{ border-bottom: 1px solid #1b2432; padding: 7px 10px; text-align: left; }}
 th {{ color:#6b7a8f; font-weight: 600; }}
 code {{ color:#ffd166; }}
</style></head><body>
<h1><span>Podium</span> viewer builds</h1>
<p>Every deployed build, immutable and side by side — old builds are
never deleted. <a href="../">latest 2-D</a> &middot;
<a href="../3d/">latest 3-D</a></p>
<table>
<tr><th>build</th><th>date (UTC)</th><th>pages</th><th>source commit</th></tr>
{rows}
</table>
</body></html>
"""


def main() -> None:
    dry = "--dry-run" in sys.argv
    first = "--first" in sys.argv
    sha = derive_sha()
    build = next_build(first)

    live_2d = inject(read(os.path.join(VIEWER, "index.html")), build, sha, "2d")
    live_3d = inject(read(os.path.join(VIEWER, "3d", "index.html")), build, sha, "3d")
    snap_2d = snapshot_2d(live_2d)
    snap_3d = snapshot_3d(live_3d)
    for label, text in (("2d", live_2d), ("3d", live_3d),
                        ("snap2d", snap_2d), ("snap3d", snap_3d)):
        verify(text, build, sha, label)
    print(f"build {build}  sha {sha[:10]}  — all four artifacts verified in memory")
    if dry:
        print("--dry-run: nothing written")
        return

    bdir = os.path.join(BUILDS, f"b{build}")
    if os.path.exists(bdir):
        die(f"{bdir} already exists — build-number collision")
    os.makedirs(bdir)
    manifest = json.load(open(MANIFEST)) if os.path.exists(MANIFEST) else []
    manifest.append({
        "build": build, "sha": sha,
        "date": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M"),
    })
    writes = {
        os.path.join(VIEWER, "index.html"): live_2d,
        os.path.join(VIEWER, "3d", "index.html"): live_3d,
        os.path.join(bdir, "index.html"): snap_2d,
        os.path.join(bdir, "3d.html"): snap_3d,
        os.path.join(BUILDS, "index.html"): catalog(manifest),
        MANIFEST: json.dumps(manifest, indent=1),
    }
    for path, text in writes.items():
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
    # read back and re-verify what actually landed on disk
    verify(read(os.path.join(VIEWER, "index.html")), build, sha, "written 2d")
    verify(read(os.path.join(VIEWER, "3d", "index.html")), build, sha, "written 3d")
    verify(read(os.path.join(bdir, "index.html")), build, sha, "written snap2d")
    verify(read(os.path.join(bdir, "3d.html")), build, sha, "written snap3d")
    print(f"done — build {build}: live pages injected, builds/b{build}/ written, "
          "catalog regenerated.")
    print("next: git add viewer/ && commit (the deploy commit) && push.")


if __name__ == "__main__":
    main()
