"""Playwright UI test for the 3-D viewer (viewer/3d/index.html).

Checks: zero JS errors, WebGL scene actually draws (draw calls > 0 and
non-blank pixels), playback advances the sim clock, scrubbing to the end
matches the physics (range ~20 m), camera toggle switches modes.

Run:  timeout 180 ./.venv/bin/python tools/ui/test_viewer3d.py
Hygiene: one browser + one http.server, both torn down in finally
(fermi pattern); after an interrupted run: bash tools/ram_sweep.sh
"""

import pathlib
import subprocess
import sys

from playwright.sync_api import sync_playwright

ROOT = pathlib.Path(__file__).resolve().parents[2]
PORT = 8902
SHOTS = pathlib.Path("/tmp/screenshots")
SHOTS.mkdir(exist_ok=True)

server = subprocess.Popen(
    [sys.executable, "-m", "http.server", str(PORT), "--bind", "127.0.0.1"],
    cwd=ROOT / "viewer",
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)
errors: list[str] = []
failures: list[str] = []


def check(cond: bool, label: str) -> None:
    print(("PASS " if cond else "FAIL ") + label)
    if not cond:
        failures.append(label)


try:
    with sync_playwright() as p:
        browser = p.chromium.launch(
            args=["--no-sandbox", "--disable-dev-shm-usage",
                  "--use-gl=swiftshader"]
        )
        try:
            page = browser.new_page(viewport={"width": 1400, "height": 900})
            page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
            page.on("pageerror", lambda e: errors.append("PAGEERR: " + str(e)))
            page.goto(f"http://127.0.0.1:{PORT}/3d/index.html",
                      wait_until="networkidle")
            page.wait_for_timeout(1200)

            check(page.evaluate("!!window.__P3D") is True, "debug handle present")
            badge = page.get_attribute("#buildlink", "href") or ""
            btxt = page.inner_text("#buildlink")
            check("/commit/" in badge and btxt.startswith("build "),
                  f"build badge present ({btxt} -> {badge.rsplit('/', 1)[-1][:10]})")
            check(page.evaluate("window.__P3D.drawCalls()") > 0, "webgl draw calls > 0")

            # non-blank canvas: sample pixels for variety
            distinct = page.evaluate(
                """() => {
                  const c = document.getElementById('view');
                  const g = document.createElement('canvas');
                  g.width = 64; g.height = 64;
                  const ctx = g.getContext('2d');
                  ctx.drawImage(c, 0, 0, 64, 64);
                  const d = ctx.getImageData(0, 0, 64, 64).data;
                  const s = new Set();
                  for (let i = 0; i < d.length; i += 4)
                    s.add(d[i] << 16 | d[i+1] << 8 | d[i+2]);
                  return s.size;
                }"""
            )
            check(distinct > 8, f"canvas non-blank ({distinct} distinct colors)")

            # playback advances
            page.click("#play")
            page.wait_for_timeout(1500)
            t1 = page.evaluate("window.__P3D.t()")
            check(t1 > 10.0, f"playback advances (t={t1:.1f}s)")

            # scrub to the end: physics end state (range ~ 20 m)
            page.evaluate("window.__P3D.play(false)")
            page.evaluate("window.__P3D.setT(1e9)")
            rng = page.evaluate("window.__P3D.range()")
            check(abs(rng - 20.0) < 8.0, f"end range ~20 m (got {rng:.1f})")
            hud = page.inner_text("#hud")
            check("dv spent" in hud and "3.7" in hud, "hud shows full dv budget")

            # camera toggle
            check(page.evaluate("window.__P3D.mode()") == "follow", "starts in follow")
            page.click("#mode")
            check(page.evaluate("window.__P3D.mode()") == "orbit", "toggles to orbit")

            page.screenshot(path=str(SHOTS / "viewer3d.png"))
            check(len(errors) == 0, f"zero console errors ({errors[:3]})")
        finally:
            browser.close()
finally:
    server.terminate()
    server.wait(timeout=10)

if failures:
    print("FAILURES:", failures)
    sys.exit(1)
print("ALL PASS")
