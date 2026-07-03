"""Playwright UI test for the V-bar approach viewer (viewer/index.html).

Checks: page loads with zero JS errors, HUD telemetry is sane, playback
advances the sim clock, the profile scrubber seeks, the camera toggle works,
and the end state matches the physics (range ~20 m, dv ~3.76 m/s).

Run:  timeout 180 ./.venv/bin/python tools/ui/test_viewer.py
Hygiene: one browser + one http.server, both torn down in finally
(fermi pattern); after an interrupted run: bash tools/ram_sweep.sh
"""

import pathlib
import subprocess
import sys

from playwright.sync_api import sync_playwright

ROOT = pathlib.Path(__file__).resolve().parents[2]
PORT = 8901
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
            args=["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"]
        )
        try:
            page = browser.new_page(viewport={"width": 1400, "height": 900})
            page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
            page.on("pageerror", lambda e: errors.append("PAGEERR: " + str(e)))
            page.goto(f"http://127.0.0.1:{PORT}/index.html", wait_until="networkidle")
            page.wait_for_timeout(800)

            # data + boot
            meta = page.evaluate("DATA.meta")
            check(abs(meta["dv_total"] - 3.7645) < 1e-3, f"dv_total in meta = {meta['dv_total']}")
            check(page.evaluate("T.length") == 1207, "1207 samples embedded")
            scen = page.text_content("#scenario") or ""
            check("V-bar glideslope" in scen, "scenario header populated")

            # HUD at t=0
            hud = page.text_content("#hud") or ""
            check("range" in hud and "1000" in hud.replace(",", ""), "initial range ~1000 m")

            # playback advances the sim clock
            t_before = page.evaluate("simT")
            page.click("#play")
            page.wait_for_timeout(2000)
            t_after = page.evaluate("simT")
            check(t_after > t_before + 60, f"60x playback advanced simT to {t_after:.0f}s")
            page.click("#play")  # pause

            # scrubber seek to end: state matches physics
            page.evaluate("simT = TMAX; syncScrub();")
            page.wait_for_timeout(400)
            end = page.evaluate("stateAt(TMAX)")
            rng = (end[0] ** 2 + end[1] ** 2 + end[2] ** 2) ** 0.5
            check(15 < rng < 25, f"end range {rng:.1f} m (expect ~20)")
            hud_end = page.text_content("#hud") or ""
            check("3.76" in hud_end, "HUD shows full dv spent at end")

            # profile canvas seek (click at 50%)
            box = page.locator("#profile").bounding_box()
            page.mouse.click(box["x"] + box["width"] * 0.5, box["y"] + box["height"] * 0.5)
            page.wait_for_timeout(300)
            t_mid = page.evaluate("simT")
            check(abs(t_mid - 1200) < 120, f"profile click seeks to mid ({t_mid:.0f}s)")

            # camera toggle
            page.click("#mode")
            check(
                (page.text_content("#mode") or "").strip() == "camera: overview",
                "camera toggles to overview",
            )

            page.screenshot(path=str(SHOTS / "ss-podium-viewer.png"))
            print(f"screenshot: {SHOTS}/ss-podium-viewer.png")
        finally:
            browser.close()
finally:
    server.terminate()
    server.wait(timeout=10)

check(not errors, "zero JS console/page errors")
for e in errors[:10]:
    print("  JS:", e)
print("=" * 40)
print("ALL PASS" if not failures else f"{len(failures)} FAILURES")
sys.exit(0 if not failures else 1)
