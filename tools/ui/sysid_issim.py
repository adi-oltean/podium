"""System identification of the SpaceX ISS sim's controls and telemetry.

For each control button: click it, sample all HUD telemetry (RAW strings,
to expose sign characters/units), and report the measured deltas. Also
dumps initial raw telemetry and any interesting globals. Output feeds the
autopilot's hard-coded control mapping — no more in-flight sign guessing.

Run:  timeout 400 ./.venv/bin/python tools/ui/sysid_issim.py
"""

import json
import pathlib
import sys

from playwright.sync_api import sync_playwright

RAW = """() => {
  const g = sel => (document.querySelector(sel)?.textContent ?? "<none>");
  return {
    x: g("#x-range .distance"), y: g("#y-range .distance"), z: g("#z-range .distance"),
    range: g("#range .rate"), rate: g("#rate .rate"),
    pitch: g("#pitch .error"), pitchRate: g("#pitch .rate"),
    yaw: g("#yaw .error"), yawRate: g("#yaw .rate"),
    roll: g("#roll .error"), rollRate: g("#roll .rate"),
    precTrans: g("#precision-translation-status"),
    precRot: g("#precision-rotation-status"),
  };
}"""

CLICK = """(id) => {
  const el = document.getElementById(id);
  for (const t of ["pointerdown", "mousedown", "pointerup", "mouseup", "click"])
    el.dispatchEvent(new MouseEvent(t, { bubbles: true, cancelable: true, view: window }));
}"""

BUTTONS = [
    "pitch-up-button", "pitch-down-button",
    "yaw-left-button", "yaw-right-button",
    "roll-left-button", "roll-right-button",
    "translate-forward-button", "translate-backward-button",
    "translate-left-button", "translate-right-button",
    "translate-up-button", "translate-down-button",
]

out: dict = {}
with sync_playwright() as p:
    browser = p.chromium.launch(args=["--no-sandbox", "--disable-dev-shm-usage"])
    try:
        page = browser.new_page(viewport={"width": 1280, "height": 800})
        page.goto("https://iss-sim.spacex.com/", wait_until="load", timeout=90_000)
        page.wait_for_selector("#begin-button", state="visible", timeout=180_000)
        page.wait_for_timeout(2_000)
        page.click("#begin-button")
        page.wait_for_function(
            "() => (document.querySelector('#range .rate')?.textContent || '').trim().length > 0",
            timeout=90_000,
        )
        page.wait_for_timeout(6_000)  # intro settles

        out["initial_raw"] = page.evaluate(RAW)

        # interesting globals the sim may expose
        out["globals"] = page.evaluate(
            """() => Object.keys(window).filter(k =>
                 /camera|dragon|iss|craft|motion|rate|fixed|player/i.test(k)).slice(0, 40)"""
        )

        # per-button probes: 2 clicks, watch telemetry evolve for 3 s
        probes = {}
        for btn in BUTTONS:
            before = page.evaluate(RAW)
            page.evaluate(CLICK, btn)
            page.wait_for_timeout(150)
            page.evaluate(CLICK, btn)
            page.wait_for_timeout(1_500)
            mid = page.evaluate(RAW)
            page.wait_for_timeout(1_500)
            after = page.evaluate(RAW)
            probes[btn] = {"before": before, "t1.5": mid, "t3.0": after}
            # undo: 2 clicks of the opposite button (pair order in BUTTONS)
            idx = BUTTONS.index(btn)
            opp = BUTTONS[idx + 1] if idx % 2 == 0 else BUTTONS[idx - 1]
            page.evaluate(CLICK, opp)
            page.wait_for_timeout(150)
            page.evaluate(CLICK, opp)
            page.wait_for_timeout(800)
        out["probes"] = probes
        out["final_raw"] = page.evaluate(RAW)
    finally:
        browser.close()

path = pathlib.Path("tmp/ro/issim_sysid.json")
path.write_text(json.dumps(out, indent=1))
print(f"wrote {path}")
print("initial:", json.dumps(out["initial_raw"]))
print("globals:", out["globals"])
sys.exit(0)
