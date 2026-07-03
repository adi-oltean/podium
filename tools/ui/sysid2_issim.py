"""ISS-sim system identification, round 2: direct function calls + globals.

Verifies: rotation button signs vs displayed rates, displayed-rate scale per
pulse, translation body-axis -> displayed x/y/z mapping and per-pulse dv,
and pulse-size globals. Everything sampled from the sim's exact internal
state (rateRotationX/Y/Z, motionVector, camera.quaternion).

Run:  timeout 300 ./.venv/bin/python tools/ui/sysid2_issim.py
"""

import json
import pathlib

from playwright.sync_api import sync_playwright

SNAP = """() => {
  const g = sel => parseFloat((document.querySelector(sel)?.textContent || "").replace(/[^\\d.+-]/g, "")) || 0;
  const q = camera.quaternion.clone().conjugate();
  const vb = motionVector.clone().applyQuaternion(q);
  return {
    disp: { x: g("#x-range .distance"), y: g("#y-range .distance"), z: g("#z-range .distance"),
            range: g("#range .rate"), rate: g("#rate .rate"),
            pitch: g("#pitch .error"), pitchRate: g("#pitch .rate"),
            yaw: g("#yaw .error"), yawRate: g("#yaw .rate"),
            roll: g("#roll .error"), rollRate: g("#roll .rate") },
    rr: [rateRotationX, rateRotationY, rateRotationZ],
    vBody: [vb.x, vb.y, vb.z],
    vWorld: [motionVector.x, motionVector.y, motionVector.z],
    pulses: [translationPulseSize, rotationPulseSize, rateSpeedSize],
  };
}"""


def probe(page, label, call, n, settle_ms, undo, out):
    page.evaluate(f"() => {{ for (let i = 0; i < {n}; i++) {call}(); }}")
    page.wait_for_timeout(settle_ms)
    after = page.evaluate(SNAP)
    out[label] = after
    if undo:
        page.evaluate(f"() => {{ for (let i = 0; i < {n}; i++) {undo}(); }}")
        page.wait_for_timeout(800)


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
        page.wait_for_timeout(6_000)

        out["initial"] = page.evaluate(SNAP)
        probe(page, "pitchUp_x1", "pitchUp", 1, 1500, "pitchDown", out)
        probe(page, "yawLeft_x1", "yawLeft", 1, 1500, "yawRight", out)
        probe(page, "rollLeft_x1", "rollLeft", 1, 1500, "rollRight", out)
        probe(page, "fwd_x10", "translateForward", 10, 2500, "translateBackward", out)
        probe(page, "right_x10", "translateRight", 10, 2500, "translateLeft", out)
        probe(page, "up_x10", "translateUp", 10, 2500, "translateDown", out)
        page.evaluate("toggleRotation")  # just reference-check it exists
        out["final"] = page.evaluate(SNAP)
    finally:
        browser.close()

path = pathlib.Path("tmp/ro/issim_sysid2.json")
path.write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
