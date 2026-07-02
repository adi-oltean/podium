# 02 — v0 web simulation viewer on GitHub Pages

GitHub issue: https://github.com/adi-oltean/podium/issues/2

## Problem

Podium has no visual output. docs/visualization.md designs a fermi-style
zero-build viewer; v0.1 needs at least a simple simulation web page with a
public URL now that the repo is public.

## Affected Components

| Component | Status | Notes |
|-----------|--------|-------|
| `viewer/index.html` | MISS | self-contained page, embedded sim data |
| `tmp/rw/gen_viewer_data.py` | MISS | generates the embedded scenario JSON |
| GitHub Pages | MISS | serve repo root from main |
| `docs/roadmap.md`, README | PARTIAL | link the live page |

## Fix

- Scenario: 1 km V-bar glideslope approach to a 10 m hold (the repo's
  canonical example), impulses applied between segments propagated with the
  **nonlinear truth model** (podium.dynamics.nonlinear), 2 s cadence.
- `viewer/index.html`: single file, no external requests (fermi discipline).
  Canvas 2-D rendering: V-bar plane view (along-track vs radial) +
  cross-track strip, target with keep-out circle and docking axis, fading
  trail, burn glyphs, HUD (t, range, range rate, cumulative delta-v),
  play/pause + speed + scrubber. Render loop (requestAnimationFrame)
  decoupled from the playback clock (fermi pattern).
- Pages: Actions-based deploy (`.github/workflows/pages.yml`) serving the
  `viewer/` directory as the site root. URL: https://adi-oltean.github.io/podium/
  (Note: the legacy Jekyll build path wedged — the Pages site was first
  created while the repo was private, then GitHub Pages had a
  degraded-performance incident; deleted/recreated the site with
  build_type=workflow.)

## Tests

Page is static; verification is manual load + `curl` 200 check in CI is
overkill for v0. Data generator is exercised when regenerating.

## Acceptance Criteria

- [x] Page renders and plays back the approach (script smoke-tested headless)
- [x] Pages URL live (200 OK, 2026-07-02)
- [x] Roadmap/README link it

## Push/merge instructions

Single commit on main: `02 — v0 web simulation viewer (#2)`, push, enable
Pages via API, verify URL serves 200, close #2.

## Verification steps

Open the Pages URL; confirm playback, scrubber, HUD values (total delta-v
should read ~3.7-3.8 m/s, final range ~10 m).
