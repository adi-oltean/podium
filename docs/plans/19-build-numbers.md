# 19 — Build numbers + side-by-side viewer versions

GitHub issue: https://github.com/adi-oltean/podium/issues/19

## Problem

Deployed viewer pages had no identity: no way to tell which commit a
rendered page came from, and no way to keep older simulator versions
accessible after a redeploy.

## Fix (landed) — fermi tools/deploy.py discipline adapted

- Badge in both pages: `build N` linking to the exact podium commit;
  sources carry `build 0` / `HEAD` dev sentinels.
- `tools/deploy_viewer.py`: BUILD derived (max b<N> + 1, --first for
  the first), SHA derived with a clean-tree guard, marker injection
  with exactly-once counts, immutable `viewer/builds/bN/` snapshots
  (relative links rewritten for depth), regenerated catalog page +
  manifest.json, read-back re-verification, collision guard. Old
  builds never touched — deploys only add.
- Shared version-pinned vendor: `viewer/vendor/three-0.172.0/`;
  upgrades add pinned dirs, never mutate.
- Playwright suites check the badge in both pages.

## Acceptance Criteria

- [x] Both UI suites green pre-deploy (sentinels visible)
- [x] deploy --first produces b1; injected pages + snapshots verified
- [x] Live: badge links to commit; /builds/ catalog; b1 pages render
- [x] visualization.md build-discipline section

## Push/merge instructions

Two commits: sources (whose SHA ships), then the deploy commit from
running the tool; push; verify live.
