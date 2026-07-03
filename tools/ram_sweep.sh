#!/usr/bin/env bash
# WSL RAM hygiene (pattern from ../fermi — this box has ~8 GB and has
# OOM-crashed from accumulated tooling). Kills ONLY the orphaned helper
# processes this project spawns — Playwright browsers and test http.servers.
# It does NOT touch Cursor / Electron / system processes and does NOT delete
# any disk cache.
#
# Run after any Playwright run that was interrupted or timed out (a `timeout`
# SIGTERM skips `finally:` cleanup and orphans the browser), or whenever RAM
# feels tight:   bash tools/ram_sweep.sh

sweep() {   # $1 = pgrep -f regex, $2 = human label
  local pids
  pids=$(pgrep -f "$1" || true)
  if [ -n "$pids" ]; then
    echo "  $2: killing $(echo "$pids" | tr '\n' ' ')"
    kill $pids 2>/dev/null || true
    sleep 1
    pids=$(pgrep -f "$1" || true)
    [ -n "$pids" ] && kill -9 $pids 2>/dev/null || true
  else
    echo "  $2: none"
  fi
}

echo "RAM sweep — orphaned project tools only (no disk touched):"
sweep 'ms-playwright'   'Playwright browsers'
sweep 'http\.server'    'test http.server'
echo "--"
free -h | awk 'NR<=2'
