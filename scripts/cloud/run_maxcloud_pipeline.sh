#!/usr/bin/env bash
# Wait for gaps to finish, then launch Phase 2 on MaxCloudON (run on cloud VM).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data"
LOG="$DATA/maxcloud_pipeline.log"

exec >>"$LOG" 2>&1
echo "=== Pipeline watch $(date) ==="

wait_gaps() {
  while true; do
    if [ -f "$DATA/maxcloud_gaps.log" ] && grep -q "=== Done gaps" "$DATA/maxcloud_gaps.log"; then
      echo "Gaps complete (log) $(date)"
      return 0
    fi
    if [ -f "$DATA/maxcloud_gaps_e1000.log" ] && grep -q "=== Done gaps_e1000" "$DATA/maxcloud_gaps_e1000.log"; then
      echo "Gaps E1000 complete (log) $(date)"
      echo "=== Done gaps $(date) ===" >>"$DATA/maxcloud_gaps.log"
      return 0
    fi
    # Data-complete fallback: E1000 must reach 1200; E500 may be 1190+ (prior partial run).
    if python3 - <<'PY' 2>/dev/null
from pathlib import Path

def count_events(path: Path) -> int:
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.strip() and not line.lstrip().startswith("#"):
                n += 1
    return n

def energy_events(energy: int) -> int:
    d = Path("/root/GammaRayTelescope/build/data")
    files = sorted(d.glob(f"aeff_energy_{energy}GeV_shard*_run0_nt_events.csv"))
    return sum(count_events(f) for f in files)

e500, e1000 = energy_events(500), energy_events(1000)
if e1000 >= 1200 and e500 >= 1190:
    raise SystemExit(0)
raise SystemExit(1)
PY
    then
      echo "Gaps data complete (E500+E1000) $(date)"
      echo "=== Done gaps $(date) ===" >>"$DATA/maxcloud_gaps.log"
      return 0
    fi
    n="$(pgrep -c Main 2>/dev/null || echo 0)"
    echo "Waiting for gaps... Main_procs=$n $(date +%H:%M:%S)"
    sleep 30
  done
}

wait_gaps
sleep 5
pkill -9 Main 2>/dev/null || true
pkill -f run_maxcloud_campaign 2>/dev/null || true
sleep 3
echo "=== Starting Phase 2 $(date) ==="
exec bash "$ROOT/scripts/cloud/run_maxcloud_campaign.sh" phase2
