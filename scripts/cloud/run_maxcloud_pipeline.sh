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
      echo "Gaps complete $(date)"
      return 0
    fi
    # Data-complete fallback: each energy must reach 1200 real primaries (not CSV line count).
    if python3 - <<'PY' 2>/dev/null
from pathlib import Path

def count_events(path: Path) -> int:
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.strip() and not line.lstrip().startswith("#"):
                n += 1
    return n

d = Path("/root/GammaRayTelescope/build/data")
for e in (500, 1000):
    files = sorted(d.glob(f"aeff_energy_{e}GeV_shard*_run0_nt_events.csv"))
    ev = sum(count_events(f) for f in files)
    if ev < 1200:
        raise SystemExit(1)
raise SystemExit(0)
PY
    then
      echo "Gaps data complete (2400 events) $(date)"
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
echo "=== Starting Phase 2 $(date) ==="
exec bash "$ROOT/scripts/cloud/run_maxcloud_campaign.sh" phase2
