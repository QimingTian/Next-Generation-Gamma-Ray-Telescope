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
    # Data-complete fallback: 2400 events across E500+E1000 shards
    if python3 - <<'PY' 2>/dev/null
from pathlib import Path
d = Path("/root/GammaRayTelescope/build/data")
ev = 0
for e in (500, 1000):
    files = list(d.glob(f"aeff_energy_{e}GeV_shard*_run0_nt_events.csv"))
    for f in files:
        try:
            ev += max(0, sum(1 for _ in open(f)) - 1)
        except OSError:
            pass
if ev >= 2400:
    raise SystemExit(0)
raise SystemExit(1)
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
