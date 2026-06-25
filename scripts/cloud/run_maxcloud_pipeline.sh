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
    n="$(pgrep -c Main 2>/dev/null || echo 0)"
    echo "Waiting for gaps... Main_procs=$n $(date +%H:%M:%S)"
    sleep 30
  done
}

wait_gaps
sleep 5
echo "=== Starting Phase 2 $(date) ==="
exec bash "$ROOT/scripts/cloud/run_maxcloud_campaign.sh" phase2
