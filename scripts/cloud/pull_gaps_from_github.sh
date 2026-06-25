#!/usr/bin/env bash
# On Mac: pull gaps tarball from GitHub and extract for analysis.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data"

cd "$ROOT"
git pull -q

if [ ! -f "$DATA/gaps_events.tar.gz" ]; then
  echo "Missing $DATA/gaps_events.tar.gz — run push on cloud or fetch_gaps_data.sh" >&2
  exit 1
fi

mkdir -p "$DATA"
tar xzf "$DATA/gaps_events.tar.gz" -C "$DATA"
n="$(ls "$DATA"/aeff_energy_*GeV_shard*_nt_events.csv 2>/dev/null | wc -l | tr -d ' ')"
echo "Extracted $n shard CSVs into $DATA"

echo "=== Effective area analysis ==="
python3 "$ROOT/analysis/effective_area_scan.py" --data-dir "$DATA"
