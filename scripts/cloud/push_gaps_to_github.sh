#!/usr/bin/env bash
# Run on MaxCloudON: pack gaps CSV shards and push tarball to GitHub.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
DATA="$ROOT/data"
BUILD_DATA="$ROOT/build/data"

echo "=== Pack gaps event CSVs ==="
mkdir -p "$DATA"
tar czf "$DATA/gaps_events.tar.gz" -C "$BUILD_DATA" \
  $(ls "$BUILD_DATA"/aeff_energy_500GeV_shard*_nt_events.csv 2>/dev/null | xargs -n1 basename) \
  $(ls "$BUILD_DATA"/aeff_energy_1000GeV_shard*_nt_events.csv 2>/dev/null | xargs -n1 basename)
ls -lh "$DATA/gaps_events.tar.gz"

echo "=== Git push tarball ==="
git pull -q --rebase origin "$(git rev-parse --abbrev-ref HEAD)" || git pull -q
git add data/gaps_events.tar.gz
if git diff --cached --quiet; then
  echo "No changes to gaps_events.tar.gz"
  exit 0
fi
git commit -m "Update gaps A_eff simulation tarball (E500+E1000 shards)."
git push origin HEAD
