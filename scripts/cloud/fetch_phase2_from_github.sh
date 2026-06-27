#!/usr/bin/env bash
# On Mac: pull Phase 2 split tarball from GitHub and extract.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data"
ENERGY="${1:-30}"

cd "$ROOT"
git pull -q

PART_PREFIX="$DATA/phase2_E${ENERGY}_events.tar.gz.part-"
DEST="$DATA/phase2_E${ENERGY}_events.tar.gz"
parts=( "${PART_PREFIX}"* )
if [ ! -e "${parts[0]}" ]; then
  echo "Missing ${PART_PREFIX}* — run push on cloud first" >&2
  exit 1
fi
cat "${PART_PREFIX}"* > "$DEST"
mkdir -p "$DATA"
tar xzf "$DEST" -C "$DATA"
nev="$(ls "$DATA"/direction_scan_E${ENERGY}_*_nt_events.csv 2>/dev/null | wc -l | tr -d ' ')"
echo "Extracted $nev event CSVs for E=${ENERGY} GeV into $DATA"
