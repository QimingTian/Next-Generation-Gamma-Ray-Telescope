#!/usr/bin/env bash
set -euo pipefail
ENERGY="${1:-30}"
ROOT=/root/GammaRayTelescope
cd "$ROOT"
git pull -q || true

BD="$ROOT/build/data"
OUT="$ROOT/data/phase2_E${ENERGY}_events.tar.gz"
PART="$ROOT/data/phase2_E${ENERGY}_events.tar.gz.part-"
mkdir -p "$ROOT/data"
cd "$BD"
shopt -s nullglob
files=(direction_scan_E${ENERGY}_*_nt_*.csv)
echo "files: ${#files[@]}"
tar czf "$OUT" "${files[@]}"
ls -lh "$OUT"
rm -f "${PART}"*
split -b 95m "$OUT" "$PART"
ls -lh "${PART}"*

git add "$OUT" "${PART}"*
git commit -m "Add Phase 2 E${ENERGY} GeV direction_scan data (split tarball)." || true
git push origin HEAD
echo DONE
