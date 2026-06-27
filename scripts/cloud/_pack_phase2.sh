#!/usr/bin/env bash
set -euo pipefail
ENERGY="${1:-}"
ROOT=/root/GammaRayTelescope
BD="$ROOT/build/data"
OUT="$ROOT/data/phase2_events.tar.gz"
mkdir -p "$ROOT/data"
cd "$BD"
shopt -s nullglob
if [ -n "$ENERGY" ]; then
  files=(direction_scan_E${ENERGY}_*_nt_*.csv)
else
  files=(direction_scan_*_nt_*.csv)
fi
if ((${#files[@]}==0)); then
  echo "No CSV files in $BD" >&2
  ls "$BD"/direction_scan_* 2>/dev/null | head -5 || true
  exit 1
fi
tar czf "$OUT" "${files[@]}"
ls -lh "$OUT"
echo "packed ${#files[@]} files"
