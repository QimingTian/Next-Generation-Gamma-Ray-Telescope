#!/usr/bin/env bash
# Copy completed Geant4 outputs from build dirs into data/.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/data"
for dir in "$ROOT/build/data" "$ROOT/build-hadronic/data"; do
  [ -d "$dir" ] || continue
  cp -f "$dir"/*_run* "$DATA/" 2>/dev/null || true
done
echo "Synced campaign data to $DATA"
