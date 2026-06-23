#!/usr/bin/env bash
# Finish remaining campaigns after energy_scan (plane source, sync, analysis).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
DATA="$ROOT/data"
source "$ROOT/setup_env.sh"
source "$ROOT/.venv/bin/activate"

copy_macros() { cp -f "$ROOT/macros/"*.mac "$BUILD/macros/"; }

run_plane() {
  export G4SEED="${G4SEED:-123456790}" G4WRITE_PHOTONS=0 G4OUTPUT_TAG=plane_source
  export G4MACRO=plane_source.mac G4FILTER=on
  copy_macros
  echo "--- plane_source (1000 primaries) ---"
  (cd "$BUILD" && ./Main macros/plane_source.mac)
  cp -f "$BUILD/data/plane_source_"* "$DATA/" 2>/dev/null || true
}

sync_and_analyze() {
  bash "$ROOT/scripts/sync_campaign_data.sh"
  python "$ROOT/analysis/energy_calibration.py" --data-dir "$DATA"
  python "$ROOT/analysis/energy_resolution.py" --tag energy_resolution_215_on
  python "$ROOT/analysis/energy_resolution_filter_compare.py" --data-dir "$DATA"
  python "$ROOT/analysis/shower_length.py" --data-dir "$DATA" --tag energy_scan_on
  python "$ROOT/analysis/peak_factor_stats.py" --tag peak_factor 2>/dev/null || true
  python "$ROOT/analysis/saturation.py"
  python "$ROOT/analysis/filter_folding.py" --tag spectrum_photons 2>/dev/null || \
    python "$ROOT/analysis/filter_folding.py" --tag peak_factor 2>/dev/null || true
  python "$ROOT/analysis/effective_area.py" --tag plane_source --data-dir "$DATA"
  python "$ROOT/analysis/background_rejection.py" --data-dir "$DATA"
  python "$ROOT/analysis/angular_resolution.py" --tag direction_scan --train-route-b
  cp -f "$ROOT/figures/"*.png "$ROOT/paper/"
}

case "${1:-all}" in
  plane) run_plane ;;
  analyze) sync_and_analyze ;;
  all) run_plane && sync_and_analyze ;;
  *) echo "Usage: $0 [plane|analyze|all]" ;;
esac
