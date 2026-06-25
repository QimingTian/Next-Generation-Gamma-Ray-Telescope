#!/usr/bin/env bash
# Run remaining campaigns after energy_scan_on (skip duplicate energy_scan_off if filter compare exists).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"

source "$ROOT/setup_env.sh"
source "$ROOT/.venv/bin/activate"

export G4SEED="${G4SEED:-123456789}"
export G4WRITE_PHOTONS=0

copy_macros() {
  cp -f "$ROOT/macros/"*.mac "$BUILD/macros/"
  cp -f "$ROOT/macros/"*.mac "$BUILD_HAD/macros/"
}

run() {
  local tag="$1" macro="$2" filter="${3:-on}"
  export G4OUTPUT_TAG="$tag" G4MACRO="$macro" G4FILTER="$filter"
  echo "--- $tag filter=$filter ---"
  (cd "$BUILD" && ./Main "macros/$macro")
  cp -f "$BUILD/data/${tag}_"* "$DATA/" 2>/dev/null || true
}

python "$ROOT/scripts/generate_energy_campaign.py" --cal-events 10 --res-events 40
copy_macros

# Optional full filter-off scan (skip if SKIP_ENERGY_SCAN_OFF=1 and filter compare exists)
if [[ "${SKIP_ENERGY_SCAN_OFF:-1}" == "1" ]] && [[ -f "$DATA/energy_resolution_100_off_run0_nt_events.csv" ]]; then
  echo "Skipping energy_scan_off (filter fraction from resolution campaigns)"
else
  run energy_scan_off energy_scan.mac off
fi

run peak_factor peak_factor.mac on
run plane_source plane_source.mac on

export G4WRITE_PHOTONS=1
run spectrum_photons peak_factor.mac on

bash "$ROOT/scripts/sync_campaign_data.sh"
python "$ROOT/analysis/energy_calibration.py" --data-dir "$DATA"
python "$ROOT/analysis/energy_resolution.py" --tag energy_resolution_215_on
python "$ROOT/analysis/energy_resolution_filter_compare.py" --data-dir "$DATA"
python "$ROOT/analysis/shower_length.py" --data-dir "$DATA" --tag energy_scan_on
python "$ROOT/analysis/peak_factor_stats.py" --tag peak_factor
python "$ROOT/analysis/saturation.py"
python "$ROOT/analysis/filter_folding.py" --tag spectrum_photons
python "$ROOT/analysis/effective_area.py" --tag plane_source --data-dir "$DATA"
python "$ROOT/analysis/background_rejection.py" --data-dir "$DATA"
python "$ROOT/analysis/angular_resolution.py" --tag direction_scan --train-route-b
cp -f "$ROOT/figures/"*.png "$ROOT/paper/"
echo "Remaining campaigns complete."
