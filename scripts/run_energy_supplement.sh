#!/usr/bin/env bash
# Re-run energy + downstream campaigns after filter fix (SiPMSD per-photon T(λ)).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"

source "$ROOT/setup_env.sh"
cd "$BUILD" && cmake --build . -j4
cd "$BUILD_HAD" && cmake --build . -j4

export G4SEED="${G4SEED:-123456789}"
export G4WRITE_PHOTONS=0

copy_macros() {
  cp -f "$ROOT/macros/"*.mac "$BUILD/macros/"
  cp -f "$ROOT/macros/"*.mac "$BUILD_HAD/macros/"
}

run_macro() {
  local tag="$1" macro="$2" filter="${3:-on}"
  export G4OUTPUT_TAG="$tag" G4MACRO="$macro" G4FILTER="$filter"
  echo "--- tag=$tag filter=$filter ---"
  (cd "$BUILD" && ./Main "macros/$macro")
  cp -f "$BUILD/data/${tag}_"* "$DATA/" 2>/dev/null || true
}

run_had_macro() {
  local tag="$1" macro="$2"
  export G4OUTPUT_TAG="$tag" G4MACRO="$macro" G4FILTER=on
  (cd "$BUILD_HAD" && ./MainHad "macros/$macro")
  cp -f "$BUILD_HAD/data/${tag}_"* "$DATA/" 2>/dev/null || true
}

cd "$ROOT"
source .venv/bin/activate
python scripts/generate_energy_campaign.py --cal-events 10 --res-events 40
copy_macros

run_macro energy_resolution_100_on energy_resolution_100.mac on
run_macro energy_resolution_100_off energy_resolution_100.mac off
run_macro energy_resolution_215_on energy_resolution_215.mac on
run_macro energy_resolution_215_off energy_resolution_215.mac off
run_macro energy_scan_on energy_scan.mac on
run_macro energy_scan_off energy_scan.mac off
run_macro peak_factor peak_factor.mac on
run_macro plane_source plane_source.mac on
run_had_macro proton_background proton_background.mac

export G4WRITE_PHOTONS=1
run_macro spectrum_photons peak_factor.mac on

python analysis/energy_calibration.py --data-dir "$DATA"
python analysis/energy_resolution.py --tag energy_resolution_215_on
python analysis/energy_resolution_filter_compare.py --data-dir "$DATA"
python analysis/shower_length.py --data-dir "$DATA" --tag energy_scan_on
python analysis/peak_factor_stats.py --tag peak_factor
python analysis/saturation.py
python analysis/filter_folding.py --tag spectrum_photons
python analysis/effective_area.py --tag plane_source --data-dir "$DATA"
python analysis/background_rejection.py --data-dir "$DATA"
cp -f figures/*.png paper/

echo "Supplementary energy campaign complete."
