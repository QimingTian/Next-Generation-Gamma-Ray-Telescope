#!/usr/bin/env bash
# Full paper-quality campaign with Geant4 multithreading (all CPU cores).
set -euo pipefail

if [ "${CAFFEINATED:-}" != "1" ] && command -v caffeinate >/dev/null; then
  export CAFFEINATED=1
  exec caffeinate -dimsu "$0" "$@"
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"
LOG="$DATA/paper_quality.log"

mkdir -p "$DATA"
exec > >(tee -a "$LOG") 2>&1

echo "=== Paper-quality campaign $(date) ==="

source "$ROOT/setup_env.sh"

# Use all logical cores for each Geant4 job (Main.cc reads G4NUM_THREADS).
export G4NUM_THREADS="${G4NUM_THREADS:-$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 8)}"
export G4MT=1
export G4SEED="${G4SEED:-123456789}"
export G4WRITE_PHOTONS=0

echo "G4NUM_THREADS=$G4NUM_THREADS"

copy_macros() {
  mkdir -p "$BUILD/macros" "$BUILD_HAD/macros"
  cp -f "$ROOT/macros/"*.mac "$BUILD/macros/"
  cp -f "$ROOT/macros/"*.mac "$BUILD_HAD/macros/"
}

run_macro() {
  local tag="$1" macro="$2" filter="${3:-on}"
  export G4OUTPUT_TAG="$tag" G4MACRO="$macro" G4FILTER="$filter"
  echo "--- $(date) tag=$tag macro=$macro filter=$filter threads=$G4NUM_THREADS ---"
  (cd "$BUILD" && ./Main "macros/$macro")
  cp -f "$BUILD/data/${tag}_"* "$DATA/" 2>/dev/null || true
}

run_had_macro() {
  local tag="$1" macro="$2"
  export G4OUTPUT_TAG="$tag" G4MACRO="$macro" G4FILTER=on
  echo "--- $(date) HADRONIC tag=$tag macro=$macro threads=$G4NUM_THREADS ---"
  (cd "$BUILD_HAD" && ./MainHad "macros/$macro")
  cp -f "$BUILD_HAD/data/${tag}_"* "$DATA/" 2>/dev/null || true
}

cd "$BUILD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=OFF >/dev/null && cmake --build . -j"${G4NUM_THREADS}"
mkdir -p "$BUILD_HAD"
cd "$BUILD_HAD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=ON >/dev/null && cmake --build . -j"${G4NUM_THREADS}"

cd "$ROOT"
source .venv/bin/activate
python scripts/generate_energy_campaign.py --cal-events 15 --res-events 60
copy_macros

# Skip campaigns already complete with full statistics (direction scan done).
echo "=== Energy resolution (filter on/off, 60 events) ==="
run_macro energy_resolution_100_on energy_resolution_100.mac on
run_macro energy_resolution_100_off energy_resolution_100.mac off
run_macro energy_resolution_215_on energy_resolution_215.mac on
run_macro energy_resolution_215_off energy_resolution_215.mac off

echo "=== Energy calibration (50--1000 GeV, filter on/off, 15 events/energy) ==="
run_macro energy_scan_on energy_scan.mac on
run_macro energy_scan_off energy_scan.mac off

echo "=== Peak factor @ 1 TeV (30 events) ==="
run_macro peak_factor peak_factor.mac on

echo "=== Effective area (5000 plane primaries @ 200 GeV) ==="
run_macro plane_source plane_source.mac on

echo "=== Proton background (hadronic) ==="
run_had_macro proton_background proton_background.mac

echo "=== Wavelength spectrum (photons ON, 1 run) ==="
export G4WRITE_PHOTONS=1
run_macro spectrum_photons peak_factor.mac on
export G4WRITE_PHOTONS=0

echo "=== Analysis ==="
python analysis/energy_calibration.py --data-dir "$DATA"
python analysis/energy_resolution.py --tag energy_resolution_215_on
python analysis/energy_resolution_filter_compare.py --data-dir "$DATA"
python analysis/shower_length.py --data-dir "$DATA" --tag energy_scan_on
python analysis/peak_factor_stats.py --tag peak_factor
python analysis/saturation.py
python analysis/filter_folding.py --tag spectrum_photons
python analysis/density_maps.py --tag peak_factor
python analysis/angular_resolution.py --tag direction_scan --train-route-b
python analysis/effective_area.py --tag plane_source --data-dir "$DATA"
python analysis/background_rejection.py --data-dir "$DATA"
cp -f figures/*.png paper/

echo "=== Paper-quality campaign finished $(date) ==="
