#!/usr/bin/env bash
# Overnight campaign: direction recon + energy + A_eff + background
set -euo pipefail

if [ "${CAFFEINATED:-}" != "1" ] && command -v caffeinate >/dev/null; then
  export CAFFEINATED=1
  exec caffeinate -dimsu "$0" "$@"
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"
LOG="$ROOT/data/overnight.log"

mkdir -p "$DATA"
exec > >(tee -a "$LOG") 2>&1

echo "=== Overnight campaign started $(date) ==="

source "$ROOT/setup_env.sh"

export G4SEED="${G4SEED:-123456789}"
export G4WRITE_PHOTONS=0

copy_macros() {
  mkdir -p "$BUILD/macros" "$BUILD_HAD/macros"
  cp -f "$ROOT/macros/"*.mac "$BUILD/macros/"
  cp -f "$ROOT/macros/"*.mac "$BUILD_HAD/macros/"
}

run_macro() {
  local tag="$1"
  local macro="$2"
  local filter="${3:-on}"
  export G4OUTPUT_TAG="$tag"
  export G4MACRO="$macro"
  export G4FILTER="$filter"
  echo "--- $(date) tag=$tag macro=$macro filter=$filter ---"
  (cd "$BUILD" && ./Main "macros/$macro")
  cp -f "$BUILD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  if ! ls "$DATA/${tag}_run0_nt_events.csv" >/dev/null 2>&1; then
    echo "ERROR: no output for tag=$tag" >&2
    return 1
  fi
}

run_had_macro() {
  local tag="$1"
  local macro="$2"
  export G4OUTPUT_TAG="$tag"
  export G4MACRO="$macro"
  export G4FILTER=on
  echo "--- $(date) HADRONIC tag=$tag macro=$macro ---"
  (cd "$BUILD_HAD" && ./MainHad "macros/$macro")
  cp -f "$BUILD_HAD/data/${tag}_"* "$DATA/" 2>/dev/null || true
}

# Build standard binary
cd "$BUILD"
cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=OFF >/dev/null
cmake --build . -j4

# Build hadronic binary
mkdir -p "$BUILD_HAD"
cd "$BUILD_HAD"
cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=ON >/dev/null
cmake --build . -j4

# Python env + filter table
cd "$ROOT"
python3 -m venv .venv 2>/dev/null || true
source .venv/bin/activate
pip install -q -r requirements.txt scikit-learn
python analysis/Filter-Design.py
python scripts/generate_direction_scan.py --events 12
python scripts/generate_energy_campaign.py --cal-events 10 --res-events 40
copy_macros

echo "=== Priority 1: direction scan ==="
run_macro direction_scan direction_scan.mac on

echo "=== Priority 2: energy resolution (filter on/off) ==="
run_macro energy_resolution_100_on energy_resolution_100.mac on
run_macro energy_resolution_100_off energy_resolution_100.mac off
run_macro energy_resolution_215_on energy_resolution_215.mac on
run_macro energy_resolution_215_off energy_resolution_215.mac off

echo "=== Priority 3: energy calibration ==="
run_macro energy_scan_on energy_scan.mac on
run_macro energy_scan_off energy_scan.mac off

echo "=== Priority 4: peak factor ==="
run_macro peak_factor peak_factor.mac on

echo "=== Priority 5: effective area ==="
run_macro plane_source plane_source.mac on

echo "=== Priority 6: proton background ==="
run_had_macro proton_background proton_background.mac

echo "=== Priority 7: one photons-on run for spectrum ==="
export G4WRITE_PHOTONS=1
run_macro spectrum_photons peak_factor.mac on
export G4WRITE_PHOTONS=0

echo "=== Analysis ==="
python analysis/energy_calibration.py --data-dir "$DATA"
python analysis/energy_resolution.py --tag energy_resolution_100_on
python analysis/energy_resolution.py --tag energy_resolution_215_on
python analysis/energy_resolution_filter_compare.py --data-dir "$DATA"
python analysis/shower_length.py --data-dir "$DATA"
python analysis/saturation.py
python analysis/filter_folding.py --tag spectrum_photons
python analysis/density_maps.py --tag peak_factor
python analysis/angular_resolution.py --tag direction_scan --train-route-b
python analysis/effective_area.py --tag plane_source
python analysis/background_rejection.py --gamma-tag direction_scan --proton-tag proton_background

# Copy figures to paper/
cp -f figures/*.png paper/ 2>/dev/null || true

echo "=== Overnight campaign finished $(date) ==="
