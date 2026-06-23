#!/usr/bin/env bash
# Run Geant4 simulation campaigns and copy results to data/
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
DATA="$ROOT/data"

source "$ROOT/setup_env.sh"
cd "$BUILD"
cmake .. >/dev/null
cmake --build . -j4

mkdir -p "$DATA"
export G4SEED="${G4SEED:-123456789}"

run_macro() {
  local tag="$1"
  local macro="$2"
  export G4OUTPUT_TAG="$tag"
  export G4MACRO="$macro"
  echo "=== Running $macro (tag=$tag) ==="
  ./Main "macros/$macro"
  cp -f data/"${tag}"_* "$DATA/" 2>/dev/null || true
}

run_macro energy_scan energy_scan.mac
run_macro energy_resolution energy_resolution.mac
run_macro angle_scan angle_scan.mac
run_macro peak_factor peak_factor.mac

echo "=== Running analysis ==="
cd "$ROOT"
source .venv/bin/activate 2>/dev/null || python3 -m venv .venv && source .venv/bin/activate
pip install -q -r requirements.txt

python analysis/Filter-Design.py
python analysis/density_maps.py --tag peak_factor
python analysis/energy_calibration.py --data-dir "$DATA"
python analysis/energy_resolution.py --tag energy_resolution
python analysis/shower_length.py --data-dir "$DATA"
python analysis/filter_folding.py --tag peak_factor
python analysis/saturation.py

echo "Campaign complete. Figures in figures/, summaries in analysis/output/"
