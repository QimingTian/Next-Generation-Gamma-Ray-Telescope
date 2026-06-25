#!/usr/bin/env bash
# Paper-quality campaign: parallel serial Geant4 jobs across all CPU cores.
set -euo pipefail

if [ "${CAFFEINATED:-}" != "1" ] && command -v caffeinate >/dev/null; then
  export CAFFEINATED=1
  exec caffeinate -dimsu "$0" "$@"
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"
LOG="$DATA/paper_quality_parallel.log"
NCPU="$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 8)"

mkdir -p "$DATA"
exec > >(tee -a "$LOG") 2>&1

echo "=== Paper-quality PARALLEL campaign $(date) ==="
echo "CPU cores for job pool: $NCPU (expect ~${NCPU}00% total CPU while pool is full)"

source "$ROOT/setup_env.sh"
export G4MT=0
export G4SEED="${G4SEED:-123456789}"

cd "$BUILD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=OFF >/dev/null && cmake --build . -j"$NCPU"
mkdir -p "$BUILD_HAD" "$BUILD/macros/parallel" "$BUILD_HAD/macros/parallel"
cd "$BUILD_HAD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=ON >/dev/null && cmake --build . -j"$NCPU"

cd "$ROOT"
source .venv/bin/activate
python scripts/generate_energy_campaign.py --cal-events 15 --res-events 60
python scripts/generate_parallel_macros.py --cal-events 15 --peak-events 30 --plane-events 5000 --proton-events 50
cp -f "$ROOT/macros/"*.mac "$BUILD/macros/"
cp -f "$ROOT/macros/parallel/"*.mac "$BUILD/macros/parallel/"
cp -f "$ROOT/macros/parallel/"*.mac "$BUILD_HAD/macros/parallel/"
chmod +x "$ROOT/scripts/run_geant4_shard.sh"

run_pool() {
  local mode="$1"  # em | hadronic
  local manifest="$2"
  local write_photons="${3:-0}"
  local idx=0 tag macro filter shard_tag seed

  while IFS='|' read -r tag macro filter; do
    [ -z "$tag" ] && continue
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$NCPU" ]; do
      wait -n 2>/dev/null || sleep 1
    done
    idx=$((idx + 1))
    shard_tag="${tag}_run${idx}"
    seed=$((G4SEED + idx))
    echo ">> launch $shard_tag ($macro) [$mode] $(date +%H:%M:%S)"
    (
      export G4WRITE_PHOTONS="$write_photons"
      if [ "$mode" = hadronic ]; then
        source "$ROOT/setup_env.sh"
        export G4MT=0 G4OUTPUT_TAG="$shard_tag" G4MACRO="$macro" G4FILTER="$filter" G4SEED="$seed"
        export G4WRITE_PHOTONS="$write_photons"
        (cd "$BUILD_HAD" && ./MainHad "macros/parallel/$macro")
        cp -f "$BUILD_HAD/data/${shard_tag}_"* "$DATA/" 2>/dev/null || true
      else
        bash "$ROOT/scripts/run_geant4_shard.sh" "$shard_tag" "$macro" "$filter" "$seed"
      fi
    ) &
  done < "$manifest"
  wait
}

MANIFEST="$ROOT/macros/parallel/manifest.txt"

echo "=== Energy resolution (skip if complete) ==="
for spec in "energy_resolution_100_on:energy_resolution_100.mac:on" \
            "energy_resolution_100_off:energy_resolution_100.mac:off" \
            "energy_resolution_215_on:energy_resolution_215.mac:on" \
            "energy_resolution_215_off:energy_resolution_215.mac:off"; do
  tag="${spec%%:*}"; rest="${spec#*:}"; macro="${rest%%:*}"; filt="${rest##*:}"
  if [ ! -f "$DATA/${tag}_run0_nt_events.csv" ]; then
    cp -f "$ROOT/macros/$macro" "$BUILD/macros/"
    bash "$ROOT/scripts/run_geant4_shard.sh" "$tag" "$macro" "$filt" "$G4SEED"
  else
    echo "skip $tag"
  fi
done

echo "=== Parallel energy calibration (12 jobs) ==="
grep '^energy_scan' "$MANIFEST" > "$DATA/_pool_energy.txt"
run_pool em "$DATA/_pool_energy.txt"

echo "=== Parallel peak factor (6 shards @ 1 TeV, 30 events total) ==="
grep '^peak_factor' "$MANIFEST" > "$DATA/_pool_peak.txt"
run_pool em "$DATA/_pool_peak.txt"

echo "=== Parallel plane source (10 shards, 5000 primaries total) ==="
grep '^plane_source' "$MANIFEST" > "$DATA/_pool_plane.txt"
run_pool em "$DATA/_pool_plane.txt"

echo "=== Parallel proton background (3 energies, hadronic) ==="
grep '^proton_background' "$MANIFEST" > "$DATA/_pool_proton.txt"
run_pool hadronic "$DATA/_pool_proton.txt"

echo "=== Spectrum photons (1 TeV, photons ON) ==="
export G4WRITE_PHOTONS=1
bash "$ROOT/scripts/run_geant4_shard.sh" spectrum_photons spectrum_photons.mac on "$((G4SEED + 999))"
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

echo "=== Finished $(date) ==="
