#!/usr/bin/env bash
# ACD-on gamma/proton validation (~390 primaries, ~6–8 h on 12 cores).
# If local_gaps is running, defaults to 2 parallel shards to avoid CPU contention.
set -euo pipefail

if [ "${CAFFEINATED:-}" != "1" ] && command -v caffeinate >/dev/null; then
  export CAFFEINATED=1
  exec caffeinate -dimsu "$0" "$@"
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"
LOG="$DATA/acd_campaign.log"
NCPU="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 12)"

if pgrep -f 'run_local_gaps.sh|macros/gaps/' >/dev/null 2>&1; then
  NCPU="${G4NUM_THREADS:-2}"
  echo "local_gaps detected — limiting ACD pool to $NCPU cores"
fi

mkdir -p "$DATA"
exec > >(tee -a "$LOG") 2>&1

echo "=== ACD campaign $(date) ==="
echo "CPU pool: $NCPU"

source "$ROOT/setup_env.sh"
export G4MT=0 G4WRITE_PHOTONS=0 G4ACD=on G4SEED="${G4SEED:-525252}"

cd "$BUILD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=OFF >/dev/null && cmake --build . -j"$NCPU"
mkdir -p "$BUILD_HAD" "$BUILD/macros/acd" "$BUILD_HAD/macros/acd"
cd "$BUILD_HAD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=ON >/dev/null && cmake --build . -j"$NCPU"

cd "$ROOT"
python scripts/generate_acd_campaign.py --gamma-events 50 --gamma-tilt-events 30 --proton-events 50
cp -f "$ROOT/macros/acd/"*.mac "$BUILD/macros/acd/"
cp -f "$ROOT/macros/acd/"*.mac "$BUILD_HAD/macros/acd/"

run_pool() {
  local mode="$1"
  local manifest="$2"
  local idx=0 tag macro filter

  while IFS='|' read -r tag macro filter _mode; do
    [ -z "$tag" ] && continue
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$NCPU" ]; do
      wait -n 2>/dev/null || sleep 1
    done
    idx=$((idx + 1))
    shard_tag="${tag}_run${idx}"
    seed=$((G4SEED + idx))
    echo ">> launch $shard_tag ($macro) [$mode] $(date +%H:%M:%S)"
    (
      export G4WRITE_PHOTONS=0 G4ACD=on
      if [ "$mode" = hadronic ]; then
        source "$ROOT/setup_env.sh"
        export G4MT=0 G4OUTPUT_TAG="$shard_tag" G4MACRO="acd/$macro" G4FILTER="$filter" G4SEED="$seed"
        (cd "$BUILD_HAD" && ./MainHad "macros/acd/$macro")
        cp -f "$BUILD_HAD/data/${shard_tag}_"* "$DATA/" 2>/dev/null || true
      else
        source "$ROOT/setup_env.sh"
        export G4MT=0 G4OUTPUT_TAG="$shard_tag" G4MACRO="acd/$macro" G4FILTER="$filter" G4SEED="$seed"
        (cd "$BUILD" && ./Main "macros/acd/$macro")
        cp -f "$BUILD/data/${shard_tag}_"* "$DATA/" 2>/dev/null || true
      fi
    ) &
  done < "$manifest"
  wait
}

MANIFEST="$ROOT/macros/acd/manifest.txt"
grep '|em$' "$MANIFEST" > "$DATA/_pool_acd_em.txt"
grep '|hadronic$' "$MANIFEST" > "$DATA/_pool_acd_had.txt"

echo "=== Gamma (EM, ACD on) ==="
run_pool em "$DATA/_pool_acd_em.txt"

echo "=== Proton (hadronic, ACD on) ==="
run_pool hadronic "$DATA/_pool_acd_had.txt"

echo "=== Analysis ==="
source "$ROOT/.venv/bin/activate" 2>/dev/null || true
python analysis/acd_veto_scan.py --data-dir "$DATA"
python analysis/background_rejection.py \
  --gamma-tag acd_gamma_100 \
  --proton-tag acd_proton_100 \
  --data-dir "$DATA" \
  --combined

echo "=== Finished $(date) ==="
