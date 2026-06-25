#!/usr/bin/env bash
# Paper-gap supplements sized for ~12 h on 12-core Mac (~6750 primaries).
set -euo pipefail

if [ "${CAFFEINATED:-}" != "1" ] && command -v caffeinate >/dev/null; then
  export CAFFEINATED=1
  exec caffeinate -dimsu "$0" "$@"
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"
LOG="$DATA/local_gaps.log"
NCPU="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 12)"

mkdir -p "$DATA"
exec > >(tee -a "$LOG") 2>&1

echo "=== Local gaps campaign $(date) ==="
echo "CPU pool: $NCPU"

source "$ROOT/setup_env.sh"
export G4MT=0 G4WRITE_PHOTONS=0 G4SEED="${G4SEED:-424242}"

cd "$BUILD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=OFF >/dev/null && cmake --build . -j"$NCPU"
mkdir -p "$BUILD_HAD" "$BUILD/macros/gaps" "$BUILD_HAD/macros/gaps"
cd "$BUILD_HAD" && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=ON >/dev/null && cmake --build . -j"$NCPU"

cd "$ROOT"
python scripts/generate_gaps_campaign.py --aeff-energy-events 1200 --aeff-theta-events 600 --proton-events 50
cp -f "$ROOT/macros/gaps/"*.mac "$BUILD/macros/gaps/"
cp -f "$ROOT/macros/gaps/"*.mac "$BUILD_HAD/macros/gaps/"
chmod +x "$ROOT/scripts/run_geant4_shard.sh"

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
    if compgen -G "$DATA/${shard_tag}_"*_nt_events.csv >/dev/null 2>&1; then
      echo ">> skip $shard_tag (already complete) $(date +%H:%M:%S)"
      continue
    fi
    echo ">> launch $shard_tag ($macro) [$mode] $(date +%H:%M:%S)"
    (
      export G4WRITE_PHOTONS=0
      if [ "$mode" = hadronic ]; then
        source "$ROOT/setup_env.sh"
        export G4MT=0 G4OUTPUT_TAG="$shard_tag" G4MACRO="gaps/$macro" G4FILTER="$filter" G4SEED="$seed"
        (cd "$BUILD_HAD" && ./MainHad "macros/gaps/$macro")
        cp -f "$BUILD_HAD/data/${shard_tag}_"* "$DATA/" 2>/dev/null || true
      else
        bash "$ROOT/scripts/run_geant4_shard.sh" "$shard_tag" "gaps/$macro" "$filter" "$seed"
      fi
    ) &
  done < "$manifest"
  wait
}

MANIFEST="$ROOT/macros/gaps/manifest.txt"
grep '|em$' "$MANIFEST" > "$DATA/_pool_gaps_em.txt"
grep '|hadronic$' "$MANIFEST" > "$DATA/_pool_gaps_had.txt"

echo "=== A_eff energy + theta (EM) ==="
run_pool em "$DATA/_pool_gaps_em.txt"

echo "=== Proton supplement (hadronic) ==="
run_pool hadronic "$DATA/_pool_gaps_had.txt"

echo "=== Analysis ==="
source "$ROOT/.venv/bin/activate" 2>/dev/null || true
python analysis/effective_area_scan.py --data-dir "$DATA"
python analysis/background_rejection.py --data-dir "$DATA"
cp -f figures/effective_area_scan.png figures/background_rejection.png paper/ 2>/dev/null || true

echo "=== Finished $(date) ==="
