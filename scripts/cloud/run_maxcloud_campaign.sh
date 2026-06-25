#!/usr/bin/env bash
# Run cloud campaigns on MaxCloudON (serial shards, N parallel = nproc).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
MODE="${1:-gaps}"
NCPU="$(nproc)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"
LOG="$DATA/maxcloud_${MODE}.log"

source "$ROOT/setup_env.sh"
export G4MT=0 G4WRITE_PHOTONS=0 G4ACD=on G4SEED="${G4SEED:-525252}"

mkdir -p "$DATA" "$BUILD/data" "$BUILD_HAD/data"
exec > >(tee -a "$LOG") 2>&1

echo "=== MaxCloud campaign $MODE $(date) ncpu=$NCPU ==="

run_shard() {
  local tag="$1" macro="$2" hadronic="$3" seed="$4"
  if [ "$hadronic" = 1 ]; then
    cp -f "$ROOT/macros/$macro" "$BUILD_HAD/macros/" 2>/dev/null || cp -f "$ROOT/macros/gaps/$macro" "$BUILD_HAD/macros/gaps/"
    (cd "$BUILD_HAD" && G4OUTPUT_TAG="$tag" G4SEED="$seed" ./MainHad "macros/$macro")
    cp -f "$BUILD_HAD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  else
    cp -f "$ROOT/macros/$macro" "$BUILD/macros/" 2>/dev/null || cp -f "$ROOT/macros/acd/$macro" "$BUILD/macros/acd/"
    (cd "$BUILD" && G4OUTPUT_TAG="$tag" G4SEED="$seed" ./Main "macros/$macro")
    cp -f "$BUILD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  fi
}

run_pool() {
  local -a jobs=()
  local idx=0
  while IFS='|' read -r tag macro had; do
    [ -z "$tag" ] && continue
    idx=$((idx + 1))
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$NCPU" ]; do
      wait -n 2>/dev/null || sleep 1
    done
    echo ">> $tag $macro"
    run_shard "$tag" "$macro" "$had" $((G4SEED + idx)) &
  done
  wait
}

case "$MODE" in
  gaps)
    python3 "$ROOT/scripts/generate_gaps_campaign.py" --aeff-energy-events 1200 --aeff-theta-events 600 --proton-events 50
    run_pool <<EOF
aeff_energy_500GeV|gaps/aeff_E500GeV.mac|0
aeff_energy_1000GeV|gaps/aeff_E1000GeV.mac|0
EOF
    ;;
  phase2)
    python3 "$ROOT/scripts/generate_direction_scan.py" --events 500
    # Parallel shards: one macro per run file in macros/parallel if generated
    if [ -d "$ROOT/macros/parallel" ]; then
      idx=0
      for mac in "$ROOT/macros/parallel"/*.mac; do
        idx=$((idx + 1))
        base=$(basename "$mac" .mac)
        echo "aeff_phase2_${base}|parallel/${base}.mac|0"
      done | run_pool
    else
      echo "Run: python3 scripts/generate_direction_scan.py with shard generator first" >&2
      exit 1
    fi
    ;;
  acd)
    python3 "$ROOT/scripts/generate_acd_campaign.py"
    grep '|em$' "$ROOT/macros/acd/manifest.txt" | while IFS='|' read -r tag macro _ _; do
      echo "${tag}|acd/${macro}|0"
    done | run_pool
    grep '|hadronic$' "$ROOT/macros/acd/manifest.txt" | while IFS='|' read -r tag macro _ _; do
      echo "${tag}|acd/${macro}|1"
    done | run_pool
    ;;
  *)
    echo "Usage: $0 {gaps|phase2|acd}" >&2
    exit 1
    ;;
esac

echo "=== Done $MODE $(date) ==="
