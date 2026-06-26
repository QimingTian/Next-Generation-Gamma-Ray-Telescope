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
# Serial shards: thread env from setup_env.sh triggers optical segfault on Linux.
unset G4FORCENUMBEROFTHREADS G4NUM_THREADS

mkdir -p "$DATA" "$BUILD/data" "$BUILD_HAD/data"
exec > >(tee -a "$LOG") 2>&1

echo "=== MaxCloud campaign $MODE $(date) ncpu=$NCPU ==="

run_shard() {
  local tag="$1" macro="$2" hadronic="$3" seed="$4"
  local src="$ROOT/macros/$macro"
  if [ "$hadronic" = 1 ]; then
    mkdir -p "$BUILD_HAD/macros/$(dirname "$macro")"
    cp -f "$src" "$BUILD_HAD/macros/$macro"
    (cd "$BUILD_HAD" && G4OUTPUT_TAG="$tag" G4SEED="$seed" ./MainHad "macros/$macro")
    cp -f "$BUILD_HAD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  else
    mkdir -p "$BUILD/macros/$(dirname "$macro")"
    cp -f "$src" "$BUILD/macros/$macro"
    (cd "$BUILD" && G4OUTPUT_TAG="$tag" G4SEED="$seed" ./Main "macros/$macro")
    cp -f "$BUILD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  fi
}

count_csv_events() {
  python3 - "$1" <<'PY'
import sys
from pathlib import Path
p = Path(sys.argv[1])
if not p.is_file():
    print(0)
    raise SystemExit(0)
n = sum(1 for ln in p.read_text().splitlines() if ln.strip() and not ln.lstrip().startswith("#"))
print(n)
PY
}

shard_complete() {
  local tag="$1" need="$2" f n
  for f in "$BUILD/data/${tag}_"*"_nt_events.csv" "$DATA/${tag}_"*"_nt_events.csv"; do
    [ -f "$f" ] || continue
    n="$(count_csv_events "$f")"
    if [ "$n" -ge "$need" ]; then
      return 0
    fi
  done
  return 1
}

expected_beam() {
  grep '/run/beamOn' "$ROOT/macros/$1" | awk '{print $2}'
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

run_pool_skip_complete() {
  local idx=0 skipped=0
  while IFS='|' read -r tag macro had; do
    [ -z "$tag" ] && continue
    local need
    need="$(expected_beam "$macro")"
    if shard_complete "$tag" "$need"; then
      echo ">> skip $tag ($need events)"
      skipped=$((skipped + 1))
      continue
    fi
    idx=$((idx + 1))
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$NCPU" ]; do
      wait -n 2>/dev/null || sleep 1
    done
    echo ">> $tag $macro"
    run_shard "$tag" "$macro" "$had" $((G4SEED + idx)) &
  done
  wait
  echo "Skipped $skipped complete shards"
}

case "$MODE" in
  gaps|gaps_e1000)
    export G4ACD=off  # ACD geometry segfaults on Linux; not needed for A_eff(E)
    if [ "$MODE" = "gaps_e1000" ]; then
      python3 "$ROOT/scripts/generate_gaps_campaign.py" \
        --cloud-shard --shard-count "$NCPU" --aeff-energy-events 1200 \
        --cloud-energies 1000
    else
      python3 "$ROOT/scripts/generate_gaps_campaign.py" \
        --cloud-shard --shard-count "$NCPU" --aeff-energy-events 1200
    fi
    grep '|em$' "$ROOT/macros/gaps/cloud_manifest.txt" | while IFS='|' read -r tag macro _ _; do
      echo "${tag}|gaps/${macro}|0"
    done | run_pool
    ;;
  phase2)
    export G4ACD=off
    export G4WRITE_PHOTONS=0
    python3 "$ROOT/scripts/generate_direction_scan.py" \
      --cloud-shard --events-per-bin 500 --subshards 4 --n-phi 24
    grep '|em$' "$ROOT/macros/phase2/cloud_manifest.txt" | while IFS='|' read -r tag macro _ _; do
      echo "${tag}|${macro}|0"
    done | run_pool_skip_complete
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
    echo "Usage: $0 {gaps|gaps_e1000|phase2|acd}" >&2
    exit 1
    ;;
esac

echo "=== Done $MODE $(date) ==="
