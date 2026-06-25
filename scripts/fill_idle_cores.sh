#!/usr/bin/env bash
# Launch queued shards when cores are idle (skip if output already complete).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BUILD_HAD="$ROOT/build-hadronic"
DATA="$ROOT/data"
NCPU="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 12)"
MAX_NEW="${1:-6}"

source "$ROOT/setup_env.sh"
export G4MT=0 G4WRITE_PHOTONS=0 G4FILTER=on G4ACD=on G4SEED="${G4SEED:-525252}"

mkdir -p "$BUILD/macros/acd" "$BUILD_HAD/macros/acd" "$BUILD_HAD/macros/gaps" "$DATA"

running_jobs() { ps aux | grep -E '[./]Main( Had)? ' | grep -v grep | wc -l | tr -d ' '; }

shard_done() {
  local tag="$1" need="$2"
  local f
  for f in "$DATA/${tag}_run"*"_nt_events.csv" "$BUILD/data/${tag}_run"*"_nt_events.csv" \
           "$BUILD_HAD/data/${tag}_run"*"_nt_events.csv"; do
    [ -f "$f" ] || continue
    local n
    n=$(grep -c '^[0-9]' "$f" 2>/dev/null || echo 0)
    if [ "$n" -ge "$need" ]; then
      return 0
    fi
  done
  return 1
}

shard_running() {
  local macro="$1"
  ps aux | grep -F "$macro" | grep -E '[./]Main' | grep -qv grep
}

launch_em() {
  local tag="$1" macro="$2" seed="$3"
  shard_running "$macro" && return 0
  shard_done "$tag" 1 && echo "skip (done): $tag" && return 0
  while [ "$(running_jobs)" -ge "$NCPU" ]; do sleep 2; done
  cp -f "$ROOT/macros/acd/$macro" "$BUILD/macros/acd/"
  echo ">> EM $tag ($macro) seed=$seed"
  (
    export G4OUTPUT_TAG="$tag" G4SEED="$seed"
    cd "$BUILD" && ./Main "macros/acd/$macro"
    cp -f "$BUILD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  ) &
}

launch_had_acd() {
  local tag="$1" macro="$2" seed="$3"
  shard_running "$macro" && return 0
  shard_done "$tag" 1 && echo "skip (done): $tag" && return 0
  while [ "$(running_jobs)" -ge "$NCPU" ]; do sleep 2; done
  cp -f "$ROOT/macros/acd/$macro" "$BUILD_HAD/macros/acd/"
  echo ">> HAD $tag ($macro) seed=$seed"
  (
    export G4OUTPUT_TAG="$tag" G4SEED="$seed"
    cd "$BUILD_HAD" && ./MainHad "macros/acd/$macro"
    cp -f "$BUILD_HAD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  ) &
}

launch_had_gaps() {
  local tag="$1" macro="$2" seed="$3"
  shard_running "$macro" && return 0
  shard_done "$tag" 1 && echo "skip (done): $tag" && return 0
  while [ "$(running_jobs)" -ge "$NCPU" ]; do sleep 2; done
  cp -f "$ROOT/macros/gaps/$macro" "$BUILD_HAD/macros/gaps/"
  echo ">> HAD $tag ($macro) seed=$seed"
  (
    export G4OUTPUT_TAG="$tag" G4SEED="$seed"
    cd "$BUILD_HAD" && ./MainHad "macros/gaps/$macro"
    cp -f "$BUILD_HAD/data/${tag}_"* "$DATA/" 2>/dev/null || true
  ) &
}

echo "=== fill_idle_cores $(date) — running $(running_jobs)/$NCPU ==="

queued=0
# ACD gamma (EM) — only if not already running via acd_campaign
launch_em acd_gamma_100_theta45_run6 gamma_100GeV_theta45deg.mac $((G4SEED + 106)) && queued=$((queued + 1))

# ACD protons
launch_had_acd acd_proton_100GeV_run7 proton_E100GeV.mac $((G4SEED + 107)) && queued=$((queued + 1))
launch_had_acd acd_proton_300GeV_run8 proton_E300GeV.mac $((G4SEED + 108)) && queued=$((queued + 1))
launch_had_acd acd_proton_1000GeV_run9 proton_E1000GeV.mac $((G4SEED + 109)) && queued=$((queued + 1))

# Gaps proton supplement (early — run_local_gaps will skip if done)
launch_had_gaps proton_supplement_run1 proton_supp_100GeV.mac $((G4SEED + 201)) && queued=$((queued + 1))
launch_had_gaps proton_supplement_run2 proton_supp_300GeV.mac $((G4SEED + 202)) && queued=$((queued + 1))
launch_had_gaps proton_supplement_run3 proton_supp_1000GeV.mac $((G4SEED + 203)) && queued=$((queued + 1))

echo "Launched up to $queued new shards (cap $MAX_NEW). Still running: $(running_jobs)"
