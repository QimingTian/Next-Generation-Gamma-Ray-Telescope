#!/usr/bin/env bash
# Run one Geant4 shard (serial process). Args: TAG MACRO FILTER [SEED]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${BUILD_DIR:-$ROOT/build}"
TAG="$1"
MACRO="$2"
FILTER="${3:-on}"
SEED="${4:-${G4SEED:-123456789}}"

source "$ROOT/setup_env.sh"
export G4MT=0
export G4WRITE_PHOTONS="${G4WRITE_PHOTONS:-0}"
export G4OUTPUT_TAG="$TAG"
export G4MACRO="$MACRO"
export G4FILTER="$FILTER"
export G4SEED="$SEED"

mkdir -p "$BUILD/macros/parallel" "$BUILD/data" "$ROOT/data"

if [ -f "$ROOT/macros/parallel/$MACRO" ]; then
  cp -f "$ROOT/macros/parallel/$MACRO" "$BUILD/macros/parallel/"
  MACRO_PATH="macros/parallel/$MACRO"
elif [ -f "$ROOT/macros/gaps/$MACRO" ]; then
  mkdir -p "$BUILD/macros/gaps"
  cp -f "$ROOT/macros/gaps/$MACRO" "$BUILD/macros/gaps/"
  MACRO_PATH="macros/gaps/$MACRO"
elif [[ "$MACRO" == gaps/* ]] && [ -f "$ROOT/macros/$MACRO" ]; then
  mkdir -p "$BUILD/macros/gaps"
  cp -f "$ROOT/macros/$MACRO" "$BUILD/macros/gaps/${MACRO#gaps/}"
  MACRO_PATH="macros/gaps/${MACRO#gaps/}"
elif [ -f "$ROOT/macros/$MACRO" ]; then
  cp -f "$ROOT/macros/$MACRO" "$BUILD/macros/"
  MACRO_PATH="macros/$MACRO"
else
  echo "Macro not found: $MACRO" >&2
  exit 1
fi

(
  cd "$BUILD"
  ./Main "$MACRO_PATH"
)
cp -f "$BUILD/data/${TAG}_"* "$ROOT/data/" 2>/dev/null || true
