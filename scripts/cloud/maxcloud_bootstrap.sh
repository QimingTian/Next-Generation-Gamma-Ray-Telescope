#!/usr/bin/env bash
# First-time setup on MaxCloudON Linux VM (Ubuntu 22.04+).
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -qq
sudo apt-get install -y -qq \
  build-essential cmake git rsync python3 python3-pip python3-venv \
  libexpat1-dev qtbase5-dev libxmu-dev libxi-dev libglu1-mesa-dev \
  libxerces-c-dev

REPO="${1:-$HOME/GammaRayTelescope}"
if [ ! -d "$REPO/.git" ] && [ ! -f "$REPO/CMakeLists.txt" ]; then
  echo "Clone or rsync the repo to $REPO first, then re-run." >&2
  exit 1
fi

cd "$REPO"
if [ ! -f deps/geant4-install/bin/geant4.sh ]; then
  echo "Geant4 not found — building (see scripts/cloud/build_geant4.sh)..."
  bash "$REPO/scripts/cloud/build_geant4.sh"
fi

# Mac-only paths in setup_env.sh; harmless on Linux if dirs missing.
if [ "$(uname -s)" = "Linux" ]; then
  export CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-}"
fi

python3 -m venv .venv
source .venv/bin/activate
pip install -q -r analysis/requirements.txt 2>/dev/null || pip install -q numpy pandas matplotlib scikit-learn

mkdir -p build build-hadronic
source setup_env.sh 2>/dev/null || true
export G4MT=0 G4WRITE_PHOTONS=0 G4ACD=on

if [ -f deps/geant4-install/bin/geant4.sh ]; then
  (cd build && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=OFF && cmake --build . -j"$(nproc)")
  (cd build-hadronic && cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=ON && cmake --build . -j"$(nproc)" MainHad)
  echo "Build OK: build/Main, build-hadronic/MainHad"
else
  echo "Skipped Geant4 build — upload deps first."
fi

echo "Bootstrap done. Next: bash scripts/cloud/run_maxcloud_phase2.sh"
