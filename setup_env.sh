#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/deps/geant4-install/bin/geant4.sh"
export Geant4_DIR="${SCRIPT_DIR}/deps/geant4-install/lib/cmake/Geant4"
export CMAKE_PREFIX_PATH="/opt/homebrew/opt/expat:/opt/homebrew/opt/qt@5:${CMAKE_PREFIX_PATH:-}"
export G4SEED="${G4SEED:-123456789}"
