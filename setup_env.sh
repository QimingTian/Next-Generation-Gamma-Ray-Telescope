#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/deps/geant4-install/bin/geant4.sh"
export Geant4_DIR="${SCRIPT_DIR}/deps/geant4-install/lib/cmake/Geant4"
export CMAKE_PREFIX_PATH="/opt/homebrew/opt/expat:/opt/homebrew/opt/qt@5:${CMAKE_PREFIX_PATH:-}"
export G4SEED="${G4SEED:-123456789}"

# Geant4 honours G4FORCENUMBEROFTHREADS over SetNumberOfThreads(); use all logical cores.
unset G4FORCE_RUN_MANAGER_TYPE
unset G4FORCENUMBEROFTHREADS
_NCPU="$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 8)"
export G4FORCENUMBEROFTHREADS="${G4NUM_THREADS:-$_NCPU}"
export G4NUM_THREADS="${G4NUM_THREADS:-$G4FORCENUMBEROFTHREADS}"
export G4ACD="${G4ACD:-on}"
