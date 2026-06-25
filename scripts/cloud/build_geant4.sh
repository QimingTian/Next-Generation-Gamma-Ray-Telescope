#!/usr/bin/env bash
# Build Geant4 11.3.2 on MaxCloudON Linux VM (headless, no Qt).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DEPS="$ROOT/deps"
VER=11.3.2
TARBALL="$DEPS/geant4-v${VER}.tar.gz"
SRC="$DEPS/geant4-${VER}"
BUILD="$DEPS/geant4-build"
INSTALL="$DEPS/geant4-install"
LOG="$ROOT/data/geant4_build.log"

mkdir -p "$ROOT/data" "$DEPS" "$BUILD" "$INSTALL"
exec > >(tee -a "$LOG") 2>&1

echo "=== Geant4 $VER build $(date) ncpu=$(nproc) ==="

if [ -f "$INSTALL/bin/geant4.sh" ]; then
  echo "Already installed at $INSTALL"
  exit 0
fi

if [ ! -d "$SRC" ]; then
  if [ ! -f "$TARBALL" ]; then
    wget -q "https://github.com/Geant4/geant4/archive/refs/tags/v${VER}.tar.gz" -O "$TARBALL"
  fi
  tar xzf "$TARBALL" -C "$DEPS"
fi

cmake -S "$SRC" -B "$BUILD" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL" \
  -DGEANT4_USE_QT=OFF \
  -DGEANT4_USE_OPENGL_X11=OFF \
  -DGEANT4_INSTALL_DATA=ON \
  -DGEANT4_BUILD_MULTITHREADED=ON \
  -DGEANT4_USE_SYSTEM_EXPAT=ON

cmake --build "$BUILD" -j"$(nproc)"
cmake --install "$BUILD"

echo "GEANT4_BUILD_DONE $(date)" > "$ROOT/data/geant4_build.done"
echo "=== Done ==="
