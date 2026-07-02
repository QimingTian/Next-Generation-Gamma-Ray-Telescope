#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TECTONIC_DIR="${SCRIPT_DIR}/deps/tectonic"
TECTONIC_BIN="${TECTONIC_DIR}/bin/tectonic"
TECTONIC_VERSION="0.16.9"
ARCH="$(uname -m)"

case "${ARCH}" in
  arm64|aarch64)
    TECTONIC_URL="https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.16.9/tectonic-${TECTONIC_VERSION}-aarch64-apple-darwin.tar.gz"
    ;;
  x86_64)
    TECTONIC_URL="https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.16.9/tectonic-${TECTONIC_VERSION}-x86_64-apple-darwin.tar.gz"
    ;;
  *)
    echo "Unsupported architecture: ${ARCH}" >&2
    echo "Install MacTeX or BasicTeX manually, then run paper/compile_pdflatex.sh" >&2
    exit 1
    ;;
esac

if [[ -x "${TECTONIC_BIN}" ]]; then
  echo "Tectonic already installed: $("${TECTONIC_BIN}" --version)"
  exit 0
fi

mkdir -p "${TECTONIC_DIR}/bin"
TMP="${TECTONIC_DIR}/tectonic.tar.gz"

echo "Downloading Tectonic ${TECTONIC_VERSION} for ${ARCH}..."
curl -fsSL -o "${TMP}" "${TECTONIC_URL}"
tar -xzf "${TMP}" -C "${TECTONIC_DIR}"
mv "${TECTONIC_DIR}/tectonic" "${TECTONIC_BIN}"
chmod +x "${TECTONIC_BIN}"
rm -f "${TMP}"

echo "Installed: $("${TECTONIC_BIN}" --version)"
echo "Compile the paper with: ${SCRIPT_DIR}/paper/compile.sh"
