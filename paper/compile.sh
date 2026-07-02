#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TECTONIC="${REPO_ROOT}/deps/tectonic/bin/tectonic"

if [[ ! -x "${TECTONIC}" ]]; then
  echo "Tectonic not found at ${TECTONIC}" >&2
  echo "Run: ${REPO_ROOT}/setup_tex.sh" >&2
  exit 1
fi

cd "${SCRIPT_DIR}"
export PATH="${REPO_ROOT}/deps/tectonic/bin:${PATH}"

echo "Compiling Paper.tex with Tectonic..."
"${TECTONIC}" -X compile --keep-logs --keep-intermediates Paper.tex

if [[ -f Paper.pdf ]]; then
  echo "Success: ${SCRIPT_DIR}/Paper.pdf"
else
  echo "Compilation finished but Paper.pdf not found." >&2
  exit 1
fi
