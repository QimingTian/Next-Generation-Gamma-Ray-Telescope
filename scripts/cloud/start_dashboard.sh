#!/usr/bin/env bash
# Local progress dashboard → http://127.0.0.1:8765/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
exec python3 scripts/cloud/progress_server.py
