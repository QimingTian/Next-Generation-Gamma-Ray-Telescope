#!/usr/bin/env bash
# Pull Phase 2 direction_scan CSVs from MaxCloudON via Parallels Windows SSH.
# Usage: fetch_phase2_data.sh [--energy 30|100|300]   (default: all energies on disk)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data"
CLOUD_IP="${GRT_CLOUD_IP:-10.101.58.1}"
PRL_VM="${GRT_PARALLELS_VM:-Windows 11}"
SSH_KEY="${GRT_SSH_KEY_WIN:-C:\\Windows\\Temp\\maxcloud_ssh\\id_ed25519}"
REMOTE_TAR="/root/GammaRayTelescope/data/phase2_events.tar.gz"
WIN_DEST='C:\Mac\Home\Documents\GammaRayTelescope\data\phase2_events.tar.gz'

ENERGY_FILTER=""
if [ "${1:-}" = "--energy" ] && [ -n "${2:-}" ]; then
  ENERGY_FILTER="$2"
fi

if [ -n "$ENERGY_FILTER" ]; then
  PATTERN="direction_scan_E${ENERGY_FILTER}_*_nt_*.csv"
  echo "=== Pack Phase 2 E=${ENERGY_FILTER} GeV on cloud ==="
else
  PATTERN="direction_scan_*_nt_*.csv"
  echo "=== Pack all Phase 2 direction_scan CSVs on cloud ==="
fi

PACK_CMD="cd /root/GammaRayTelescope/build/data && shopt -s nullglob
files=($PATTERN)
if ((\${#files[@]}==0)); then echo 'No files match $PATTERN' >&2; exit 1; fi
tar czf $REMOTE_TAR \${files[@]}
ls -lh $REMOTE_TAR
echo file_count=\${#files[@]}"

prlctl exec "$PRL_VM" cmd /c "ssh -o BatchMode=yes -i $SSH_KEY root@$CLOUD_IP bash -lc \"$PACK_CMD\""

echo "=== SCP to Mac (via Windows share) ==="
prlctl exec "$PRL_VM" cmd /c "scp -o BatchMode=yes -i $SSH_KEY root@$CLOUD_IP:$REMOTE_TAR $WIN_DEST"

echo "=== Extract into data/ ==="
mkdir -p "$DATA"
tar xzf "$DATA/phase2_events.tar.gz" -C "$DATA"
nev="$(ls "$DATA"/direction_scan_*_nt_events.csv 2>/dev/null | wc -l | tr -d ' ')"
nsm="$(ls "$DATA"/direction_scan_*_nt_sipm_summary.csv 2>/dev/null | wc -l | tr -d ' ')"
echo "Extracted $nev event CSVs, $nsm sipm_summary CSVs -> $DATA"
