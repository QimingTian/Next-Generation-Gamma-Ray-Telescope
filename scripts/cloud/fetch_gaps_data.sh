#!/usr/bin/env bash
# Pull gaps A_eff CSV shards from MaxCloudON via Parallels Windows SSH.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DATA="$ROOT/data"
CLOUD_IP="${GRT_CLOUD_IP:-10.101.58.1}"
PRL_VM="${GRT_PARALLELS_VM:-Windows 11}"
SSH_KEY="${GRT_SSH_KEY_WIN:-C:\\Windows\\Temp\\maxcloud_ssh\\id_ed25519}"
REMOTE_TAR="/root/GammaRayTelescope/data/gaps_events.tar.gz"
WIN_DEST='C:\Mac\Home\Documents\GammaRayTelescope\data\gaps_events.tar.gz'

echo "=== Pack gaps events on cloud ==="
prlctl exec "$PRL_VM" cmd /c "ssh -o BatchMode=yes -i $SSH_KEY root@$CLOUD_IP \"cd /root/GammaRayTelescope/build/data && tar czf $REMOTE_TAR aeff_energy_500GeV_shard*_nt_events.csv aeff_energy_1000GeV_shard*_nt_events.csv && ls -lh $REMOTE_TAR\""

echo "=== SCP to Mac (via Windows share) ==="
prlctl exec "$PRL_VM" cmd /c "scp -o BatchMode=yes -i $SSH_KEY root@$CLOUD_IP:$REMOTE_TAR $WIN_DEST"

echo "=== Extract into data/ ==="
mkdir -p "$DATA"
tar xzf "$DATA/gaps_events.tar.gz" -C "$DATA"
n="$(ls "$DATA"/aeff_energy_*GeV_shard*_nt_events.csv 2>/dev/null | wc -l | tr -d ' ')"
echo "Extracted $n event CSV files to $DATA"
