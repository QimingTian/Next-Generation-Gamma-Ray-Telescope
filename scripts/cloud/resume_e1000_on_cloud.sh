#!/usr/bin/env bash
# Stop Phase 2 and start E1000 gaps补跑 on MaxCloudON (via Parallels Windows SSH).
set -euo pipefail
CLOUD_IP="${GRT_CLOUD_IP:-10.101.58.1}"
PRL_VM="${GRT_PARALLELS_VM:-Windows 11}"
SSH_KEY="${GRT_SSH_KEY_WIN:-C:\\Windows\\Temp\\maxcloud_ssh\\id_ed25519}"

REMOTE='pkill -9 Main 2>/dev/null || true
pkill -f run_maxcloud_pipeline 2>/dev/null || true
pkill -f "run_maxcloud_campaign.sh phase2" 2>/dev/null || true
sleep 2
cd /root/GammaRayTelescope
git pull -q
chmod +x scripts/cloud/*.sh scripts/cloud/*.py
nohup bash scripts/cloud/run_maxcloud_campaign.sh gaps_e1000 \
  > data/maxcloud_gaps_e1000_nohup.log 2>&1 &
sleep 8
echo "Main_procs=$(pgrep -c Main || echo 0)"
head -3 data/maxcloud_gaps_e1000.log 2>/dev/null || true'

echo "=== Stop Phase 2, launch gaps_e1000 on $CLOUD_IP ==="
prlctl exec "$PRL_VM" cmd /c "ssh -o BatchMode=yes -i $SSH_KEY root@$CLOUD_IP \"$REMOTE\""
