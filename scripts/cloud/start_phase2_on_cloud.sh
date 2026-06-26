#!/usr/bin/env bash
# Start Phase 2 on MaxCloudON (skips completed direction_scan shards).
set -euo pipefail
CLOUD_IP="${GRT_CLOUD_IP:-10.101.58.1}"
PRL_VM="${GRT_PARALLELS_VM:-Windows 11}"
SSH_KEY="${GRT_SSH_KEY_WIN:-C:\\Windows\\Temp\\maxcloud_ssh\\id_ed25519}"

REMOTE='cd /root/GammaRayTelescope
git fetch -q origin
git reset -q --hard origin/main
chmod +x scripts/cloud/*.sh scripts/cloud/*.py
pkill -9 Main 2>/dev/null || true
pkill -f run_maxcloud_campaign 2>/dev/null || true
sleep 2
nohup bash scripts/cloud/run_maxcloud_campaign.sh phase2 \
  > data/maxcloud_phase2_nohup.log 2>&1 </dev/null &
sleep 10
echo Main_procs=$(pgrep -c Main || echo 0)
head -3 data/maxcloud_phase2.log 2>/dev/null || true'

echo "=== Launch Phase 2 on $CLOUD_IP ==="
prlctl exec "$PRL_VM" cmd /c "ssh -o BatchMode=yes -i $SSH_KEY root@$CLOUD_IP bash -lc \"$REMOTE\""
