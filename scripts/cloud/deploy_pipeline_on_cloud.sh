#!/usr/bin/env bash
# Watch gaps completion on cloud, then auto-launch Phase 2 (resume incomplete shards).
set -euo pipefail
CLOUD_IP="${GRT_CLOUD_IP:-10.101.58.1}"
PRL_VM="${GRT_PARALLELS_VM:-Windows 11}"
SSH_KEY="${GRT_SSH_KEY_WIN:-C:\\Windows\\Temp\\maxcloud_ssh\\id_ed25519}"

REMOTE='cd /root/GammaRayTelescope
git fetch -q origin
git reset -q --hard origin/main
chmod +x scripts/cloud/*.sh scripts/cloud/*.py
pkill -f run_maxcloud_pipeline 2>/dev/null || true
nohup bash scripts/cloud/run_maxcloud_pipeline.sh > data/maxcloud_pipeline_nohup.log 2>&1 </dev/null &
sleep 5
tail -3 data/maxcloud_pipeline.log 2>/dev/null || true'

echo "=== Deploy pipeline watcher (gaps -> Phase 2) on $CLOUD_IP ==="
prlctl exec "$PRL_VM" cmd /c "ssh -o BatchMode=yes -i $SSH_KEY root@$CLOUD_IP bash -lc \"$REMOTE\""
