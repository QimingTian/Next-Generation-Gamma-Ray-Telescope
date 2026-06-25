# Run inside Parallels Windows VM (VPN Connected). Uses Mac shared folder Z:\ or C:\Mac\Home.
param(
  [string]$CloudIP = "10.101.58.1",
  [string]$CloudUser = "root",
  [string]$KeyPath = "$env:TEMP\maxcloud_ssh\id_ed25519",
  [string]$MacRepo = "C:\Mac\Home\Documents\GammaRayTelescope"
)

$ErrorActionPreference = "Stop"
$ssh = "ssh -o StrictHostKeyChecking=accept-new -o BatchMode=yes -i `"$KeyPath`""
$scp = "scp -o StrictHostKeyChecking=accept-new -o BatchMode=yes -i `"$KeyPath`""

function Invoke-Remote([string]$Cmd) {
  Invoke-Expression "$ssh ${CloudUser}@${CloudIP} `"$Cmd`""
}

Write-Host "=== Ping cloud ==="
ping -n 2 $CloudIP
if ($LASTEXITCODE -ne 0) { throw "Cannot reach $CloudIP — reconnect MaxCloudON VPN (Sync Manager -> Connected)." }

Write-Host "=== Ensure repo on cloud ==="
Invoke-Remote "test -f ~/GammaRayTelescope/CMakeLists.txt || git clone --depth 1 https://github.com/QimingTian/Next-Generation-Gamma-Ray-Telescope.git ~/GammaRayTelescope"

Write-Host "=== Sync cloud scripts from Mac (small) ==="
Invoke-Remote "mkdir -p ~/GammaRayTelescope/scripts/cloud ~/GammaRayTelescope/data"
Invoke-Expression "$scp -r `"$MacRepo\scripts\cloud\*`" ${CloudUser}@${CloudIP}:~/GammaRayTelescope/scripts/cloud/"

Write-Host "=== Geant4 build (background if not installed) ==="
Invoke-Remote "bash ~/GammaRayTelescope/scripts/cloud/build_geant4.sh"
if ($LASTEXITCODE -ne 0) {
  Write-Host "Geant4 build running or failed — check: ssh ... 'tail -20 ~/GammaRayTelescope/data/geant4_build.log'"
}

Write-Host "=== Bootstrap project ==="
Invoke-Remote "cd ~/GammaRayTelescope && bash scripts/cloud/maxcloud_bootstrap.sh"

Write-Host "=== Start gaps campaign (E500+E1000) ==="
Invoke-Remote "cd ~/GammaRayTelescope && nohup bash scripts/cloud/run_maxcloud_campaign.sh gaps > data/maxcloud_gaps_nohup.log 2>&1 &"

Write-Host "=== Deploy started. Monitor: ==="
Write-Host "  ssh ${CloudUser}@${CloudIP} 'tail -f ~/GammaRayTelescope/data/maxcloud_gaps.log'"
