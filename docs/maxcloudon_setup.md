# MaxCloudON setup (GammaRayTelescope)

## Network note

VPN runs in your **Windows VM** (Sync Manager Connected). Your **Mac is not on the VPN** unless you add a public IP or route through Windows.

**Workflow:** SSH from **Windows VM** → MaxCloudON Linux VM.

## Dashboard checklist

1. [acc.maxcloudon.com](https://acc.maxcloudon.com/login) → **Virtual Machines**
2. Create VM if needed: **Linux Ubuntu 22.04**, 128 vCPU / 64 GB plan
3. **Power ON** (no grey bar on top)
4. Note **internal IP** and **username/password** on dashboard
5. **Profile → SSH Public Key** → paste key from Windows VM:
   ```powershell
   # In Windows VM PowerShell (if no key yet):
   ssh-keygen -t ed25519 -N '""' -f $env:USERPROFILE\.ssh\id_ed25519
   type $env:USERPROFILE\.ssh\id_ed25519.pub
   ```
6. **Reboot VM** after adding SSH key (MaxCloudON FAQ)

## Connect from Windows VM

```powershell
ssh root@<INTERNAL_IP>
# or ubuntu@<INTERNAL_IP> — use dashboard credentials
```

## Upload project (from Mac → cloud via Windows)

Option A — **git** (if repo is on GitHub):
```bash
git clone <your-repo-url> ~/GammaRayTelescope
```

Option B — **rsync from Mac through Windows** (two hops) or zip upload.

Option C — **rsync deps + code from Mac** (recommended, includes Geant4 install tree):

On **Mac** (after VPN path works from Windows, run rsync *on Windows* with WinSCP, or from Mac if public IP):

```bash
# If you get public IP $IP:
rsync -avz --progress \
  --exclude '.git' --exclude 'build' --exclude 'build-hadronic' \
  ~/Documents/GammaRayTelescope/ root@$IP:~/GammaRayTelescope/
```

Geant4 `deps/` is large (~several GB) but saves 1–2 h compile on cloud.

## Bootstrap & run

```bash
cd ~/GammaRayTelescope
bash scripts/cloud/maxcloud_bootstrap.sh
bash scripts/cloud/run_maxcloud_campaign.sh gaps   # E500+E1000
bash scripts/cloud/run_maxcloud_campaign.sh phase2 # 216k direction scan
```

## Mac SSH public key (optional, for public IP later)

```
ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIJSAOYGVx8MrR30lh6CD5j6vl6vCkVzZiNUcZyMQT7lc qtian.28@pomfret.org
```
