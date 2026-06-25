#!/usr/bin/env python3
"""Localhost dashboard for MaxCloudON Geant4 campaign progress."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

PORT = int(os.environ.get("GRT_DASHBOARD_PORT", "8765"))
CLOUD_IP = os.environ.get("GRT_CLOUD_IP", "10.101.58.1")
CLOUD_USER = os.environ.get("GRT_CLOUD_USER", "root")
PRL_VM = os.environ.get("GRT_PARALLELS_VM", "Windows 11")
SSH_KEY = os.environ.get(
    "GRT_SSH_KEY_WIN",
    r"C:\Windows\Temp\maxcloud_ssh\id_ed25519",
)
ROOT = Path(__file__).resolve().parents[2]

HTML = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>GammaRayTelescope — Cloud Progress</title>
  <style>
    :root {
      --bg: #0f1419;
      --card: #1a2332;
      --text: #e6edf3;
      --muted: #8b949e;
      --accent: #58a6ff;
      --ok: #3fb950;
      --warn: #d29922;
      --err: #f85149;
    }
    * { box-sizing: border-box; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--text);
      margin: 0;
      padding: 24px;
      line-height: 1.5;
    }
    h1 { font-size: 1.35rem; margin: 0 0 4px; }
    .sub { color: var(--muted); font-size: 0.9rem; margin-bottom: 20px; }
    .grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 16px;
      margin-bottom: 20px;
    }
    .card {
      background: var(--card);
      border-radius: 10px;
      padding: 16px 18px;
      border: 1px solid #30363d;
    }
    .card h2 { font-size: 0.85rem; color: var(--muted); margin: 0 0 10px; text-transform: uppercase; letter-spacing: .04em; }
    .big { font-size: 2rem; font-weight: 600; }
    .bar-wrap {
      background: #21262d;
      border-radius: 6px;
      height: 12px;
      overflow: hidden;
      margin: 8px 0;
    }
    .bar { height: 100%; background: linear-gradient(90deg, #238636, var(--ok)); transition: width .4s; }
    .row { display: flex; justify-content: space-between; font-size: 0.9rem; }
    .pill {
      display: inline-block;
      padding: 2px 10px;
      border-radius: 999px;
      font-size: 0.8rem;
      font-weight: 600;
    }
    .pill.ok { background: #23863633; color: var(--ok); }
    .pill.err { background: #f8514933; color: var(--err); }
    .pill.wait { background: #d2992233; color: var(--warn); }
    pre {
      background: #0d1117;
      border-radius: 8px;
      padding: 12px;
      font-size: 0.75rem;
      overflow-x: auto;
      color: var(--muted);
      max-height: 200px;
    }
    #err { color: var(--err); margin-top: 12px; display: none; }
  </style>
</head>
<body>
  <h1>☁️ MaxCloudON Campaign Monitor</h1>
  <p class="sub">Auto-refresh every 5s · <span id="updated">—</span></p>
  <p><span id="conn" class="pill wait">connecting…</span> <span id="campaign-label" class="pill wait">—</span></p>

  <div class="grid">
    <div class="card">
      <h2 id="overall-title">Active campaign</h2>
      <div class="big" id="total-pct">—</div>
      <div class="bar-wrap"><div class="bar" id="total-bar" style="width:0%"></div></div>
      <div class="row"><span id="total-events">—</span><span id="campaign-state">—</span></div>
    </div>
    <div class="card">
      <h2>CPU</h2>
      <div class="big" id="main-procs">—</div>
      <div class="row"><span id="load">load —</span><span id="ncpu">— cores</span></div>
      <div class="bar-wrap"><div class="bar" id="cpu-bar" style="width:0%; background: linear-gradient(90deg,#1f6feb,var(--accent))"></div></div>
    </div>
    <div class="card" id="card-gaps">
      <h2>Gaps · E500</h2>
      <div class="big" id="e500-pct">—</div>
      <div class="bar-wrap"><div class="bar" id="e500-bar" style="width:0%"></div></div>
      <div class="row"><span id="e500-ev">—</span><span id="e500-sh">—</span></div>
    </div>
    <div class="card" id="card-gaps2">
      <h2>Gaps · E1000</h2>
      <div class="big" id="e1000-pct">—</div>
      <div class="bar-wrap"><div class="bar" id="e1000-bar" style="width:0%"></div></div>
      <div class="row"><span id="e1000-ev">—</span><span id="e1000-sh">—</span></div>
    </div>
    <div class="card" id="card-phase2">
      <h2>Phase 2 · direction scan</h2>
      <div class="big" id="p2-pct">—</div>
      <div class="bar-wrap"><div class="bar" id="p2-bar" style="width:0%; background: linear-gradient(90deg,#8957e5,#a371f7)"></div></div>
      <div class="row"><span id="p2-ev">—</span><span id="p2-jobs">—</span></div>
    </div>
  </div>

  <div class="card">
    <h2>Log tail</h2>
    <pre id="log">—</pre>
  </div>
  <div id="err"></div>

<script>
async function refresh() {
  try {
    const r = await fetch('/api/status');
    const d = await r.json();
    const errEl = document.getElementById('err');
    errEl.style.display = 'none';

    if (!d.ok) {
      document.getElementById('conn').className = 'pill err';
      document.getElementById('conn').textContent = d.error || 'offline';
      return;
    }

    document.getElementById('conn').className = 'pill ok';
    document.getElementById('conn').textContent = d.cloud.hostname + ' · connected';
    document.getElementById('updated').textContent = d.cloud.time_utc || d.fetched_local;

    const c = d.cloud;
    const camp = c.campaign || 'gaps';
    const label = document.getElementById('campaign-label');
    label.textContent = camp;
    label.className = 'pill ' + (camp.startsWith('phase2') ? 'ok' : (camp === 'gaps' ? 'wait' : 'ok'));

    document.getElementById('overall-title').textContent =
      camp.startsWith('phase2') ? 'Phase 2 overall' : 'Gaps overall';

    document.getElementById('total-pct').textContent = c.total_pct + '%';
    document.getElementById('total-bar').style.width = Math.min(100, c.total_pct) + '%';
    document.getElementById('total-events').textContent =
      c.total_events.toLocaleString() + ' / ' + c.total_target.toLocaleString() + ' events';
    document.getElementById('campaign-state').textContent =
      c.campaign_done ? '✅ done' : (c.main_procs > 0 ? '▶ running' : '⏸ idle');

    document.getElementById('main-procs').textContent = c.main_procs + ' Main';
    document.getElementById('load').textContent =
      'load ' + (c.load_avg || []).join(' / ');
    document.getElementById('ncpu').textContent = c.ncpu + ' cores · ' + c.cpu_util_pct + '%';
    document.getElementById('cpu-bar').style.width = Math.min(100, c.cpu_util_pct) + '%';

    const gaps = c.gaps || {};
    for (const [e, ids] of [[500, ['e500']], [1000, ['e1000']]]) {
      const row = (gaps.energies || []).find(x => x.energy_GeV === e) || {};
      document.getElementById(ids[0] + '-pct').textContent = (row.pct || 0) + '%';
      document.getElementById(ids[0] + '-bar').style.width = (row.pct || 0) + '%';
      document.getElementById(ids[0] + '-ev').textContent =
        (row.events_done || 0) + ' / ' + (row.events_target || 1200) + ' events';
      document.getElementById(ids[0] + '-sh').textContent =
        (row.shards_done || 0) + ' / ' + (row.shards_total || 128) + ' shards';
    }

    const p2 = c.phase2 || {};
    document.getElementById('p2-pct').textContent = (p2.pct || 0) + '%';
    document.getElementById('p2-bar').style.width = Math.min(100, p2.pct || 0) + '%';
    document.getElementById('p2-ev').textContent =
      (p2.events_done || 0).toLocaleString() + ' / ' + (p2.events_target || 216000).toLocaleString();
    document.getElementById('p2-jobs').textContent =
      (p2.jobs_started || 0) + ' / ' + (p2.jobs_total || 1728) + ' jobs';

    document.getElementById('card-phase2').style.opacity = camp.startsWith('phase2') ? '1' : '0.55';

    document.getElementById('log').textContent =
      (c.log_tail && c.log_tail.length) ? c.log_tail.join('\\n') : '(no log yet)';
  } catch (e) {
    const errEl = document.getElementById('err');
    errEl.style.display = 'block';
    errEl.textContent = 'Fetch error: ' + e;
    document.getElementById('conn').className = 'pill err';
    document.getElementById('conn').textContent = 'error';
  }
}
refresh();
setInterval(refresh, 5000);
</script>
</body>
</html>
"""


def fetch_cloud_status() -> dict:
    remote = (
        "python3 /root/GammaRayTelescope/scripts/cloud/cloud_status.py "
        "2>/dev/null || python3 ~/GammaRayTelescope/scripts/cloud/cloud_status.py"
    )
    remote_escaped = remote.replace('"', '\\"')
    win_cmd = (
        f'ssh -o ConnectTimeout=20 -o BatchMode=yes -i {SSH_KEY} '
        f'{CLOUD_USER}@{CLOUD_IP} "{remote_escaped}"'
    )
    try:
        proc = subprocess.run(
            ["prlctl", "exec", PRL_VM, "cmd", "/c", win_cmd],
            capture_output=True,
            text=True,
            timeout=90,
        )
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "SSH timeout (VPN down or VM overloaded?)"}
    except FileNotFoundError:
        return {"ok": False, "error": "prlctl not found — need Parallels Windows VM"}

    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()[:400]
        return {"ok": False, "error": err or f"ssh exit {proc.returncode}"}

    text = proc.stdout.strip()
    if not text:
        return {"ok": False, "error": "empty response from cloud"}

    try:
        return {"ok": True, "cloud": json.loads(text.splitlines()[-1])}
    except json.JSONDecodeError as exc:
        return {"ok": False, "error": f"bad json: {exc}", "raw": text[:500]}


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt: str, *args) -> None:
        pass

    def do_GET(self) -> None:
        path = urllib.parse.urlparse(self.path).path
        if path == "/api/status":
            from datetime import datetime

            payload = fetch_cloud_status()
            if payload.get("ok"):
                payload["fetched_local"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S local")
            body = json.dumps(payload).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)
            return

        if path in ("/", "/index.html"):
            body = HTML.encode()
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(body)
            return

        self.send_error(404)


def main() -> None:
    print(f"Dashboard: http://127.0.0.1:{PORT}/")
    print(f"Cloud: {CLOUD_USER}@{CLOUD_IP} via Parallels '{PRL_VM}'")
    print("Ctrl+C to stop")
    server = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
        sys.exit(0)


if __name__ == "__main__":
    main()
