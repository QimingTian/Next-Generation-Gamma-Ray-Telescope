#!/usr/bin/env python3
"""Emit JSON status for cloud campaigns (run on MaxCloudON Linux VM)."""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(os.environ.get("GRT_ROOT", Path.home() / "GammaRayTelescope"))
BUILD_DATA = ROOT / "build" / "data"
DATA = ROOT / "data"
PHASE2_TOTAL = 216_000


def count_events_in_files(paths: list[Path]) -> int:
    """Count data rows in Geant4 wcsv ntuples (skip # comment header lines)."""
    total = 0
    for p in paths:
        try:
            with p.open() as fh:
                total += sum(
                    1 for line in fh if line.strip() and not line.lstrip().startswith("#")
                )
        except OSError:
            continue
    return total


def shard_progress(energy: int, shards_total: int = 128) -> dict:
    prefix = f"aeff_energy_{energy}GeV_shard"
    files = sorted(BUILD_DATA.glob(f"{prefix}*_nt_events.csv"))
    if not files:
        files = sorted(DATA.glob(f"{prefix}*_nt_events.csv"))
    events = min(1200, count_events_in_files(files))
    shards_started = len(files)
    return {
        "energy_GeV": energy,
        "shards_total": shards_total,
        "shards_started": shards_started,
        "shards_done": shards_started,
        "events_done": events,
        "events_target": 1200,
        "pct": round(100.0 * events / 1200, 1),
    }


def phase2_progress() -> dict:
    files = sorted(BUILD_DATA.glob("direction_scan_*_run0_nt_events.csv"))
    if not files:
        files = sorted(DATA.glob("direction_scan_*_run0_nt_events.csv"))
    events = count_events_in_files(files)
    jobs_total = 1728
    jobs_started = len(files)
    log_path = DATA / "maxcloud_phase2.log"
    done = False
    if log_path.is_file():
        try:
            done = "=== Done phase2" in log_path.read_text(errors="replace")
        except OSError:
            pass
    return {
        "jobs_total": jobs_total,
        "jobs_started": jobs_started,
        "events_done": events,
        "events_target": PHASE2_TOTAL,
        "pct": round(min(100.0, 100.0 * events / PHASE2_TOTAL), 2),
        "campaign_done": done,
    }


def active_campaign() -> str:
    if (DATA / "maxcloud_phase2.log").is_file():
        try:
            text = (DATA / "maxcloud_phase2.log").read_text(errors="replace")
            if "=== MaxCloud campaign phase2" in text and "=== Done phase2" not in text:
                return "phase2"
            if "=== Done phase2" in text:
                return "phase2_done"
        except OSError:
            pass
    if (DATA / "maxcloud_gaps.log").is_file():
        try:
            text = (DATA / "maxcloud_gaps.log").read_text(errors="replace")
            if "=== Done gaps" not in text and "=== MaxCloud campaign gaps" in text:
                return "gaps"
        except OSError:
            pass
    # Heuristic: running Main macro path beats stale log files
    try:
        out = subprocess.check_output(["pgrep", "-a", "Main"], text=True, stderr=subprocess.DEVNULL)
        if "gaps/" in out:
            if "1000GeV" in out or "E1000" in out:
                return "gaps_e1000"
            return "gaps"
        if "phase2/" in out:
            return "phase2"
    except (subprocess.CalledProcessError, OSError):
        pass
    if (DATA / "maxcloud_gaps_e1000.log").is_file():
        try:
            text = (DATA / "maxcloud_gaps_e1000.log").read_text(errors="replace")
            if "=== MaxCloud campaign gaps_e1000" in text and "=== Done gaps_e1000" not in text:
                return "gaps_e1000"
        except OSError:
            pass
    return "idle"


def log_tail(name: str) -> list[str]:
    path = DATA / f"maxcloud_{name}.log"
    if not path.is_file():
        return []
    try:
        return path.read_text(errors="replace").splitlines()[-8:]
    except OSError:
        return []


def main() -> None:
    main_procs = 0
    try:
        out = subprocess.check_output(["pgrep", "-c", "Main"], text=True).strip()
        main_procs = int(out) if out else 0
    except (subprocess.CalledProcessError, ValueError):
        main_procs = 0

    load = os.getloadavg() if hasattr(os, "getloadavg") else (0.0, 0.0, 0.0)
    ncpu = os.cpu_count() or 1
    campaign = active_campaign()

    gaps_log = DATA / "maxcloud_gaps.log"
    gaps_done = False
    segfaults = 0
    for log in (gaps_log, DATA / "maxcloud_phase2.log"):
        if log.is_file():
            try:
                text = log.read_text(errors="replace")
                segfaults += len(re.findall(r"Segmentation fault", text))
                if "=== Done gaps" in text:
                    gaps_done = True
            except OSError:
                pass

    gaps = {
        "energies": [shard_progress(e) for e in (500, 1000)],
        "total_events": sum(shard_progress(e)["events_done"] for e in (500, 1000)),
        "total_target": 2400,
        "campaign_done": gaps_done
        or all(shard_progress(e)["events_done"] >= 1200 for e in (500, 1000)),
    }
    gaps["total_pct"] = round(min(100.0, 100.0 * gaps["total_events"] / 2400), 1)

    p2 = phase2_progress()

    if campaign.startswith("phase2"):
        total_pct = p2["pct"]
        total_events = p2["events_done"]
        total_target = p2["events_target"]
        campaign_done = p2["campaign_done"]
        tail = log_tail("phase2")
    elif campaign == "gaps_e1000":
        e1000 = next(x for x in gaps["energies"] if x["energy_GeV"] == 1000)
        total_pct = e1000["pct"]
        total_events = e1000["events_done"]
        total_target = e1000["events_target"]
        campaign_done = e1000["events_done"] >= 1200
        tail = log_tail("gaps_e1000")
    else:
        total_pct = gaps["total_pct"]
        total_events = gaps["total_events"]
        total_target = gaps["total_target"]
        campaign_done = gaps["campaign_done"]
        tail = log_tail("gaps")

    status = {
        "ok": True,
        "time_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
        "hostname": os.uname().nodename,
        "ncpu": ncpu,
        "load_avg": [round(x, 2) for x in load],
        "main_procs": main_procs,
        "cpu_util_pct": round(min(100.0, 100.0 * load[0] / ncpu), 1),
        "campaign": campaign,
        "campaign_done": campaign_done,
        "segfaults": segfaults,
        "gaps": gaps,
        "phase2": p2,
        "total_events": total_events,
        "total_target": total_target,
        "total_pct": total_pct,
        "log_tail": tail,
    }
    sys.stdout.write(json.dumps(status) + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
