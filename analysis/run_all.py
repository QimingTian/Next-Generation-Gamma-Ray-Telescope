#!/usr/bin/env python3
"""Run all analysis scripts after a Geant4 campaign."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
if not any(DATA.glob("*_nt_events.csv")):
    DATA = ROOT / "build" / "data"

scripts = [
    [sys.executable, "analysis/Filter-Design.py"],
    [sys.executable, "analysis/density_maps.py", "--tag", "peak_factor"],
    [sys.executable, "analysis/energy_calibration.py", "--data-dir", str(DATA)],
    [sys.executable, "analysis/energy_resolution.py", "--tag", "energy_resolution_215_on"],
    [sys.executable, "analysis/energy_resolution_filter_compare.py", "--data-dir", str(DATA)],
    [sys.executable, "analysis/shower_length.py", "--data-dir", str(DATA)],
    [sys.executable, "analysis/peak_factor_stats.py", "--tag", "peak_factor"],
    [sys.executable, "analysis/filter_folding.py", "--tag", "peak_factor"],
    [sys.executable, "analysis/saturation.py"],
    [sys.executable, "analysis/angular_resolution.py", "--tag", "direction_scan", "--train-route-b"],
    [sys.executable, "analysis/effective_area.py", "--tag", "plane_source", "--data-dir", str(DATA)],
    [sys.executable, "analysis/background_rejection.py", "--data-dir", str(DATA)],
]

if __name__ == "__main__":
    for cmd in scripts:
        print(">>", " ".join(cmd))
        subprocess.check_call(cmd, cwd=ROOT)
