"""Load and merge Geant4 run CSVs with unique event IDs."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import parse_geant4_csv


def truth_direction_direction_scan(run_index: int, phi_deg: float = 0.0) -> np.ndarray:
    """Truth momentum for direction_scan.mac layout (GPS rot1 about Y)."""
    energies = [30, 100, 300]
    thetas = [0, 5, 10, 20, 30, 45]
    configs = [(e, t) for e in energies for t in thetas]
    if run_index >= len(configs):
        return np.array([0.0, 0.0, -1.0])
    _, theta = configs[run_index]
    if theta == 0:
        return np.array([0.0, 0.0, -1.0])
    th = np.radians(theta)
    ph = np.radians(phi_deg)
    # Rotation Ry(-theta) then Rz(phi) on (0,0,-1).
    x = -np.sin(th) * np.cos(ph)
    y = -np.sin(th) * np.sin(ph)
    z = -np.cos(th)
    v = np.array([x, y, z], dtype=float)
    return v / np.linalg.norm(v)


def run_index_from_path(path: Path) -> int:
    m = re.search(r"_run(\d+)_", path.name)
    return int(m.group(1)) if m else 0


def load_runs(data_dir: Path, tag: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load all runs for a tag; offset EventIDs; attach run_index column."""
    summary_files = sorted(data_dir.glob(f"{tag}_run*_nt_sipm_summary.csv"))
    event_files = sorted(data_dir.glob(f"{tag}_run*_nt_events.csv"))
    if summary_files:
        return _merge_run_pairs(summary_files, event_files, tag)

    # Phase 2 cloud shards: direction_scan_E{E}_t{θ}_p{φ}_s{shard}_run0_nt_*.csv
    if tag == "direction_scan":
        phase2 = load_phase2_shards(data_dir)
        if phase2 is not None:
            return phase2

    raise FileNotFoundError(f"No sipm_summary for tag {tag} in {data_dir}")


def load_phase2_shards(
    data_dir: Path, energy_gev: float | None = None
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Load Phase 2 sharded direction_scan CSVs; truth from events DirX/Y/Z."""
    if energy_gev is not None:
        glob_pat = f"direction_scan_E{energy_gev:g}_*_run*_nt_sipm_summary.csv"
    else:
        glob_pat = "direction_scan_E*_t*_p*_s*_run*_nt_sipm_summary.csv"
    summary_files = sorted(data_dir.glob(glob_pat))
    if not summary_files:
        return None

    event_files = []
    for sf in summary_files:
        ef = Path(str(sf).replace("_nt_sipm_summary.csv", "_nt_events.csv"))
        if not ef.is_file():
            raise FileNotFoundError(f"Missing events CSV for {sf.name}")
        event_files.append(ef)

    summaries, events = [], []
    offset = 0
    for sf, ef in zip(summary_files, event_files):
        s = parse_geant4_csv(sf)
        e = parse_geant4_csv(ef)
        if offset > 0:
            s = s.copy()
            e = e.copy()
            s["EventID"] = s["EventID"] + offset
            e["EventID"] = e["EventID"] + offset
        summaries.append(s)
        events.append(e)
        offset += int(max(s["EventID"].max(), e["EventID"].max(), 0)) + 1

    return pd.concat(summaries, ignore_index=True), pd.concat(events, ignore_index=True)


def _merge_run_pairs(
    summary_files: list[Path], event_files: list[Path], tag: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries, events = [], []
    offset = 0
    for sf, ef in zip(summary_files, event_files):
        s = parse_geant4_csv(sf)
        e = parse_geant4_csv(ef)
        run_idx = run_index_from_path(sf)
        s = s.copy()
        e = e.copy()
        e["RunIndex"] = run_idx
        if offset > 0:
            s["EventID"] = s["EventID"] + offset
            e["EventID"] = e["EventID"] + offset
        if tag == "direction_scan":
            truth = truth_direction_direction_scan(run_idx)
            e["DirX"] = truth[0]
            e["DirY"] = truth[1]
            e["DirZ"] = truth[2]
        summaries.append(s)
        events.append(e)
        offset += int(max(s["EventID"].max(), e["EventID"].max(), 0)) + 1

    return pd.concat(summaries, ignore_index=True), pd.concat(events, ignore_index=True)
