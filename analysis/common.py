"""Shared utilities for Geant4 CSV ntuple analysis."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
FIGURES_DIR = ROOT / "figures"
OUTPUT_DIR = ROOT / "analysis" / "output"

GRID = 18
SIPM_PER_FACE = GRID * GRID
TOTAL_SIPMS = SIPM_PER_FACE * 6
PDE = 0.25
FILTER_CUTOFF_NM = 175.0


def ensure_dirs() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_geant4_csv(path: Path) -> pd.DataFrame:
    """Load Geant4 tools::wcsv ntuple (comment header + comma rows)."""
    columns: list[str] = []
    rows: list[list[str]] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#column"):
                col = line.split(maxsplit=2)[-1]
                columns.append(col.split()[-1] if " " in col else col)
            elif line.startswith("#"):
                continue
            else:
                rows.append(line.split(","))

    if not columns:
        raise ValueError(f"No columns found in {path}")

    df = pd.DataFrame(rows, columns=columns)
    for col in df.columns:
        if col == "Process":
            continue
        if re.search(r"int", col, re.I):
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        else:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def find_run_files(tag: str | None = None, run_id: int = 0) -> dict[str, Path]:
    """Locate ntuple CSV files under data/ or build/data/."""
    patterns = []
    if tag:
        patterns.append(f"{tag}_run{run_id}_nt_*.csv")
    patterns.append(f"*_run{run_id}_nt_*.csv")

    search_roots = [DATA_DIR, ROOT / "build" / "data"]
    for root in search_roots:
        if not root.exists():
            continue
        for pattern in patterns:
            hits = sorted(root.glob(pattern))
            if hits:
                out: dict[str, Path] = {}
                for p in hits:
                    key = p.stem.split("_nt_")[-1]
                    out[key] = p
                return out
    raise FileNotFoundError(f"No run CSV files found for tag={tag!r}, run_id={run_id}")


def face_index(sipm_id: int) -> int:
    return int(sipm_id) // SIPM_PER_FACE


def face_local_id(sipm_id: int) -> int:
    return int(sipm_id) % SIPM_PER_FACE


def face_grid_counts(summary: pd.DataFrame, event_id: int | None = None) -> dict[int, pd.DataFrame]:
    """Return {face_index: 18x18 grid} of photon counts."""
    df = summary if event_id is None else summary[summary["EventID"] == event_id]
    grouped = df.groupby(["SiPMID"], as_index=False)["PhotonCount"].sum()
    faces: dict[int, pd.DataFrame] = {}
    for fid in range(6):
        grid = grouped[grouped["SiPMID"].between(fid * SIPM_PER_FACE, (fid + 1) * SIPM_PER_FACE - 1)]
        arr = [[0] * GRID for _ in range(GRID)]
        for _, row in grid.iterrows():
            local = face_local_id(int(row["SiPMID"]))
            r, c = divmod(local, GRID)
            arr[r][c] = int(row["PhotonCount"])
        faces[fid] = pd.DataFrame(arr)
    return faces


def peak_factor(summary: pd.DataFrame, event_id: int | None = None) -> float:
    df = summary if event_id is None else summary[summary["EventID"] == event_id]
    counts = df.groupby("SiPMID")["PhotonCount"].sum()
    if counts.empty:
        return float("nan")
    return float(counts.max() / counts.mean())


def load_meta(run_dir: Path) -> dict:
    meta_path = run_dir / "meta.json"
    if meta_path.exists():
        return json.loads(meta_path.read_text())
    return {}


def save_summary(name: str, payload: dict) -> Path:
    ensure_dirs()
    path = OUTPUT_DIR / f"{name}.json"
    path.write_text(json.dumps(payload, indent=2))
    return path
