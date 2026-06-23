#!/usr/bin/env python3
"""Convolve simulated Cherenkov wavelength spectrum with filter transmission."""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    FIGURES_DIR,
    FILTER_CUTOFF_NM,
    OUTPUT_DIR,
    PDE,
    ensure_dirs,
    find_run_files,
    parse_geant4_csv,
    save_summary,
)

ROOT = Path(__file__).resolve().parents[1]


def load_filter_transmission() -> pd.DataFrame:
    """Load filter T(λ) from Filter-Design module or cached CSV."""
    cached = OUTPUT_DIR / "filter_transmission.csv"
    if cached.exists():
        return pd.read_csv(cached)

    spec = importlib.util.spec_from_file_location(
        "filter_design", ROOT / "analysis" / "Filter-Design.py"
    )
    if spec and spec.loader:
        mod = importlib.util.module_from_spec(spec)
        sys.modules["filter_design"] = mod
        spec.loader.exec_module(mod)
        if hasattr(mod, "main"):
            mod.main()

    if cached.exists():
        return pd.read_csv(cached)
    raise FileNotFoundError("Run Filter-Design.py first to produce filter_transmission.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=None)
    parser.add_argument("--run-id", type=int, default=0)
    args = parser.parse_args()

    ensure_dirs()
    files = find_run_files(args.tag, args.run_id)
    photons = parse_geant4_csv(files["photons"])
    cherenkov = photons[photons["Process"] == "Cerenkov"]["Wavelength_nm"].astype(float)

    hist, edges = np.histogram(cherenkov, bins=80, range=(160, 600))
    centers = 0.5 * (edges[:-1] + edges[1:])
    spec = pd.DataFrame({"wavelength_nm": centers, "weight": hist / hist.sum()})

    filt = load_filter_transmission()
    wl_col = "wavelength_nm" if "wavelength_nm" in filt.columns else filt.columns[0]
    t_col = "transmission" if "transmission" in filt.columns else filt.columns[1]
    t_interp = np.interp(spec["wavelength_nm"], filt[wl_col], filt[t_col], left=0, right=0)

    total = spec["weight"].sum()
    passed = (spec["weight"] * t_interp).sum()
    cutoff_frac = spec.loc[spec["wavelength_nm"] < FILTER_CUTOFF_NM, "weight"].sum() / total
    filter_fraction = passed / total if total else 0.0

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(spec["wavelength_nm"], spec["weight"], width=4, alpha=0.6, label="Cherenkov spectrum")
    ax.plot(filt[wl_col], filt[t_col] * spec["weight"].max(), "r-", lw=2, label="Filter T(λ)")
    ax.axvline(FILTER_CUTOFF_NM, color="k", ls=":", label="175 nm scint cutoff")
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("Normalized weight")
    ax.set_title(f"Filter folding: {100 * filter_fraction:.1f}% transmission")
    ax.legend()
    fig.tight_layout()
    out = FIGURES_DIR / "filter_folding.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    events = parse_geant4_csv(files["events"])
    e_mean = float(events["E_primary_GeV"].mean())
    n_mean = float(events["N_Cherenkov"].mean())
    photons_per_gev = n_mean / e_mean if e_mean else float("nan")
    filtered_per_gev = photons_per_gev * filter_fraction * PDE / PDE  # PDE already in sim hits

    save_summary(
        "filter_folding",
        {
            "filter_fraction": float(filter_fraction),
            "scint_fraction_below_175nm": float(cutoff_frac),
            "photons_per_GeV_total": float(photons_per_gev),
            "photons_per_GeV_filtered": float(filtered_per_gev),
            "figure": str(out),
        },
    )
    print(f"Filter fraction = {100 * filter_fraction:.1f}%; filtered yield ≈ {filtered_per_gev:.0f} ph/GeV")


if __name__ == "__main__":
    main()
