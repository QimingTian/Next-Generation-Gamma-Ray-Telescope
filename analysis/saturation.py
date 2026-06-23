#!/usr/bin/env python3
"""SiPM saturation energy vs array scale (occupancy model)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import FIGURES_DIR, OUTPUT_DIR, PDE, TOTAL_SIPMS, ensure_dirs, save_summary

ROOT = Path(__file__).resolve().parents[1]


def e_sat(n_cells: float, photons_per_gev: float, peak_factor: float, eta: float = 0.2) -> float:
    """Saturation energy [GeV] from paper occupancy model."""
    n_sipm = TOTAL_SIPMS
    n_avg = photons_per_gev / n_sipm
    n_peak = peak_factor * n_avg
    mu_max = -n_cells * np.log(1 - eta)
    return mu_max / (PDE * n_peak)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration-json", type=Path,
                        default=ROOT / "analysis" / "output" / "energy_calibration.json")
    parser.add_argument("--density-json", type=Path,
                        default=ROOT / "analysis" / "output" / "peak_factor_stats.json")
    args = parser.parse_args()

    ensure_dirs()

    cal = json.loads(args.calibration_json.read_text()) if args.calibration_json.exists() else {}
    dens = json.loads(args.density_json.read_text()) if args.density_json.exists() else {}
    photons_per_gev = cal.get("photons_per_GeV_filter_on", cal.get("photons_per_GeV_total", 9800.0))
    peak_factor = dens.get("peak_factor_mean", dens.get("peak_factor", 2.2))

    cell_densities = [14400, 100000, 440000]
    array_scales = [(18, 18, 6), (40, 40, 6)]
    labels = ["18×18×6 baseline", "40×40×6 scaled"]

    fig, ax = plt.subplots(figsize=(9, 6))
    e_sat_baseline = {}
    for label, cells in zip(["14.4k", "100k", "440k"], cell_densities):
        e = e_sat(cells, photons_per_gev, peak_factor) / 1000.0
        e_sat_baseline[label] = e
        ax.axhline(e, ls="--" if label != "14.4k" else "-", label=f"{label} cells: {e:.2f} TeV")

    ax.set_xlabel("Array configuration")
    ax.set_ylabel(r"$E_{\rm sat}$ [TeV]")
    ax.set_title(f"Readout saturation (peak factor={peak_factor:.2f}×, {photons_per_gev:.0f} ph/GeV)")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "saturation_energy.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    summary = {
        "photons_per_GeV": photons_per_gev,
        "peak_factor": peak_factor,
        "E_sat_TeV_14400_cells": e_sat_baseline["14.4k"],
        "E_sat_TeV_100k_cells": e_sat_baseline["100k"],
        "E_sat_TeV_440k_cells": e_sat_baseline["440k"],
        "figure": str(out),
    }
    save_summary("saturation", summary)
    print(f"E_sat (14.4k cells) = {e_sat_baseline['14.4k']:.2f} TeV; wrote {out}")


if __name__ == "__main__":
    main()
