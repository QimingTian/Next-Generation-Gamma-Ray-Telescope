#!/usr/bin/env python3
"""N_photon vs E_primary linear calibration with filter on/off comparison."""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

from common import FIGURES_DIR, OUTPUT_DIR, ensure_dirs, parse_geant4_csv, save_summary

ROOT = Path(__file__).resolve().parents[1]


def load_energy_scan(data_dir: Path, tag: str | None = None) -> pd.DataFrame:
    """Aggregate events from energy_scan run files."""
    if tag:
        pattern = str(data_dir / f"{tag}_run*_nt_events.csv")
    else:
        pattern = str(data_dir / "energy_scan_run*_nt_events.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        files = sorted(glob.glob(str(data_dir / "energy_scan_on_run*_nt_events.csv")))
    if not files:
        raise FileNotFoundError(f"No events CSV found for calibration tag={tag}")

    frames = [parse_geant4_csv(Path(f)) for f in files]
    events = pd.concat(frames, ignore_index=True)
    events["N_total"] = events["N_Cherenkov"].astype(float) + events["N_Scint"].astype(float)
    grouped = (
        events.groupby("E_primary_GeV", as_index=False)
        .agg(N_total_mean=("N_total", "mean"), N_total_std=("N_total", "std"), n=("EventID", "count"))
        .sort_values("E_primary_GeV")
    )
    return grouped


def fit_line(grouped: pd.DataFrame) -> dict:
    x = grouped["E_primary_GeV"].values
    y = grouped["N_total_mean"].values
    fit = linregress(x, y)
    return {
        "slope": float(fit.slope),
        "intercept": float(fit.intercept),
        "r2": float(fit.rvalue ** 2),
        "x": x,
        "y": y,
        "yerr": grouped["N_total_std"].values,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--tag", default=None, help="e.g. energy_scan_on; default tries on/off pair")
    args = parser.parse_args()

    ensure_dirs()
    data_dir = args.data_dir

    try:
        on = load_energy_scan(data_dir, "energy_scan_on")
        off = load_energy_scan(data_dir, "energy_scan_off")
        dual = True
    except FileNotFoundError:
        grouped = load_energy_scan(data_dir, args.tag)
        on_fit = fit_line(grouped)
        off_fit = None
        dual = False

    fig, ax = plt.subplots(figsize=(8, 6))
    summary: dict = {}

    if dual:
        on_fit = fit_line(on)
        off_fit = fit_line(off)
        filter_frac = on_fit["slope"] / off_fit["slope"] if off_fit["slope"] else float("nan")

        ax.errorbar(on_fit["x"], on_fit["y"], yerr=on_fit["yerr"], fmt="o", capsize=3, label="Filter ON")
        ax.errorbar(off_fit["x"], off_fit["y"], yerr=off_fit["yerr"], fmt="s", capsize=3, label="Filter OFF")
        xx = np.linspace(max(1, on_fit["x"].min() * 0.8), on_fit["x"].max() * 1.05, 100)
        ax.plot(xx, on_fit["slope"] * xx + on_fit["intercept"], "r-",
                label=f"ON: N = {on_fit['slope']:.0f}×E (R²={on_fit['r2']:.4f})")
        ax.plot(xx, off_fit["slope"] * xx + off_fit["intercept"], "g--",
                label=f"OFF: N = {off_fit['slope']:.0f}×E (R²={off_fit['r2']:.4f})")
        summary = {
            "photons_per_GeV_filter_on": on_fit["slope"],
            "photons_per_GeV_filter_off": off_fit["slope"],
            "filter_fraction_slope_ratio": float(filter_frac),
            "intercept_on": on_fit["intercept"],
            "intercept_off": off_fit["intercept"],
            "r_squared_on": on_fit["r2"],
            "r_squared_off": off_fit["r2"],
        }
        slope = on_fit["slope"]
        intercept = on_fit["intercept"]
        r2 = on_fit["r2"]
    else:
        ax.errorbar(on_fit["x"], on_fit["y"], yerr=on_fit["yerr"], fmt="o", capsize=3, label="simulation")
        xx = np.linspace(max(1, on_fit["x"].min() * 0.8), on_fit["x"].max() * 1.05, 100)
        ax.plot(xx, on_fit["slope"] * xx + on_fit["intercept"], "r-",
                label=f"N = {on_fit['slope']:.0f}×E (R²={on_fit['r2']:.4f})")
        slope, intercept, r2 = on_fit["slope"], on_fit["intercept"], on_fit["r2"]
        summary = {
            "photons_per_GeV_total": slope,
            "intercept": intercept,
            "r_squared": r2,
        }

    ax.set_xlabel("Primary energy [GeV]")
    ax.set_ylabel("Photons detected")
    ax.set_title("Energy calibration: photon yield vs energy")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "energy_calibration.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    summary["figure"] = str(out)
    save_summary("energy_calibration", summary)
    print(f"Calibration slope (filter ON) = {slope:.1f} photons/GeV; wrote {out}")


if __name__ == "__main__":
    main()
