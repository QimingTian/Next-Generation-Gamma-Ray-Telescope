#!/usr/bin/env python3
"""Electromagnetic shower longitudinal length vs primary energy (L90)."""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, parse_geant4_csv, save_summary


def saturation_model(e, linf, c, k, p):
    return linf * (1 - np.exp(-c * np.power(e, k))) ** p


def load_events(data_dir: Path, tag: str | None) -> pd.DataFrame:
    if tag:
        pattern = str(data_dir / f"{tag}_run*_nt_events.csv")
    else:
        pattern = str(data_dir / "energy_scan_*_run*_nt_events.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        pattern = str(data_dir / "*_nt_events.csv")
        files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No events CSV in {data_dir} (tag={tag})")
    return pd.concat([parse_geant4_csv(Path(f)) for f in files], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--tag", default="energy_scan_on",
                        help="Campaign tag prefix (e.g. energy_scan_on)")
    args = parser.parse_args()

    ensure_dirs()
    events = load_events(args.data_dir, args.tag)
    grouped = (
        events.groupby("E_primary_GeV", as_index=False)
        .agg(
            L_mm=("L90_mm", "mean"),
            L_std=("L90_mm", "std"),
            n_trunc=("L90_truncated", "sum"),
        )
        .sort_values("E_primary_GeV")
    )
    fit_group = grouped[grouped["n_trunc"] < grouped["n_trunc"].max()].copy()
    if fit_group.empty:
        fit_group = grouped

    x = fit_group["E_primary_GeV"].values
    y = fit_group["L_mm"].values

    popt, _ = curve_fit(
        saturation_model,
        x,
        y,
        p0=[y.max() * 1.2, 0.01, 0.5, 1.0],
        maxfev=20000,
        bounds=([0, 0, 0, 0.1], [np.inf, np.inf, 3, 5]),
    )
    linf, c, k, p = popt

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(
        grouped["E_primary_GeV"],
        grouped["L_mm"],
        yerr=grouped["L_std"],
        fmt="o",
        capsize=3,
        label="simulation",
    )
    xx = np.linspace(x.min(), x.max() * 1.05, 200)
    ax.plot(xx, saturation_model(xx, *popt), "r-", label="saturation fit")
    ax.set_xlabel("Primary energy [GeV]")
    ax.set_ylabel("Shower L90 [mm]")
    ax.set_title(
        rf"$L(E) = L_\infty (1 - e^{{-c E^{{{k:.2f}}}}})^{{{p:.2f}}}$, "
        rf"$L_\infty={linf:.0f}$ mm"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "shower_length.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    grouped.to_csv(OUTPUT_DIR / "shower_length_points.csv", index=False)
    save_summary(
        "shower_length",
        {
            "L_inf_mm": float(linf),
            "c": float(c),
            "k": float(k),
            "p": float(p),
            "tag": args.tag,
            "figure": str(out),
        },
    )
    print(f"Shower length fit L_inf={linf:.0f} mm; wrote {out}")


if __name__ == "__main__":
    main()
