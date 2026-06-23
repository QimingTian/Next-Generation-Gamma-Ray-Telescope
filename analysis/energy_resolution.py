#!/usr/bin/env python3
"""Energy resolution from photon-count fluctuations at fixed energy."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import FIGURES_DIR, OUTPUT_DIR, ensure_dirs, find_run_files, parse_geant4_csv, save_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="energy_resolution")
    parser.add_argument("--run-id", type=int, default=0)
    parser.add_argument("--calibration-slope", type=float, default=None,
                        help="photons/GeV; auto from mean if omitted")
    args = parser.parse_args()

    ensure_dirs()
    files = find_run_files(args.tag, args.run_id)
    events = parse_geant4_csv(files["events"])
    events["N_total"] = events["N_Cherenkov"] + events["N_Scint"]

    e0 = float(events["E_primary_GeV"].iloc[0])
    slope = args.calibration_slope or (events["N_total"].mean() / e0)
    events["E_reco"] = events["N_total"] / slope

    sigma_e = events["E_reco"].std(ddof=1)
    mean_e = events["E_reco"].mean()
    rel_res = sigma_e / mean_e if mean_e else float("nan")
    poisson_limit = 1.0 / np.sqrt(events["N_total"].mean()) if events["N_total"].mean() else float("nan")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].hist(events["E_reco"], bins=min(20, len(events)), color="steelblue", edgecolor="k")
    axes[0].axvline(e0, color="r", ls="--", label=f"true E = {e0:.0f} GeV")
    axes[0].set_xlabel("Reconstructed energy [GeV]")
    axes[0].set_ylabel("Events")
    axes[0].set_title(f"σ_E/E = {100 * rel_res:.2f}%")
    axes[0].legend()

    axes[1].scatter(events["N_total"], events["E_reco"], alpha=0.7)
    axes[1].set_xlabel("Photon count")
    axes[1].set_ylabel("Reconstructed E [GeV]")
    axes[1].set_title(f"Poisson limit ≈ {100 * poisson_limit:.2f}%")
    fig.tight_layout()
    out = FIGURES_DIR / "energy_resolution.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    events[["EventID", "N_total", "E_reco"]].to_csv(OUTPUT_DIR / "energy_resolution_events.csv", index=False)
    summary = {
        "primary_energy_GeV": e0,
        "n_events": int(len(events)),
        "calibration_slope": float(slope),
        "sigma_E_over_E": float(rel_res),
        "sigma_E_over_E_percent": float(100 * rel_res),
        "poisson_limit_percent": float(100 * poisson_limit),
        "figure": str(out),
    }
    save_summary("energy_resolution", summary)
    print(f"σ_E/E = {100 * rel_res:.2f}% @ {e0:.0f} GeV; wrote {out}")


if __name__ == "__main__":
    main()
