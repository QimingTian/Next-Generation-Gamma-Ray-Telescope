#!/usr/bin/env python3
"""Compare energy resolution with filter on vs off."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, parse_geant4_csv, save_summary


def load_events(data_dir: Path, tag: str) -> pd.DataFrame:
    files = sorted(data_dir.glob(f"{tag}_run*_nt_events.csv"))
    if not files:
        raise FileNotFoundError(tag)
    events = pd.concat([parse_geant4_csv(f) for f in files], ignore_index=True)
    events["N_total"] = events["N_Cherenkov"] + events["N_Scint"]
    return events


def resolution(events: pd.DataFrame) -> dict:
    e0 = float(events["E_primary_GeV"].iloc[0])
    slope = events["N_total"].mean() / e0
    e_reco = events["N_total"] / slope
    sigma = e_reco.std(ddof=1)
    mean = e_reco.mean()
    rel = sigma / mean if mean else float("nan")
    poisson = 1.0 / np.sqrt(events["N_total"].mean()) if events["N_total"].mean() else float("nan")
    return {
        "primary_energy_GeV": e0,
        "n_events": int(len(events)),
        "calibration_slope": float(slope),
        "sigma_E_over_E_percent": float(100 * rel),
        "poisson_limit_percent": float(100 * poisson),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    args = parser.parse_args()

    ensure_dirs()
    pairs = [
        ("energy_resolution_100_on", "energy_resolution_100_off", 100),
        ("energy_resolution_215_on", "energy_resolution_215_off", 215),
    ]
    rows = []
    for on_tag, off_tag, energy in pairs:
        try:
            on = load_events(args.data_dir, on_tag)
            off = load_events(args.data_dir, off_tag)
        except FileNotFoundError:
            continue
        r_on = resolution(on)
        r_off = resolution(off)
        filter_frac = r_on["calibration_slope"] / r_off["calibration_slope"] if r_off["calibration_slope"] else float("nan")
        rows.append({"energy_GeV": energy, **{f"{k}_on": v for k, v in r_on.items()}, **{f"{k}_off": v for k, v in r_off.items()}, "filter_fraction": float(filter_frac)})

    if not rows:
        print("No filter on/off pairs found; run overnight energy resolution campaign first.")
        return

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / "energy_resolution_filter_compare.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.arange(len(df))
    width = 0.35
    ax.bar(x - width / 2, df["sigma_E_over_E_percent_on"], width, label="Filter ON")
    ax.bar(x + width / 2, df["sigma_E_over_E_percent_off"], width, label="Filter OFF")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{e:.0f} GeV" for e in df["energy_GeV"]])
    ax.set_ylabel(r"$\sigma_E/E$ [%]")
    ax.set_title("Energy resolution: filter ON vs OFF (photon counting)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "energy_resolution_filter_compare.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    save_summary("energy_resolution_filter_compare", {"rows": rows, "figure": str(out)})
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
