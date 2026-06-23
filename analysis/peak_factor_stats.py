#!/usr/bin/env python3
"""Peak factor averaged over many events at fixed energy."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, parse_geant4_csv, peak_factor, save_summary
from recon.io import load_runs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="peak_factor")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    args = parser.parse_args()

    ensure_dirs()
    summary, events = load_runs(args.data_dir, args.tag)
    e0 = float(events["E_primary_GeV"].iloc[0])
    pfs = []
    for eid in summary["EventID"].unique():
        pf = peak_factor(summary, int(eid))
        if np.isfinite(pf):
            pfs.append(pf)

    if not pfs:
        raise RuntimeError(f"No peak factors for tag {args.tag}")

    arr = np.array(pfs)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(arr, bins=min(20, len(arr)), color="steelblue", edgecolor="k")
    ax.axvline(arr.mean(), color="r", ls="--", label=f"mean = {arr.mean():.2f}×")
    ax.set_xlabel("Peak factor (max/mean SiPM occupancy)")
    ax.set_ylabel("Events")
    ax.set_title(f"Peak factor @ {e0:.0f} GeV (n={len(arr)})")
    ax.legend()
    fig.tight_layout()
    out = FIGURES_DIR / "peak_factor_hist.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    save_summary(
        "peak_factor_stats",
        {
            "primary_energy_GeV": e0,
            "n_events": int(len(arr)),
            "peak_factor_mean": float(arr.mean()),
            "peak_factor_std": float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
            "figure": str(out),
        },
    )
    print(f"Peak factor mean = {arr.mean():.2f}±{arr.std(ddof=1):.2f} @ {e0:.0f} GeV")


if __name__ == "__main__":
    main()
