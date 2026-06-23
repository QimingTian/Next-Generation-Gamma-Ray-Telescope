#!/usr/bin/env python3
"""Gamma vs proton separation using shower shape observables."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, save_summary
from recon.io import load_runs
from recon.route_a import face_totals


def shape_features(summary: pd.DataFrame, events: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for eid in events["EventID"].unique():
        ev = events[events["EventID"] == eid].iloc[0]
        ftot = face_totals(summary, int(eid))
        total = ftot.sum()
        if total <= 0:
            continue
        frac = ftot / total
        rows.append(
            {
                "EventID": int(eid),
                "E_primary_GeV": float(ev["E_primary_GeV"]),
                "label": label,
                "face_anisotropy": float(frac.max() / (frac.mean() + 1e-9)),
                "n_photons": float(summary[summary["EventID"] == eid]["PhotonCount"].sum()),
                "L90_mm": float(ev.get("L90_mm", np.nan)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma-tag", default="direction_scan")
    parser.add_argument("--proton-tag", default="proton_background")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    args = parser.parse_args()

    ensure_dirs()

    def load_tag(tag: str, label: str) -> pd.DataFrame:
        try:
            summary, events = load_runs(args.data_dir, tag)
        except FileNotFoundError:
            return pd.DataFrame()
        feats = shape_features(summary, events, label)
        return feats

    gamma = load_tag(args.gamma_tag, "gamma")
    proton = load_tag(args.proton_tag, "proton")
    if gamma.empty and proton.empty:
        raise FileNotFoundError("Need gamma and/or proton campaign data")

    combined = pd.concat([gamma, proton], ignore_index=True)
    combined.to_csv(OUTPUT_DIR / "background_features.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    for label, grp in combined.groupby("label"):
        ax.scatter(
            grp["face_anisotropy"],
            grp["n_photons"],
            alpha=0.5,
            s=20,
            label=label,
        )
    ax.set_xlabel("Face anisotropy (max/mean fraction)")
    ax.set_ylabel("Detected photons")
    ax.set_title("Gamma vs proton shower shape")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "background_rejection.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    rejection = {}
    if not gamma.empty and not proton.empty:
        g_med = gamma["face_anisotropy"].median()
        p_med = proton["face_anisotropy"].median()
        rejection = {
            "gamma_median_anisotropy": float(g_med),
            "proton_median_anisotropy": float(p_med),
            "separation_ratio": float(p_med / g_med) if g_med else float("nan"),
        }

    save_summary("background_rejection", {"figure": str(out), **rejection})
    print(f"Wrote {out}; rejection stats: {rejection}")


if __name__ == "__main__":
    main()
