#!/usr/bin/env python3
"""A_eff vs energy and incidence angle from gap campaigns."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, parse_geant4_csv, save_summary

THROW_AREA_M2 = 1.6 * 1.6
EDEP_THRESHOLD_MEV = 1.0


def aeff_from_tag(data_dir: Path, tag_prefix: str) -> dict | None:
    files = sorted(data_dir.glob(f"{tag_prefix}*_run*_nt_events.csv"))
    if not files:
        return None
    events = pd.concat([parse_geant4_csv(f) for f in files], ignore_index=True)
    n_thrown = len(events)
    n_int = int((events["Edep_MeV"] > EDEP_THRESHOLD_MEV).sum())
    return {
        "tag": tag_prefix,
        "n_thrown": n_thrown,
        "n_interacting": n_int,
        "A_eff_m2": float((n_int / n_thrown) * THROW_AREA_M2) if n_thrown else float("nan"),
        "E_primary_GeV": float(events["E_primary_GeV"].iloc[0]) if n_thrown else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    args = parser.parse_args()

    ensure_dirs()
    rows: list[dict] = []

    # Existing vertical 200 GeV campaign
    base = aeff_from_tag(args.data_dir, "plane_source")
    if base:
        base["theta_deg"] = 0.0
        base["source"] = "plane_source"
        rows.append(base)

    for e in (100, 300, 500, 1000):
        r = aeff_from_tag(args.data_dir, f"aeff_energy_{e:g}GeV")
        if r:
            r["theta_deg"] = 0.0
            r["source"] = "aeff_energy"
            rows.append(r)

    for theta in (15, 30, 45):
        r = aeff_from_tag(args.data_dir, f"aeff_theta_{theta:g}deg")
        if r:
            r["theta_deg"] = float(theta)
            r["source"] = "aeff_theta"
            rows.append(r)

    if not rows:
        raise FileNotFoundError("No A_eff campaign data found")

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / "effective_area_scan.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    vert = df[df["theta_deg"] == 0].sort_values("E_primary_GeV")
    if not vert.empty:
        axes[0].plot(vert["E_primary_GeV"], vert["A_eff_m2"], "o-", color="steelblue")
        axes[0].set_xlabel("Energy [GeV]")
        axes[0].set_ylabel(r"$A_{\rm eff}$ [m$^2$]")
        axes[0].set_title(r"$A_{\rm eff}(E)$, vertical")
        axes[0].grid(True, alpha=0.3)

    tilt = df[(df["source"] == "aeff_theta") | ((df["source"] == "plane_source") & (df["theta_deg"] == 0))]
    tilt = tilt.sort_values("theta_deg")
    if len(tilt) > 1:
        axes[1].plot(tilt["theta_deg"], tilt["A_eff_m2"], "s-", color="coral")
        axes[1].set_xlabel(r"Incidence angle $\theta$ [deg]")
        axes[1].set_ylabel(r"$A_{\rm eff}$ [m$^2$]")
        axes[1].set_title(r"$A_{\rm eff}(\theta)$ @ 200 GeV")
        axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    out = FIGURES_DIR / "effective_area_scan.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    save_summary(
        "effective_area_scan",
        {"n_points": len(df), "figure": str(out), "points": df.to_dict(orient="records")},
    )
    print(df.to_string(index=False))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
