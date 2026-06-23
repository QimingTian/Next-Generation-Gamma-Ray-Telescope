#!/usr/bin/env python3
"""Effective area from plane-source Monte Carlo."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, parse_geant4_csv, save_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="plane_source")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--throw-area-m2", type=float, default=(1.6 * 1.6))
    parser.add_argument("--edep-threshold-mev", type=float, default=1.0)
    args = parser.parse_args()

    ensure_dirs()
    files = sorted(args.data_dir.glob(f"{args.tag}_run*_nt_events.csv"))
    if not files:
        raise FileNotFoundError(f"No events for tag {args.tag}")

    events = parse_geant4_csv(files[0])
    n_thrown = len(events)
    n_interact = int((events["Edep_MeV"] > args.edep_threshold_mev).sum())
    a_eff = (n_interact / n_thrown) * args.throw_area_m2 if n_thrown else float("nan")

    energy = float(events["E_primary_GeV"].iloc[0]) if n_thrown else float("nan")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["Thrown", "Interacting"], [n_thrown, n_interact], color=["gray", "steelblue"])
    ax.set_ylabel("Events")
    ax.set_title(f"A_eff = {a_eff:.3f} m² @ {energy:.0f} GeV")
    fig.tight_layout()
    out = FIGURES_DIR / "effective_area.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    save_summary(
        "effective_area",
        {
            "primary_energy_GeV": energy,
            "n_thrown": n_thrown,
            "n_interacting": n_interact,
            "throw_area_m2": args.throw_area_m2,
            "A_eff_m2": float(a_eff),
            "figure": str(out),
        },
    )
    print(f"A_eff = {a_eff:.4f} m² ({n_interact}/{n_thrown})")


if __name__ == "__main__":
    main()
