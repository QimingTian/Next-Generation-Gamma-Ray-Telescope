#!/usr/bin/env python3
"""Self-veto fraction vs energy/zenith from ACD-on gamma runs."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, parse_geant4_csv, save_summary


def tag_meta(tag: str) -> dict:
    m_e = re.search(r"E(\d+(?:\.\d+)?)GeV", tag)
    m_t = re.search(r"theta(\d+(?:\.\d+)?)", tag)
    return {
        "energy_GeV": float(m_e.group(1)) if m_e else np.nan,
        "theta_deg": float(m_t.group(1)) if m_t else 0.0,
    }


def load_acd_gamma(data_dir: Path) -> pd.DataFrame:
    rows = []
    for ef in sorted(data_dir.glob("acd_gamma_*_run*_nt_events.csv")):
        tag = ef.name.split("_run")[0]
        meta = tag_meta(tag)
        df = parse_geant4_csv(ef)
        for _, ev in df.iterrows():
            rows.append(
                {
                    "tag": tag,
                    "energy_GeV": meta["energy_GeV"],
                    "theta_deg": meta["theta_deg"],
                    "acd_veto": int(ev.get("ACD_veto", 0)),
                    "acd_n_tiles": float(ev.get("ACD_n_tiles", 0)),
                    "acd_edep_MeV": float(ev.get("ACD_edep_MeV", 0)),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    args = parser.parse_args()

    ensure_dirs()
    df = load_acd_gamma(args.data_dir)
    if df.empty:
        # Fall back to smoke-test naming
        for ef in sorted(args.data_dir.glob("acd_smoke_run0_nt_events.csv")):
            smoke = parse_geant4_csv(ef)
            smoke["energy_GeV"] = 100.0
            smoke["theta_deg"] = 0.0
            smoke["tag"] = "acd_smoke"
            df = smoke.rename(
                columns={"ACD_veto": "acd_veto", "ACD_n_tiles": "acd_n_tiles", "ACD_edep_MeV": "acd_edep_MeV"}
            )
            break
    if df.empty:
        raise FileNotFoundError("No ACD gamma CSVs found")

    by_e = (
        df.groupby("energy_GeV", dropna=False)
        .agg(veto_frac=("acd_veto", "mean"), n=("acd_veto", "count"), median_edep=("acd_edep_MeV", "median"))
        .reset_index()
    )
    by_theta = (
        df[df["energy_GeV"] == 100]
        .groupby("theta_deg")
        .agg(veto_frac=("acd_veto", "mean"), n=("acd_veto", "count"))
        .reset_index()
    )

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    if not by_e.empty and by_e["energy_GeV"].notna().any():
        axes[0].plot(by_e["energy_GeV"], by_e["veto_frac"], "o-", label="veto fraction")
        axes[0].set_xscale("log")
        axes[0].set_xlabel("Energy (GeV)")
        axes[0].set_ylabel("ACD veto fraction (γ)")
        axes[0].set_title("Gamma self-veto vs E")
        axes[0].grid(True, alpha=0.3)
    if not by_theta.empty:
        axes[1].plot(by_theta["theta_deg"], by_theta["veto_frac"], "s-", color="C1")
        axes[1].set_xlabel("Zenith θ (deg)")
        axes[1].set_ylabel("ACD veto fraction (γ @ 100 GeV)")
        axes[1].set_title("Gamma self-veto vs θ")
        axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    out = FIGURES_DIR / "acd_self_veto_scan.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    summary = {
        "figure": str(out),
        "by_energy": by_e.to_dict(orient="records"),
        "by_theta_100GeV": by_theta.to_dict(orient="records"),
    }
    save_summary("acd_self_veto_scan", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
