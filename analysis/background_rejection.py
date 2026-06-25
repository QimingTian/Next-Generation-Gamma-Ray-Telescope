#!/usr/bin/env python3
"""Gamma vs proton separation using shower shape + ACD observables."""
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
                "acd_veto": int(ev.get("ACD_veto", 0)),
                "acd_n_tiles": int(ev.get("ACD_n_tiles", 0)),
                "acd_edep_MeV": float(ev.get("ACD_edep_MeV", 0.0)),
            }
        )
    return pd.DataFrame(rows)


def roc_curve(y_true: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """y_true=1 for signal (gamma keep). Higher score → more proton-like → reject."""
    order = np.argsort(scores)
    scores = scores[order]
    y_true = y_true[order]
    n_sig = (y_true == 1).sum()
    n_bkg = (y_true == 0).sum()
    tpr, fpr = [1.0], [0.0]
    sig_kept = n_sig
    bkg_kept = n_bkg
    for i in range(len(scores)):
        if y_true[i] == 1:
            sig_kept -= 1
        else:
            bkg_kept -= 1
        tpr.append(sig_kept / n_sig if n_sig else 0.0)
        fpr.append(bkg_kept / n_bkg if n_bkg else 0.0)
    return np.array(fpr), np.array(tpr)


def combined_reject(
    anisotropy: np.ndarray, acd_edep: np.ndarray, acd_n_tiles: np.ndarray, tau: float, edep_cut: float
) -> np.ndarray:
    """Reject proton-like LXe pattern OR strong ACD charged deposit."""
    return (anisotropy > tau) | ((acd_edep > edep_cut) & (acd_n_tiles >= 3))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma-tag", default="direction_scan")
    parser.add_argument("--proton-tag", default="proton_background")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--combined", action="store_true", help="Plot shape+ACD combined ROC")
    args = parser.parse_args()

    ensure_dirs()

    def load_tag(tag: str, label: str) -> pd.DataFrame:
        try:
            summary, events = load_runs(args.data_dir, tag)
        except FileNotFoundError:
            return pd.DataFrame()
        return shape_features(summary, events, label)

    gamma = load_tag(args.gamma_tag, "gamma")
    proton = load_tag(args.proton_tag, "proton")
    proton_supp = load_tag("proton_supplement", "proton")
    if not proton_supp.empty:
        proton = pd.concat([proton, proton_supp], ignore_index=True)

    # ACD campaign: merge all acd_gamma_* / acd_proton_* if primary tags missing
    if gamma.empty:
        parts = []
        for tag in sorted({p.name.split("_run")[0] for p in args.data_dir.glob("acd_gamma_*_run*_nt_events.csv")}):
            parts.append(load_tag(tag, "gamma"))
        if parts:
            gamma = pd.concat(parts, ignore_index=True)
    if proton.empty:
        parts = []
        for tag in sorted({p.name.split("_run")[0] for p in args.data_dir.glob("acd_proton_*_run*_nt_events.csv")}):
            parts.append(load_tag(tag, "proton"))
        if parts:
            proton = pd.concat(parts, ignore_index=True)

    # Smoke-test fallback (run0=γ, run1/2=p)
    if gamma.empty or (proton.empty and args.combined):
        from common import parse_geant4_csv

        def load_smoke_run(run_idx: int, label: str) -> pd.DataFrame:
            efs = sorted(args.data_dir.glob(f"acd_smoke_run{run_idx}_nt_events.csv"))
            sfs = sorted(args.data_dir.glob(f"acd_smoke_run{run_idx}_nt_sipm_summary.csv"))
            if not efs or not sfs:
                return pd.DataFrame()
            return shape_features(parse_geant4_csv(sfs[0]), parse_geant4_csv(efs[0]), label)

        if gamma.empty:
            gamma = load_smoke_run(0, "gamma")
        if proton.empty and args.combined:
            proton = pd.concat(
                [load_smoke_run(1, "proton"), load_smoke_run(2, "proton")], ignore_index=True
            )

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

    rejection: dict = {"figure": str(out)}
    if not gamma.empty and not proton.empty:
        g_med = gamma["face_anisotropy"].median()
        p_med = proton["face_anisotropy"].median()
        rejection.update(
            {
                "gamma_median_anisotropy": float(g_med),
                "proton_median_anisotropy": float(p_med),
                "separation_ratio": float(p_med / g_med) if g_med else float("nan"),
            }
        )
        if "acd_veto" in combined.columns:
            g_pass = float((gamma["acd_veto"] == 0).mean())
            p_rej = float((proton["acd_veto"] == 1).mean())
            rejection["acd_only_gamma_pass"] = g_pass
            rejection["acd_only_proton_reject"] = p_rej

        if args.combined and "acd_edep_MeV" in combined.columns:
            y = np.concatenate([np.ones(len(gamma)), np.zeros(len(proton))])
            shape_score = np.concatenate([gamma["face_anisotropy"].values, proton["face_anisotropy"].values])
            acd_edep = np.concatenate([gamma["acd_edep_MeV"].values, proton["acd_edep_MeV"].values])
            acd_tiles = np.concatenate([gamma["acd_n_tiles"].values, proton["acd_n_tiles"].values])

            fpr_s, tpr_s = roc_curve(y, shape_score)
            combined_score = shape_score + 0.002 * acd_edep + 0.05 * acd_tiles
            fpr_c, tpr_c = roc_curve(y, combined_score)

            fig2, ax2 = plt.subplots(figsize=(6, 5))
            ax2.plot(fpr_s, tpr_s, label="LXe shape only")
            ax2.plot(fpr_c, tpr_c, label="Shape + ACD edep/tiles")
            ax2.set_xlabel("Proton acceptance (background leak)")
            ax2.set_ylabel("Gamma acceptance (signal efficiency)")
            ax2.set_title("Combined veto ROC")
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            fig2.tight_layout()
            roc_out = FIGURES_DIR / "acd_combined_roc.png"
            fig2.savefig(roc_out, dpi=150)
            plt.close(fig2)
            rejection["combined_roc_figure"] = str(roc_out)

            # Report best combined cut: edep=10 MeV + anisotropy sweep
            best = {"gamma_eff": 0.0, "proton_reject": 0.0, "tau": 0.0, "edep_cut_MeV": 10.0}
            for tau in np.linspace(gamma["face_anisotropy"].quantile(0.1), proton["face_anisotropy"].quantile(0.9), 40):
                for edep_cut in (5.0, 10.0, 20.0, 50.0):
                    g_keep = ~combined_reject(
                        gamma["face_anisotropy"].values,
                        gamma["acd_edep_MeV"].values,
                        gamma["acd_n_tiles"].values,
                        tau,
                        edep_cut,
                    )
                    p_rej = combined_reject(
                        proton["face_anisotropy"].values,
                        proton["acd_edep_MeV"].values,
                        proton["acd_n_tiles"].values,
                        tau,
                        edep_cut,
                    )
                    g_eff = float(g_keep.mean())
                    p_rej_frac = float(p_rej.mean())
                    if p_rej_frac > best["proton_reject"] and g_eff >= 0.5 * best.get("gamma_eff", 0):
                        if g_eff + p_rej_frac > best["gamma_eff"] + best["proton_reject"]:
                            best = {
                                "gamma_eff": g_eff,
                                "proton_reject": p_rej_frac,
                                "tau": float(tau),
                                "edep_cut_MeV": edep_cut,
                            }
            rejection["combined_best_cut"] = best

    save_summary("background_rejection", rejection)
    print(f"Wrote {out}; rejection stats: {rejection}")


if __name__ == "__main__":
    main()
