#!/usr/bin/env python3
"""Angular resolution (sigma68) and reconstruction figures."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES_DIR, OUTPUT_DIR, ROOT, ensure_dirs, save_summary
from recon.io import load_runs
from recon.route_a import opening_angle_deg, reconstruct_route_a, truth_direction_from_event
from recon.route_b import build_feature_dataset, predict_all, train_and_predict


def sigma68(angles_deg: np.ndarray) -> float:
    if len(angles_deg) == 0:
        return float("nan")
    return float(np.percentile(angles_deg, 68))


def plot_recon_cube(ax, vertex_mm, direction, truth, title: str) -> None:
    half = 0.5
    # Cube wireframe in meters for display.
    for a, b in [(0, 1), (1, 2), (0, 2)]:
        pass
    verts = np.array(
        [
            [-half, -half, -half],
            [half, -half, -half],
            [half, half, -half],
            [-half, half, -half],
            [-half, -half, half],
            [half, -half, half],
            [half, half, half],
            [-half, half, half],
        ]
    )
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    for i, j in edges:
        ax.plot3D(*zip(verts[i], verts[j]), color="0.7", lw=0.8)

    v = vertex_mm / 1000.0
    d_reco = direction
    d_truth = truth
    scale = 0.8
    ax.plot3D([v[0], v[0] + scale * d_reco[0]], [v[1], v[1] + scale * d_reco[1]],
              [v[2], v[2] + scale * d_reco[2]], "r-", lw=2, label="reco")
    ax.plot3D([v[0], v[0] + scale * d_truth[0]], [v[1], v[1] + scale * d_truth[1]],
              [v[2], v[2] + scale * d_truth[2]], "b--", lw=2, label="truth")
    ax.scatter([v[0]], [v[1]], [v[2]], c="k", s=20)
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_zlim(-half, half)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(title)
    ax.legend(loc="upper right")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="direction_scan")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--energy", type=float, default=None, help="Phase2: load E=30/100/300 GeV shards only")
    parser.add_argument("--train-route-b", action="store_true")

    ensure_dirs()
    data_dir = args.data_dir
    try:
        if args.energy is not None:
            from recon.io import load_phase2_shards

            loaded = load_phase2_shards(data_dir, args.energy)
            if loaded is None:
                raise FileNotFoundError(f"No Phase2 E={args.energy} GeV shards in {data_dir}")
            summary, events = loaded
        else:
            summary, events = load_runs(data_dir, args.tag)
    except FileNotFoundError:
        summary, events = load_runs(data_dir, "angle_scan")

    rows = []
    for eid in sorted(summary["EventID"].unique()):
        ev = events[events["EventID"] == eid]
        if ev.empty:
            continue
        truth = truth_direction_from_event(ev.iloc[0])
        reco = reconstruct_route_a(summary, int(eid))
        ang = opening_angle_deg(reco["direction"], truth)
        rows.append(
            {
                "EventID": int(eid),
                "E_primary_GeV": float(ev.iloc[0]["E_primary_GeV"]),
                "DirX": float(truth[0]),
                "DirY": float(truth[1]),
                "DirZ": float(truth[2]),
                "theta_deg": float(np.degrees(np.arccos(np.clip(-truth[2], -1, 1)))),
                "delta_deg_route_a": ang,
            }
        )

    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT_DIR / "angular_resolution_events.csv", index=False)

    # sigma68 vs energy and theta bins.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for energy, grp in results.groupby("E_primary_GeV"):
        axes[0].scatter(
            [energy] * len(grp),
            grp["delta_deg_route_a"],
            alpha=0.35,
            s=12,
        )
    e_stats = results.groupby("E_primary_GeV")["delta_deg_route_a"].apply(sigma68)
    axes[0].plot(e_stats.index, e_stats.values, "ro-", lw=2, label=r"$\sigma_{68}(E)$")
    axes[0].set_xlabel("Primary energy [GeV]")
    axes[0].set_ylabel(r"Opening angle [deg]")
    axes[0].set_title("Route A angular resolution")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    theta_bins = [0, 5, 15, 25, 50]
    theta_centers, sigma_vals = [], []
    for lo, hi in zip(theta_bins[:-1], theta_bins[1:]):
        mask = (results["theta_deg"] >= lo) & (results["theta_deg"] < hi)
        if mask.sum() < 3:
            continue
        theta_centers.append(0.5 * (lo + hi))
        sigma_vals.append(sigma68(results.loc[mask, "delta_deg_route_a"].to_numpy()))
    axes[1].plot(theta_centers, sigma_vals, "bs-", lw=2)
    axes[1].set_xlabel(r"Incidence $\theta$ [deg]")
    axes[1].set_ylabel(r"$\sigma_{68}$ [deg]")
    axes[1].set_title("Angular resolution vs incidence angle")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    ang_fig = FIGURES_DIR / "angular_resolution.png"
    fig.savefig(ang_fig, dpi=150)
    plt.close(fig)

    # Example recon cubes for vertical and max-tilt events.
    if len(results) >= 2:
        fig = plt.figure(figsize=(12, 5))
        examples = [
            results.iloc[results["theta_deg"].idxmin()],
            results.iloc[results["theta_deg"].idxmax()],
        ]
        for i, ex in enumerate(examples, start=1):
            eid = int(ex["EventID"])
            reco = reconstruct_route_a(summary, eid)
            truth = np.array([ex["DirX"], ex["DirY"], ex["DirZ"]])
            ax = fig.add_subplot(1, 2, i, projection="3d")
            plot_recon_cube(
                ax,
                reco["vertex_mm"],
                reco["direction"],
                truth,
                f"Event {eid}: Δ={ex['delta_deg_route_a']:.2f}°",
            )
        fig.tight_layout()
        fig.savefig(FIGURES_DIR / "recon_examples.png", dpi=150)
        plt.close(fig)

    route_b_summary = {}
    if args.train_route_b and len(results) >= 30:
        X, y, meta = build_feature_dataset(data_dir, args.tag, energy_gev=args.energy)
        y_test, y_pred, model = train_and_predict(X, y)
        deltas_b = [
            opening_angle_deg(y_pred[i], y_test[i]) for i in range(len(y_test))
        ]
        route_b_summary = {
            "sigma68_route_b_deg": sigma68(np.array(deltas_b)),
            "n_test": int(len(deltas_b)),
        }
        # Full-campaign predictions for comparison.
        y_all = predict_all(model, X)
        deltas_all = [opening_angle_deg(y_all[i], y[i]) for i in range(len(y))]
        results["delta_deg_route_b"] = deltas_all
        results.to_csv(OUTPUT_DIR / "angular_resolution_events.csv", index=False)

    summary_payload = {
        "n_events": int(len(results)),
        "sigma68_route_a_overall_deg": sigma68(results["delta_deg_route_a"].to_numpy()),
        "sigma68_route_a_vertical_deg": sigma68(
            results.loc[results["theta_deg"] < 1.0, "delta_deg_route_a"].to_numpy()
        ) if (results["theta_deg"] < 1.0).any() else float("nan"),
        "sigma68_by_energy": {str(k): float(v) for k, v in e_stats.items()},
        "figure": str(ang_fig),
        **route_b_summary,
    }
    save_summary("angular_resolution", summary_payload)
    print(json.dumps(summary_payload, indent=2))


if __name__ == "__main__":
    main()
