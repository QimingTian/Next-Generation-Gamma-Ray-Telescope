#!/usr/bin/env python3
"""Per-SiPM photon density maps for six detector faces."""
from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from common import (
    FIGURES_DIR,
    ensure_dirs,
    face_grid_counts,
    find_run_files,
    parse_geant4_csv,
    peak_factor,
    save_summary,
)

FACE_NAMES = ["+Z", "-Z", "+X", "-X", "+Y", "-Y"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=None, help="Run output tag (e.g. angle_scan)")
    parser.add_argument("--run-id", type=int, default=0)
    parser.add_argument("--six-name", default="density_maps_six_faces.png",
                        help="Output filename for the six-face panel")
    parser.add_argument("--marked-name", default="DensityMapVerticalMarked.png",
                        help="Output filename for the combined +X/-X marked map")
    parser.add_argument("--label", default="Vertical", help="Title label for the marked map")
    args = parser.parse_args()

    ensure_dirs()
    files = find_run_files(args.tag, args.run_id)
    summary = parse_geant4_csv(files["sipm_summary"])

    faces = face_grid_counts(summary)
    pf = peak_factor(summary)

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for ax, (fid, grid) in zip(axes.ravel(), faces.items()):
        im = ax.imshow(grid.values, origin="lower", cmap="viridis")
        ax.set_title(f"Face {FACE_NAMES[fid]}")
        ax.set_xlabel("SiPM column")
        ax.set_ylabel("SiPM row")
        fig.colorbar(im, ax=ax, fraction=0.046, label="photons")

    fig.suptitle(f"SiPM photon density (peak factor = {pf:.2f}×)")
    fig.tight_layout()
    out = FIGURES_DIR / args.six_name
    fig.savefig(out, dpi=150)
    plt.close(fig)

    # Marked map: sum over ±X faces for paper-style figure
    vx = faces[2].values + faces[3].values
    fig2, ax2 = plt.subplots(figsize=(8, 7))
    im2 = ax2.imshow(vx, origin="lower", cmap="hot")
    ax2.set_title(f"{args.label} density map (+X and -X faces combined)")
    ax2.set_xlabel("SiPM column")
    ax2.set_ylabel("SiPM row")
    fig2.colorbar(im2, ax=ax2, label="photons")
    fig2.tight_layout()
    out2 = FIGURES_DIR / args.marked_name
    fig2.savefig(out2, dpi=150)
    plt.close(fig2)

    save_summary("density_maps", {"peak_factor": pf, "figure": str(out), "marked_figure": str(out2)})
    print(f"Wrote {out} and {out2}; peak factor = {pf:.2f}×")


if __name__ == "__main__":
    main()
