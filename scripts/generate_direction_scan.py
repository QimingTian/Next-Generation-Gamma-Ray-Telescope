#!/usr/bin/env python3
"""Generate direction-scan Geant4 macros (single file or cloud-sharded grid)."""
from __future__ import annotations

import argparse
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MACROS = ROOT / "macros"
PHASE2 = MACROS / "phase2"


def header() -> list[str]:
    return [
        "/run/initialize",
        "",
        "/gps/particle gamma",
        "/gps/pos/type Point",
        "/gps/pos/centre 0 0 60 cm",
        "",
    ]


def gps_beam(energy: float, theta: float, phi: float) -> list[str]:
    lines = [
        f"# E={energy:g} GeV theta={theta:g} phi={phi:g}",
        f"/gps/energy {energy:g} GeV",
        "/gps/ang/type beam1d",
    ]
    if theta == 0:
        lines.append("/gps/direction 0 0 -1")
    else:
        if phi != 0:
            lines.append(f"/gps/ang/rot2 0 0 {phi:g} deg")
        lines.append(f"/gps/ang/rot1 0 {-theta:g} 0 deg")
        lines.append("/gps/direction 0 0 -1")
    return lines


def default_phis(n: int = 24) -> list[float]:
    """Uniform azimuth samples in [0, 360)."""
    return [i * 360.0 / n for i in range(n)]


def write_cloud_shards(
    energies: list[float],
    thetas: list[float],
    phis: list[float],
    events_per_bin: int,
    subshards: int,
) -> None:
    """432 bins × subshards → one serial Geant4 job each (128-way pool on cloud)."""
    PHASE2.mkdir(parents=True, exist_ok=True)
    manifest: list[str] = []
    total_events = 0

    for energy in energies:
        for theta in thetas:
            for phi in phis:
                base, rem = divmod(events_per_bin, subshards)
                for i in range(subshards):
                    n_ev = base + (1 if i < rem else 0)
                    if n_ev <= 0:
                        continue
                    name = f"dir_E{energy:g}_t{theta:g}_p{phi:g}_s{i:03d}.mac"
                    tag = f"direction_scan_E{energy:g}_t{theta:g}_p{phi:g}_s{i:03d}"
                    body = header() + gps_beam(energy, theta, phi)
                    body.append(f"/run/beamOn {n_ev}")
                    body.append("")
                    (PHASE2 / name).write_text("\n".join(body))
                    manifest.append(f"{tag}|phase2/{name}|on|em")
                    total_events += n_ev

    (PHASE2 / "cloud_manifest.txt").write_text("\n".join(manifest) + "\n")
    n_bins = len(energies) * len(thetas) * len(phis)
    print(f"Wrote {len(manifest)} jobs to {PHASE2} ({n_bins} bins × {subshards} subshards)")
    print(f"Grid: {len(energies)} E × {len(thetas)} θ × {len(phis)} φ = {n_bins} bins")
    print(f"Total primaries: {total_events}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--energies", nargs="+", type=float, default=[30, 100, 300])
    parser.add_argument("--thetas", nargs="+", type=float, default=[0, 5, 10, 20, 30, 45])
    parser.add_argument("--phis", nargs="+", type=float, default=None)
    parser.add_argument("--n-phi", type=int, default=24, help="Azimuth bins when --phis not set")
    parser.add_argument("--events", type=int, default=12, help="Events per bin (single macro mode)")
    parser.add_argument(
        "--events-per-bin",
        type=int,
        default=500,
        help="Events per (E,θ,φ) bin in --cloud-shard mode",
    )
    parser.add_argument(
        "--subshards",
        type=int,
        default=4,
        help="Split each bin into N parallel sub-macros (cloud mode)",
    )
    parser.add_argument(
        "--cloud-shard",
        action="store_true",
        help="Write macros/phase2/*.mac + cloud_manifest.txt for MaxCloudON",
    )
    parser.add_argument("--output", type=Path, default=MACROS / "direction_scan.mac")
    args = parser.parse_args()

    phis = args.phis if args.phis is not None else default_phis(args.n_phi)

    if args.cloud_shard:
        write_cloud_shards(args.energies, args.thetas, phis, args.events_per_bin, args.subshards)
        return

    lines = header()
    for energy in args.energies:
        for theta in args.thetas:
            for phi in phis:
                lines.extend(gps_beam(energy, theta, phi))
                lines.append(f"/run/beamOn {args.events}")
                lines.append("")

    args.output.write_text("\n".join(lines) + "\n")
    n_bins = len(args.energies) * len(args.thetas) * len(phis)
    print(f"Wrote {args.output} ({n_bins * args.events} events)")


if __name__ == "__main__":
    main()
