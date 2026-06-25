#!/usr/bin/env python3
"""Generate per-shard Geant4 macros for parallel execution across CPU cores."""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "macros" / "parallel"


def write_gamma_point(energy: float, events: int, path: Path) -> None:
    lines = [
        "/run/initialize", "",
        "/gps/particle gamma", f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Point", "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d", "/gps/direction 0 0 -1",
        f"/run/beamOn {events}", "",
    ]
    path.write_text("\n".join(lines) + "\n")


def write_gamma_plane(energy: float, events: int, path: Path) -> None:
    lines = [
        "/run/initialize", "",
        "/gps/particle gamma", f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Plane", "/gps/pos/shape Square",
        "/gps/pos/halfx 80 cm", "/gps/pos/halfy 80 cm",
        "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d", "/gps/direction 0 0 -1",
        f"/run/beamOn {events}", "",
    ]
    path.write_text("\n".join(lines) + "\n")


def write_proton(energy: float, events: int, path: Path) -> None:
    lines = [
        "/run/initialize", "",
        "/gps/particle proton", f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Point", "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d", "/gps/direction 0 0 -1",
        f"/run/beamOn {events}", "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cal-events", type=int, default=15)
    parser.add_argument("--peak-events", type=int, default=30)
    parser.add_argument("--plane-events", type=int, default=5000)
    parser.add_argument("--proton-events", type=int, default=50)
    parser.add_argument("--plane-shards", type=int, default=10)
    parser.add_argument("--peak-shards", type=int, default=6)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    manifest: list[str] = []

    for e in [50, 100, 200, 300, 500, 1000]:
        for filt in ("on", "off"):
            name = f"energy_{e:g}GeV_{filt}.mac"
            write_gamma_point(e, args.cal_events, OUT / name)
            manifest.append(f"energy_scan_{filt}|{name}|{filt}")

    per_peak = max(1, args.peak_events // args.peak_shards)
    for i in range(args.peak_shards):
        n = per_peak if i < args.peak_shards - 1 else args.peak_events - per_peak * (args.peak_shards - 1)
        name = f"peak_factor_{i}.mac"
        write_gamma_point(1000, n, OUT / name)
        manifest.append(f"peak_factor|{name}|on")

    per_plane = max(1, args.plane_events // args.plane_shards)
    for i in range(args.plane_shards):
        n = per_plane if i < args.plane_shards - 1 else args.plane_events - per_plane * (args.plane_shards - 1)
        name = f"plane_source_{i}.mac"
        write_gamma_plane(200, n, OUT / name)
        manifest.append(f"plane_source|{name}|on")

    for e in [100, 300, 1000]:
        name = f"proton_{e:g}GeV.mac"
        write_proton(e, args.proton_events, OUT / name)
        manifest.append(f"proton_background|{name}|on")

    # One photons-on run for spectrum figure
    write_gamma_point(1000, 1, OUT / "spectrum_photons.mac")
    manifest.append("spectrum_photons|spectrum_photons.mac|on")

    (OUT / "manifest.txt").write_text("\n".join(manifest) + "\n")
    print(f"Wrote {len(manifest)} parallel jobs to {OUT}")


if __name__ == "__main__":
    main()
