#!/usr/bin/env python3
"""Generate energy calibration and resolution macros."""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MACROS = ROOT / "macros"


def write_energy_scan(events: int, path: Path) -> None:
    energies = [50, 100, 200, 300, 500, 1000]
    lines = [
        "/run/initialize", "",
        "/gps/particle gamma", "/gps/pos/type Point", "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d", "/gps/direction 0 0 -1", "",
    ]
    for e in energies:
        lines += [f"/gps/energy {e} GeV", f"/run/beamOn {events}", ""]
    path.write_text("\n".join(lines) + "\n")


def write_energy_resolution(energy: float, events: int, path: Path) -> None:
    lines = [
        "/run/initialize", "",
        "/gps/particle gamma", f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Point", "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d", "/gps/direction 0 0 -1",
        f"/run/beamOn {events}", "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cal-events", type=int, default=15)
    parser.add_argument("--res-events", type=int, default=60)
    args = parser.parse_args()

    write_energy_scan(args.cal_events, MACROS / "energy_scan.mac")
    write_energy_resolution(100, args.res_events, MACROS / "energy_resolution_100.mac")
    write_energy_resolution(215, args.res_events, MACROS / "energy_resolution_215.mac")
    write_energy_resolution(1000, 15, MACROS / "peak_factor.mac")
    print("Wrote energy_scan.mac, energy_resolution_*.mac, peak_factor.mac")


if __name__ == "__main__":
    main()
