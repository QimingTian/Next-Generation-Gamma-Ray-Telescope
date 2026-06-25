#!/usr/bin/env python3
"""Generate macros for paper-gap supplements (A_eff E/theta, proton stats)."""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "macros" / "gaps"


def header() -> list[str]:
    return ["/run/initialize", ""]


def gamma_plane_vertical(energy: float, events: int) -> list[str]:
    return [
        "/gps/particle gamma",
        f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Plane",
        "/gps/pos/shape Square",
        "/gps/pos/halfx 80 cm",
        "/gps/pos/halfy 80 cm",
        "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d",
        "/gps/direction 0 0 -1",
        f"/run/beamOn {events}",
        "",
    ]


def gamma_plane_tilted(energy: float, theta_deg: float, events: int) -> list[str]:
    lines = [
        "/gps/particle gamma",
        f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Plane",
        "/gps/pos/shape Square",
        "/gps/pos/halfx 80 cm",
        "/gps/pos/halfy 80 cm",
        "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d",
        f"/gps/ang/rot1 0 {-theta_deg:g} 0 deg",
        "/gps/direction 0 0 -1",
        f"/run/beamOn {events}",
        "",
    ]
    return lines


def proton_vertical(energy: float, events: int) -> list[str]:
    return [
        "/gps/particle proton",
        f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Point",
        "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d",
        "/gps/direction 0 0 -1",
        f"/run/beamOn {events}",
        "",
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aeff-energy-events", type=int, default=1200)
    parser.add_argument("--aeff-theta-events", type=int, default=600)
    parser.add_argument("--proton-events", type=int, default=50)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    manifest: list[str] = []

    for e in (100, 300, 500, 1000):
        name = f"aeff_E{e:g}GeV.mac"
        lines = header() + gamma_plane_vertical(e, args.aeff_energy_events)
        (OUT / name).write_text("\n".join(lines))
        manifest.append(f"aeff_energy_{e:g}GeV|{name}|on|em")

    for theta in (15, 30, 45):
        name = f"aeff_theta{theta:g}deg.mac"
        lines = header() + gamma_plane_tilted(200, theta, args.aeff_theta_events)
        (OUT / name).write_text("\n".join(lines))
        manifest.append(f"aeff_theta_{theta:g}deg|{name}|on|em")

    for e in (100, 300, 1000):
        name = f"proton_supp_{e:g}GeV.mac"
        lines = header() + proton_vertical(e, args.proton_events)
        (OUT / name).write_text("\n".join(lines))
        manifest.append(f"proton_supplement|{name}|on|hadronic")

    (OUT / "manifest.txt").write_text("\n".join(manifest) + "\n")
    n_em = 4 * args.aeff_energy_events + 3 * args.aeff_theta_events
    n_had = 3 * args.proton_events
    print(f"Wrote {len(manifest)} jobs to {OUT}")
    print(f"EM primaries: {n_em}, hadronic: {n_had}, total: {n_em + n_had}")


if __name__ == "__main__":
    main()
