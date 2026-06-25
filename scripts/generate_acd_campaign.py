#!/usr/bin/env python3
"""Generate ACD-on gamma/proton validation macros."""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "macros" / "acd"


def header() -> list[str]:
    return ["/run/initialize", ""]


def gamma_beam(energy: float, theta_deg: float, events: int) -> list[str]:
    lines = [
        "/gps/particle gamma",
        f"/gps/energy {energy:g} GeV",
        "/gps/pos/type Point",
        "/gps/pos/centre 0 0 60 cm",
        "/gps/ang/type beam1d",
    ]
    if theta_deg != 0:
        lines.append(f"/gps/ang/rot1 0 {-theta_deg:g} 0 deg")
    lines += ["/gps/direction 0 0 -1", f"/run/beamOn {events}", ""]
    return lines


def proton_beam(energy: float, events: int) -> list[str]:
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
    parser.add_argument("--gamma-events", type=int, default=50)
    parser.add_argument("--gamma-tilt-events", type=int, default=30)
    parser.add_argument("--proton-events", type=int, default=50)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    manifest: list[str] = []

    for e in (100, 300, 500):
        name = f"gamma_E{e:g}GeV.mac"
        (OUT / name).write_text("\n".join(header() + gamma_beam(e, 0, args.gamma_events)))
        manifest.append(f"acd_gamma_{e:g}GeV|{name}|on|em")

    for theta in (15, 30, 45):
        name = f"gamma_100GeV_theta{theta:g}deg.mac"
        (OUT / name).write_text("\n".join(header() + gamma_beam(100, theta, args.gamma_tilt_events)))
        manifest.append(f"acd_gamma_100_theta{theta:g}|{name}|on|em")

    for e in (100, 300, 1000):
        name = f"proton_E{e:g}GeV.mac"
        (OUT / name).write_text("\n".join(header() + proton_beam(e, args.proton_events)))
        manifest.append(f"acd_proton_{e:g}GeV|{name}|on|hadronic")

    (OUT / "manifest.txt").write_text("\n".join(manifest) + "\n")
    print(f"Wrote {len(manifest)} macros to {OUT}")


if __name__ == "__main__":
    main()
