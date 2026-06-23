#!/usr/bin/env python3
"""Generate direction-scan Geant4 macro."""
from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MACROS = ROOT / "macros"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--energies", nargs="+", type=float, default=[30, 100, 300])
    parser.add_argument("--thetas", nargs="+", type=float, default=[0, 5, 10, 20, 30, 45])
    parser.add_argument("--phis", nargs="+", type=float, default=[0])
    parser.add_argument("--events", type=int, default=12)
    parser.add_argument("--output", type=Path, default=MACROS / "direction_scan.mac")
    args = parser.parse_args()

    lines = ["/run/initialize", "", "/gps/particle gamma", "/gps/pos/type Point",
             "/gps/pos/centre 0 0 60 cm", ""]

    for energy in args.energies:
        for theta in args.thetas:
            for phi in args.phis:
                lines.append(f"# E={energy:g} GeV theta={theta:g} phi={phi:g}")
                lines.append(f"/gps/energy {energy:g} GeV")
                lines.append("/gps/ang/type beam1d")
                if theta == 0:
                    lines.append("/gps/direction 0 0 -1")
                else:
                    if phi != 0:
                        lines.append(f"/gps/ang/rot2 0 0 {phi:g} deg")
                    # Tilt theta off vertical in the plane defined by rot2 (default X-Z at phi=0).
                    lines.append(f"/gps/ang/rot1 0 {-theta:g} 0 deg")
                    lines.append("/gps/direction 0 0 -1")
                lines.append(f"/run/beamOn {args.events}")
                lines.append("")

    args.output.write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.output} ({len(args.energies)*len(args.thetas)*len(args.phis)*args.events} events)")


if __name__ == "__main__":
    main()
