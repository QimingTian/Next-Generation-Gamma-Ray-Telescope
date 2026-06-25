#!/usr/bin/env bash
# Quick ACD validation: gamma + proton with G4ACD=on, photons off.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${BUILD_DIR:-$ROOT/build-hadronic}"

source "$ROOT/setup_env.sh"
export G4MT=0
export G4ACD=on
export G4WRITE_PHOTONS=0
export G4FILTER=on
export G4OUTPUT_TAG=acd_smoke
export G4SEED="${G4SEED:-424242}"

mkdir -p "$BUILD/data" "$ROOT/data"
cp -f "$ROOT/macros/acd_validation.mac" "$BUILD/macros/"

(
  cd "$BUILD"
  ./MainHad macros/acd_validation.mac
)
cp -f "$BUILD/data/acd_smoke_"* "$ROOT/data/" 2>/dev/null || true

(
  cd "$ROOT"
  python3 - <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, "analysis")
from common import parse_geant4_csv

root = Path("data")
files = sorted(root.glob("acd_smoke_run*_nt_events.csv"))
if not files:
    print("No ACD smoke CSV found")
    raise SystemExit(0)

rows = []
for f in files:
    df = parse_geant4_csv(f)
    run = int(f.name.split("_run")[1].split("_")[0])
    label = {0: "gamma", 1: "proton100", 2: "proton300"}.get(run, f"run{run}")
    df["label"] = label
    rows.append(df)
all_df = __import__("pandas").concat(rows, ignore_index=True)
print("=== ACD smoke summary ===")
for label, g in all_df.groupby("label"):
    veto = g["ACD_veto"].astype(bool)
    print(
        f"{label:10s} n={len(g):3d}  veto_frac={veto.mean():.2f}  "
        f"median_tiles={g.loc[veto, 'ACD_n_tiles'].median() if veto.any() else 0:.0f}  "
        f"median_edep_MeV={g.loc[veto, 'ACD_edep_MeV'].median() if veto.any() else 0:.1f}"
    )
PY
)

echo "Wrote CSVs under $ROOT/data/acd_smoke_run*"
