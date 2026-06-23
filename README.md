# Next-Generation Gamma-Ray Telescope

[![Geant4](https://img.shields.io/badge/Geant4-11.3-blue.svg)](https://geant4.web.cern.ch/)
[![Python](https://img.shields.io/badge/Python-3.8+-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Geant4-orange.svg)](http://cern.ch/geant4/license)

A publication-grade Monte Carlo framework for a **space-based liquid xenon (LXe) Cherenkov telescope** with interior SiPM readout. The project couples a rebuilt **Geant4 11.3** detector simulation, **transfer-matrix optical filter design**, and a **Python analysis pipeline** that produces all figures in the companion paper.

**Repository:** [github.com/QimingTian/Next-Generation-Gamma-Ray-Telescope](https://github.com/QimingTian/Next-Generation-Gamma-Ray-Telescope)

---

## Table of Contents

- [Scientific Overview](#scientific-overview)
- [What This Repository Contains](#what-this-repository-contains)
- [Key Performance Results](#key-performance-results)
- [Project Structure](#project-structure)
- [Detector and Physics Model](#detector-and-physics-model)
- [Optical Filter Design](#optical-filter-design)
- [Installation](#installation)
- [Building the Simulation](#building-the-simulation)
- [Running Simulations](#running-simulations)
- [Environment Variables](#environment-variables)
- [Simulation Output](#simulation-output)
- [Analysis Pipeline](#analysis-pipeline)
- [Campaign Scripts](#campaign-scripts)
- [Long-Running Jobs](#long-running-jobs)
- [Paper and Figures](#paper-and-figures)
- [Reproducibility](#reproducibility)
- [Requirements](#requirements)
- [License](#license)

---

## Scientific Overview

This instrument concept targets **indirect dark-matter searches** and high-energy astrophysics in the **50 GeV–1 TeV** band, where existing space missions (e.g. *Fermi*-LAT) lose sensitivity and ground-based IACTs have limited duty cycle and field of view.

The simulated detector is a **1 m³ cubic LXe volume** instrumented with **1,944 SiPMs** (18×18 per face) mounted on the interior surfaces. Electromagnetic showers from gamma rays produce:

- **Cherenkov light** (190–600 nm) — directional, used for tracking and energy reconstruction
- **Scintillation light** (~175 nm) — isotropic background from LXe excitation

A **CaF₂ + MgF₂ optical filter** on each SiPM face transmits Cherenkov photons while suppressing scintillation. Per-photon filter transmission is applied in the sensitive detector (`SiPMSD`), so filter ON/OFF comparisons are first-class simulation modes.

The framework supports:

| Capability | Description |
|------------|-------------|
| Energy calibration | Linear photon yield vs. primary energy |
| Energy resolution | Poisson-limited σ_E/E at 100 and 215 GeV |
| Shower profiling | Longitudinal dE/dz, L₉₀ containment depth |
| Angular reconstruction | Two reconstruction routes from SiPM count maps |
| SiPM saturation | Peak-factor–limited dynamic range |
| Effective area | Plane-source Monte Carlo |
| Background rejection | γ vs. proton face-anisotropy comparison |

Full scientific motivation, detector rationale, and figure captions are in [`paper/Paper.tex`](paper/Paper.tex).

---

## What This Repository Contains

The codebase was rebuilt in 2026 for reproducible, publication-grade results. Compared with the original prototype, the current version provides:

1. **Correct Geant4 architecture** — sensitive detectors registered in `ConstructSDandField()`, serial batch mode with merged CSV output, GPS-driven macros (no hardcoded beam).
2. **Structured ntuples** — per-photon, per-event, and per-SiPM summaries written via `G4AnalysisManager`.
3. **Runtime configuration** — filter, photon output, output tags, and random seed controlled by environment variables.
4. **Campaign automation** — shell scripts for overnight energy, direction, effective-area, and hadronic background runs.
5. **End-to-end analysis** — Python scripts from raw CSV to figures in `figures/` and `paper/`.
6. **Hadronic physics build** — optional `MainHad` binary with `FTFP_BERT` for proton background studies.

See [`docs/data_provenance.md`](docs/data_provenance.md) for a detailed comparison of pre-rebuild vs. current simulation values.

---

## Key Performance Results

Representative results from the June 2026 campaign (`G4SEED=123456789`, filter design in `analysis/output/filter_transmission.csv`):

| Quantity | Value | Notes |
|----------|-------|-------|
| Photon yield (filter OFF) | 9,781 photons/GeV | Energy calibration slope |
| Photon yield (filter ON) | 6,982 photons/GeV | Per-photon T(λ) in `SiPMSD` |
| Mean filter transmission | ~0.71 | Simulation-integrated |
| σ_E/E @ 215 GeV (filter ON) | ~0.50% | 40-event resolution macro |
| σ_E/E @ 100 GeV (filter ON) | ~0.65% | 40-event resolution macro |
| Peak factor | ~2.2× | Peak-to-mean SiPM occupancy |
| Saturation energy E_sat | ~1.63 TeV | 18×18×6, 14,400-cell SiPMs |
| L₉₀ shower depth | ~595 mm | Longitudinal profile |
| Angular σ₆₈ (vertical, Route A) | ~0.13° | Direction scan campaign |
| Effective area @ 200 GeV | ~0.88 m² | Plane-source MC |
| γ/p face-anisotropy ratio | ~2.6× | Background rejection metric |

---

## Project Structure

```
GammaRayTelescope/
├── Main.cc                     # Geant4 application entry point
├── CMakeLists.txt              # Builds Main (+ MainHad with -DWITH_HADRONIC=ON)
├── setup_env.sh                # Sources local Geant4 install and toolchain paths
│
├── include/                    # Detector and runtime headers
│   ├── DetectorConstruction.hh
│   ├── PhysicsList.hh
│   ├── SiPMSD.hh              # Optical photon sensitive detector + filter
│   ├── AnalysisManager.hh     # CSV ntuple booking
│   ├── FilterTransmission.hh  # Wavelength-dependent T(λ) lookup
│   ├── OpticalMaterials.hh    # LXe Sellmeier / absorption / Rayleigh
│   └── RuntimeConfig.hh       # G4WRITE_PHOTONS, G4FILTER flags
│
├── src/                        # Geant4 implementation (.cc)
│
├── macros/                     # GPS batch macros (copied into build dirs)
│   ├── energy_scan.mac         # 50–1000 GeV calibration scan
│   ├── energy_resolution_*.mac # Fixed-energy resolution samples
│   ├── direction_scan.mac      # Angular reconstruction training set
│   ├── peak_factor.mac         # High-statistics shower for density maps
│   ├── plane_source.mac        # Effective-area plane source
│   ├── proton_background.mac   # Hadronic background (MainHad)
│   └── quick_test.mac          # Single-event smoke test
│
├── analysis/                   # Python post-processing
│   ├── common.py               # CSV parsing, SiPM grid utilities
│   ├── Filter-Design.py        # TMM filter spectrum → filter_transmission.csv
│   ├── run_all.py              # Run full analysis chain
│   ├── recon/                  # Angular reconstruction (Route A & B)
│   └── *.py                    # Individual figure/metric scripts
│
├── scripts/                    # Campaign orchestration
│   ├── run_campaign.sh         # Core energy + angle + peak campaigns
│   ├── run_overnight.sh        # Full paper campaign (direction, A_eff, background)
│   ├── run_energy_supplement.sh
│   ├── finish_campaigns.sh
│   ├── sync_campaign_data.sh
│   ├── generate_direction_scan.py
│   ├── generate_energy_campaign.py
│   └── keep_awake.sh           # Prevent sleep during long runs
│
├── data/                       # Campaign CSV output (gitignored except meta JSON)
├── figures/                    # Generated plots (gitignored)
├── analysis/output/            # JSON summaries from analysis scripts
├── paper/                      # LaTeX manuscript + figure PNGs
├── docs/                       # Data provenance and methodology notes
│
├── build/                      # Standard EM build (gitignored)
├── build-hadronic/             # Hadronic build with MainHad (gitignored)
└── deps/                       # Local Geant4 source/build/install (gitignored)
```

---

## Detector and Physics Model

### Geometry

| Component | Specification |
|-----------|---------------|
| World | 1.2 m cube, vacuum (`G4_Galactic`) |
| Aluminum shell | 1.02 m cube, structural support |
| LXe active volume | 1.0 m cube, ρ = 2.953 g/cm³ at 165 K |
| SiPM array | 18×18 per face × 6 faces = **1,944 channels** |
| SiPM pixel | 1.0 × 1.0 cm², 1 mm thick |
| SiPM pitch | 5 cm center-to-center (~14.6% fill factor per face) |

SiPM placement, optical surfaces, and material properties are implemented in `src/DetectorConstruction.cc` and `src/OpticalMaterials.cc`.

### Physics

- **Electromagnetic:** `G4EmStandardPhysics_option4`
- **Optical:** official `G4OpticalPhysics` — Cherenkov, scintillation, Rayleigh scattering, boundary processes
- **Hadronic (optional):** `FTFP_BERT` when built with `-DWITH_HADRONIC=ON` (`MainHad`)
- **Primaries:** `G4GeneralParticleSource` via macro files

### Executables

| Binary | Build flag | Purpose |
|--------|------------|---------|
| `Main` | default (`WITH_HADRONIC=OFF`) | Gamma-ray and EM shower campaigns |
| `MainHad` | `-DWITH_HADRONIC=ON` | Proton background with hadronic physics |

---

## Optical Filter Design

The filter stack (air → LXe) is:

```
Air (n = 1.0)
  → MgF₂ AR layer 1 (~20 nm)
  → MgF₂ AR layer 2 (~20 nm)
  → CaF₂ substrate (~0.2 mm)
  → Liquid xenon
```

**Design tool:** transfer matrix method (TMM) in `analysis/Filter-Design.py`, with optional AR coating optimization in `analysis/AR-Coating-Optimization.py`.

**Simulation integration:** `FilterTransmission` loads `analysis/output/filter_transmission.csv`. In `SiPMSD::ProcessHits`, each optical photon is accepted with probability T(λ); scintillation at 175 nm is strongly suppressed while Cherenkov light at 190–600 nm passes at >95%.

Regenerate the filter table before campaigns:

```bash
source .venv/bin/activate
python analysis/Filter-Design.py
```

Outputs: `analysis/output/filter_transmission.csv`, `filter_design.png`, `filter_design.pdf`.

---

## Installation

### 1. Geant4

This project expects **Geant4 11.3.x** with datasets installed. A local install path is wired in `setup_env.sh`:

```bash
source setup_env.sh   # sets Geant4_DIR, CMAKE_PREFIX_PATH, G4SEED
```

If you build Geant4 from source, install to `deps/geant4-install/` or edit `setup_env.sh` to point at your installation.

**macOS (Homebrew alternative):**

```bash
brew install geant4 cmake qt@5 expat
```

### 2. Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Dependencies: NumPy, SciPy, Matplotlib, Pandas, scikit-learn.

---

## Building the Simulation

### Standard build (gamma / EM campaigns)

```bash
source setup_env.sh
mkdir -p build && cd build
cmake .. -DWITH_GEANT4_UIVIS=OFF
cmake --build . -j$(sysctl -n hw.ncpu)    # macOS
# cmake --build . -j$(nproc)              # Linux
```

### Hadronic build (proton background)

```bash
mkdir -p build-hadronic && cd build-hadronic
cmake .. -DWITH_GEANT4_UIVIS=OFF -DWITH_HADRONIC=ON
cmake --build . -j$(sysctl -n hw.ncpu)
```

Macros are copied into each build directory at configure time. After editing macros in `macros/`, re-run `cmake ..` or copy manually:

```bash
cp macros/*.mac build/macros/
cp macros/*.mac build-hadronic/macros/
```

---

## Running Simulations

### Quick smoke test

```bash
source setup_env.sh
cd build
export G4OUTPUT_TAG=quick_test G4FILTER=on G4WRITE_PHOTONS=0
./Main macros/quick_test.mac
ls data/quick_test_run0_nt_*.csv
```

### Single campaign macro

```bash
source setup_env.sh
cd build
export G4SEED=123456789
export G4OUTPUT_TAG=energy_scan_on
export G4MACRO=energy_scan.mac
export G4FILTER=on
export G4WRITE_PHOTONS=0
./Main macros/energy_scan.mac
cp data/energy_scan_on_* ../data/
```

### Interactive mode (with visualization)

Build with `-DWITH_GEANT4_UIVIS=ON`, then:

```bash
cd build
./Main
# In the Geant4 session:
/control/execute macros/quick_test.mac
```

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `G4SEED` | `123456789` | CLHEP random seed (set in `Main.cc`) |
| `G4OUTPUT_TAG` | `default` | Prefix for output files (`{tag}_run{N}_nt_*.csv`) |
| `G4MACRO` | — | Recorded in run metadata JSON |
| `G4FILTER` | `on` | `on`/`off` — per-photon filter rejection in `SiPMSD` |
| `G4FILTER_CSV` | `../analysis/output/filter_transmission.csv` | Filter transmission table |
| `G4WRITE_PHOTONS` | `0` | `1` to write per-photon ntuple (slower, large files) |

Filter OFF is used for scintillation-inclusive baseline runs; filter ON is the nominal instrument mode.

---

## Simulation Output

Each run produces CSV ntuples under `build/data/` (copied to `data/` by campaign scripts):

### `{tag}_run{N}_nt_events.csv`

Per-event summary:

| Column | Description |
|--------|-------------|
| `EventID` | Event index |
| `E_primary_GeV` | Primary gamma energy |
| `DirX`, `DirY`, `DirZ` | Primary direction |
| `N_Cherenkov`, `N_Scint` | Photon counts by process |
| `Edep_MeV` | Total energy deposit |
| `ShowerLength_mm` | Longitudinal shower extent |
| `L90_mm`, `L90_truncated` | 90% energy containment depth |

### `{tag}_run{N}_nt_sipm_summary.csv`

| Column | Description |
|--------|-------------|
| `EventID` | Event index |
| `SiPMID` | SiPM copy number (0–1943) |
| `PhotonCount` | Detected photons on that SiPM |

### `{tag}_run{N}_nt_photons.csv` (optional, `G4WRITE_PHOTONS=1`)

Per-photon hits: position, time, wavelength, creating process (Cherenkov / Scintillation).

### `{tag}_run{N}_meta.json`

Run provenance: macro name, filter flag, seed, timestamp, git commit (when available).

---

## Analysis Pipeline

After simulation data is in `data/`, run the full analysis chain:

```bash
source .venv/bin/activate
python analysis/run_all.py
```

Or run individual scripts:

```bash
python analysis/energy_calibration.py --data-dir data
python analysis/energy_resolution.py --tag energy_resolution_215_on
python analysis/energy_resolution_filter_compare.py --data-dir data
python analysis/shower_length.py --data-dir data --tag energy_scan_on
python analysis/saturation.py
python analysis/density_maps.py --tag peak_factor
python analysis/angular_resolution.py --tag direction_scan --train-route-b
python analysis/effective_area.py --tag plane_source --data-dir data
python analysis/background_rejection.py --data-dir data
python analysis/filter_folding.py --tag spectrum_photons
```

**Outputs:**

- Figures → `figures/*.png` (copied to `paper/` for the manuscript)
- JSON summaries → `analysis/output/*.json`

### Angular reconstruction

Two routes in `analysis/recon/`:

- **Route A** — geometric luminous-axis reconstruction from face-weighted centroids (accurate for near-vertical incidence)
- **Route B** — machine-learning regression on SiPM count features (trained on `direction_scan` campaign)

---

## Campaign Scripts

| Script | Purpose |
|--------|---------|
| `scripts/run_campaign.sh` | Core campaigns: energy scan, resolution, angle scan, peak factor + basic analysis |
| `scripts/run_overnight.sh` | Full paper campaign: direction scan, filter on/off resolution, calibration, peak factor, plane source, proton background, full analysis |
| `scripts/run_energy_supplement.sh` | Re-run energy metrics after filter/SiPMSD fixes |
| `scripts/finish_campaigns.sh` | Plane source + sync + final analysis pass |
| `scripts/sync_campaign_data.sh` | Copy latest CSV from `build/` and `build-hadronic/data/` to `data/` |
| `scripts/generate_direction_scan.py` | Generate inclined-beam direction-scan macro |
| `scripts/generate_energy_campaign.py` | Generate energy-resolution macros with configurable event counts |

**Recommended full reproduction:**

```bash
bash scripts/run_overnight.sh
# or, after partial runs:
bash scripts/run_energy_supplement.sh
bash scripts/finish_campaigns.sh all
```

`run_overnight.sh` automatically wraps itself with `caffeinate` on macOS to reduce idle sleep during long runs.

---

## Long-Running Jobs

For multi-hour Geant4 campaigns with the laptop lid closed:

```bash
# 1. Plug into AC power
# 2. Start keep-awake (attaches caffeinate to running Main/MainHad PIDs)
bash scripts/keep_awake.sh start

# 3. Prevent lid-close sleep (run once in Terminal; requires password)
sudo pmset -a disablesleep 1

# Check status
bash scripts/keep_awake.sh status

# When finished, restore normal power management
bash scripts/keep_awake.sh stop
sudo pmset -a disablesleep 0
```

---

## Paper and Figures

The LaTeX manuscript lives in [`paper/Paper.tex`](paper/Paper.tex):

> *A Liquid Xenon Cherenkov Gamma-Ray Telescope: Simulation Performance of a Cubic SiPM Instrument*

Analysis scripts write PNG figures to `figures/`; campaign scripts copy them into `paper/` for inclusion in the manuscript. Regenerate figures after updating simulation data:

```bash
python analysis/run_all.py
cp figures/*.png paper/
```

---

## Reproducibility

1. Set a fixed seed: `export G4SEED=123456789`
2. Record Geant4 version (11.3.2) and git commit — stored in `*_meta.json`
3. Regenerate filter table: `python analysis/Filter-Design.py`
4. Run campaigns with the same macros and `G4FILTER` settings documented in `docs/data_provenance.md`

To compare against published numbers, see the old-vs-new table in [`docs/data_provenance.md`](docs/data_provenance.md).

---

## Requirements

### C++ / Geant4

| Component | Version |
|-----------|---------|
| Geant4 | 11.3.x (with EM and optical physics datasets) |
| CMake | ≥ 3.16 |
| C++ compiler | C++17 (GCC 9+, Clang 10+, Apple Clang) |
| Qt5 | Required if building with `WITH_GEANT4_UIVIS=ON` |

### Python

| Package | Purpose |
|---------|---------|
| numpy, scipy | Numerics, optimization |
| matplotlib | Figures |
| pandas | CSV ntuple analysis |
| scikit-learn | Route B angular reconstruction |

---

## License

Geant4 application code follows the **[Geant4 Software License](http://cern.ch/geant4/license)**. Acknowledge the Geant4 Collaboration in publications that use this simulation.

Python analysis scripts and project-specific C++ extensions are provided for research use alongside the manuscript. See individual file headers where applicable.

---

## Acknowledgments

- **Geant4 Collaboration** — Monte Carlo toolkit
- **Hamamatsu Photonics** — SiPM specifications used in saturation modeling
- Optical properties compiled from published LXe data (Aprile, Grace, Seidel, and related references cited in `paper/Paper.tex`)

---

**Last updated:** June 2026
