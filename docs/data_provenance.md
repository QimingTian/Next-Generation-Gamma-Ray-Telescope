# Data Provenance and Old vs New Comparison

This document tracks the regeneration of simulation results after the publication-grade rebuild (June 2026).

## Why data was regenerated

The previous Geant4 code could not reproduce paper figures:
- `SiPMSD` was created in `Construct()` instead of `ConstructSDandField()` (MT failure)
- Output used broken per-thread buffers; CSV contained headers only
- Primary source was hardcoded (1 GeV, 45°); paper scans required GPS macros
- Analysis scripts for density maps, energy calibration, resolution, shower length, and saturation were never committed
- Direction reconstruction was hand-drawn (TikZ); no quantitative σ68

## Old paper values (pre-rebuild, for reference)

| Quantity | Old paper value | Source |
|----------|----------------|--------|
| Total photons/GeV | 9215 | Fig. energy calibration |
| Filtered photons/GeV | 6730 | Fig. energy calibration |
| Filter fraction | 73.1% | Fig. energy calibration |
| Energy resolution @ 215 GeV | 0.80% | Section 3.5 |
| Peak factor | 3.4× | Saturation section |
| E_sat (18×18×6, 14400 cells) | 0.79 TeV | Fig. saturation |
| Cherenkov angle @ 190 nm | 51.9° | Section 3.2 |
| Angular resolution | TBD (hand-drawn) | TikZ recon figures |

## New results (June 2026 paper-quality parallel campaign)

Campaign: `bash scripts/run_paper_quality_parallel.sh` (12-core serial shards, `G4MT=0`, `caffeinate`). Analysis: `python analysis/run_all.py`.

| Quantity | New value | Notes |
|----------|-----------|-------|
| Total photons/GeV (filter OFF) | 9683 | `energy_calibration.json` slope |
| Total photons/GeV (filter ON) | 6920 | `energy_calibration.json` slope |
| Filter fraction (sim) | 0.715 | per-photon T(λ) in `SiPMSD.cc` |
| σ_E/E @ 100 GeV (ON/OFF) | 0.65% / 0.62% | 40 events each |
| σ_E/E @ 215 GeV (ON/OFF) | 0.50% / 0.55% | 40 events each |
| σ68 angular (vertical, Route A) | 0.13° | `direction_scan`, n=36 @ θ=0° |
| σ68 angular (inclined, Route A) | ≈ θ | chord underestimation (5°→5.1°, 45°→45.0°) |
| σ68 angular (Route B, held-out test) | 16.5° | GBR on 20 features, n_test=44 |
| L90 saturation depth | 666 mm | `shower_length.json`, longitudinal profile |
| Peak factor @ 1 TeV | 2.07 ± 0.45 | 30 events, `peak_factor_stats.json` |
| E_sat (18×18×6, 14400 cells) | 1.75 TeV | `saturation.json`, filter ON yield |
| A_eff @ 200 GeV | 1.04 m² | plane source: 2057/5064 on 2.56 m² throw |
| γ/p face-anisotropy ratio | 2.60× | 150 proton + direction_scan gammas |

### Simulation flags (new)

- `G4WRITE_PHOTONS=0` (default): sipm_summary + events only (fast campaigns)
- `G4FILTER=on`: per-photon filter rejection via `SiPMSD` + `FilterTransmission`
- `L90_mm`, `L90_truncated`: longitudinal dE/dz profile in events ntuple

Campaign commands:
- `bash scripts/run_paper_quality_parallel.sh` — full paper-quality suite (final)
- `bash scripts/run_overnight.sh` — direction scan + priority campaigns
- `bash scripts/run_energy_supplement.sh` — filter-fixed energy/L90/peak/A_eff campaigns
- `bash scripts/finish_campaigns.sh` — plane source + final analysis

## Reproducibility

- Set `G4SEED` environment variable before running
- Geant4 version: 11.3.2 (see `deps/geant4-install`)
- Record git commit in each campaign via `*_meta.json`
