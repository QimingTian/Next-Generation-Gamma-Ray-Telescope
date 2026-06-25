# Paper improvements tracker

Four workstreams for the LXe Cherenkov telescope paper.

---

## 1. Performance gaps — local Mac (running / done)

| ID | Item | Status |
|----|------|--------|
| G6–G9 | Energy calibration, σ_E/E, L90, E_sat, peak factor | **Done** |
| G2 | A_eff(E, θ) — was single point @ 200 GeV vertical | **Partial** — local: 100/300 GeV + θ done; **E=500/1000 GeV → cloud** |
| G3 | γ/p separation (face anisotropy only, n=150) | **Done** (proton supplement 150 events) |
| G10 | Density maps (extra θ examples) | After ML campaign |
| — | `Paper.tex` number sync | Mostly done; refresh after gaps |

**Not rerun unless geometry changes:** L90 body physics, E_geom^max, filter TMM design.

---

## 2. Angular reconstruction ML — cloud CPU (waiting)

| ID | Item | Status |
|----|------|--------|
| G1 | Phase 2 direction scan (432×500 = 216k events, full φ) | **Waiting** TensorDock / AWS / RunPod |
| G5 | Inclined σ68 (Route B test ≈16.5° on 216 events) | Blocked on G1 |
| G16 | Route C CNN (6×18×18, angular loss) | Code + train after G1 |
| G14 | SiPM density 36×36 / 40×40 | **Deferred** until G1 on 18×18 |

---

## 3. ACD anti-coincidence — Geant4 (in progress)

| ID | Item | Status |
|----|------|--------|
| G4 | Fermi-style segmented plastic scintillator + tile SD | **Done** (`G4ACD=on`) |
| G17 | ACD + LXe combined veto ROC | **Done** sim; run `background_rejection.py --combined` |
| G18 | Gamma self-veto vs E, θ | `analysis/acd_veto_scan.py` |

Design basis: Fermi-LAT ACD (segmented panels, six-face wrap). See `docs/acd_design.md`.

---

## 4. Misc / optional

| ID | Item | Status |
|----|------|--------|
| — | Cloud deploy scripts (`run_runpod_bundle.sh`) | Pending machine purchase |
| G12 | Electron background sample | Optional |
| G13 | Systematic sweeps (optics, PDE, filter ±10%) | Optional |
| G15 | Geant4 MT + optical segfault | Workaround: parallel serial |
| — | Engineering feasibility (cryo, structure) | Text only |

---

## Campaign commands

```bash
# Local gaps (A_eff + protons)
bash scripts/run_local_gaps.sh

# ACD smoke test (after rebuild)
bash scripts/run_acd_smoke.sh

# ACD validation campaign (~390 events, auto-limits to 2 cores if gaps running)
bash scripts/run_acd_campaign.sh

# Cloud Phase 2 (when ready)
bash scripts/run_runpod_bundle.sh   # TBD
```
