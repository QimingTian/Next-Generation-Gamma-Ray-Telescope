# ACD design notes (Fermi-LAT baseline)

## Fermi-LAT ACD (reference)

- Segmented **plastic scintillator** (polystyrene-family) panels on **top + four sides**
- **SiPM/PMT** readout per tile; charged-particle MIP → veto
- Segmentation reduces **self-veto** from electromagnetic shower secondaries at high energy
- Thin panels (~1 cm) minimize **gamma conversion** in veto material

## This telescope: copy vs adapt

| Aspect | Fermi-LAT | LXe cube (this sim) | Action |
|--------|-----------|---------------------|--------|
| Panel material | Plastic scintillator | Same (`G4_PLASTIC_SC_VINYLTOLUENE`) | Copy |
| Segmentation | ~50 cm tiles, many channels | **16×16 / face** (~6 cm tiles on 1 m face) | Scale to 1 m |
| Face coverage | Top + 4 sides | **All 6 faces** (cubic LXe) | **Add bottom face** |
| Readout | PMT/SiPM optical | Geant4 MVP: **edep threshold / tile** | Upgrade to optical later |
| Energy range | 100 MeV – 300 GeV | 50 GeV – TeV | Keep segmentation; monitor self-veto @ TeV |
| Companion tracker | Silicon tracker tracks | **None** — LXe Cherenkov only | Veto = ACD + LXe pattern |
| Panel thickness | ~1 cm | **1 cm** | Copy |

## Minimal improvements over Fermi (recommended)

1. **Six-face wrap** — required for cubic geometry (not optional).
2. **Finer tile grid** (16×16 vs Fermi’s coarser top modules) — same self-veto mitigation at TeV.
3. **No change** to material or thickness for first simulation pass.
4. **Future (not MVP):** optical SiPM coupling, timing coincidence, side/top threshold tuning for TeV self-veto studies.

## Simulation flags

- `G4ACD=on` (default): build segmented ACD outside Al shell
- `G4ACD=0|off`: legacy geometry without ACD (reproducibility)

## Veto logic (analysis)

- Per tile: fire if charged-energy deposit ≥ **0.15 MeV** (MIP-scale)
- Event veto if **≥1 tile** fires (`ACD_veto=1`)
- γ/p separation: reject `ACD_veto==1` for gamma sample efficiency vs proton rejection
