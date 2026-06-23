#include "SiPMSD.hh"

#include "FilterTransmission.hh"
#include "RuntimeConfig.hh"

#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <cstring>

SiPMSD::SiPMSD(const G4String& name) : G4VSensitiveDetector(name) {}

void SiPMSD::Initialize(G4HCofThisEvent* /*hce*/)
{
  fHits.clear();
  fSiPMCounts.clear();
  fNCherenkov = 0;
  fNScint = 0;
  fStoreHits = B1::RuntimeConfig::Instance().WritePhotons();
  const char* filterEnv = std::getenv("G4FILTER");
  fApplyFilter =
    !(filterEnv && (std::strcmp(filterEnv, "0") == 0 || std::strcmp(filterEnv, "off") == 0 ||
                    std::strcmp(filterEnv, "OFF") == 0 || std::strcmp(filterEnv, "false") == 0));
  if (fApplyFilter) {
    const char* csv = std::getenv("G4FILTER_CSV");
    B1::FilterTransmission::Instance().Initialize(
      csv ? G4String(csv) : G4String("../analysis/output/filter_transmission.csv"));
  }
}

G4bool SiPMSD::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/)
{
  auto* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    return false;
  }

  const G4double energy_eV = track->GetTotalEnergy() / eV;
  if (energy_eV <= 0.) return false;
  const G4double wavelength_nm = 1240.0 / energy_eV;

  if (fApplyFilter) {
    const G4double T = B1::FilterTransmission::Instance().Transmission(wavelength_nm);
    if (G4UniformRand() > T) {
      return false;
    }
  }

  const auto* creator = track->GetCreatorProcess();
  G4String processName = creator ? creator->GetProcessName() : "Unknown";

  auto* pre = step->GetPreStepPoint();
  const G4int sipmID = pre->GetTouchable()->GetCopyNumber();

  if (processName == "Cerenkov") {
    ++fNCherenkov;
  } else if (processName == "Scintillation") {
    ++fNScint;
  }

  if (fStoreHits) {
    HitData hit;
    hit.sipmID = sipmID;
    hit.position = pre->GetPosition();
    hit.time = pre->GetGlobalTime();
    hit.wavelength_nm = wavelength_nm;
    hit.processName = processName;
    fHits.push_back(hit);
  }
  fSiPMCounts[sipmID]++;

  return true;
}
