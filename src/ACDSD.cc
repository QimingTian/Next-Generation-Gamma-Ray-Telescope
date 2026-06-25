#include "ACDSD.hh"

#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"

ACDSD::ACDSD(const G4String& name) : G4VSensitiveDetector(name) {}

void ACDSD::Initialize(G4HCofThisEvent* /*hce*/)
{
  fTileEdep.clear();
  fTotalEdep = 0.;
}

G4bool ACDSD::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/)
{
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) return false;

  const auto* track = step->GetTrack();
  if (track->GetDefinition()->GetPDGCharge() == 0.) {
    return false;
  }

  const G4int tileID = step->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
  fTileEdep[tileID] += edep;
  fTotalEdep += edep;
  return true;
}

G4int ACDSD::GetNTilesHit() const
{
  G4int n = 0;
  for (const auto& entry : fTileEdep) {
    if (entry.second / MeV >= kTileThresholdMeV) {
      ++n;
    }
  }
  return n;
}
