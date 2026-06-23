#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

namespace B1
{

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fEventAction) return;

  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.) {
    fEventAction->AddEdep(edep);
  }

  const auto* pre = step->GetPreStepPoint();
  const auto* volume = pre->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (volume->GetName() == "LXe") {
    const auto pos = pre->GetPosition();
    fEventAction->UpdateShowerBounds(pos);

    if (edep > 0.) {
      // Longitudinal depth along primary direction from primary vertex.
      const auto* runManager = G4RunManager::GetRunManager();
      const auto* event = runManager->GetCurrentEvent();
      G4ThreeVector primaryDir{0, 0, 1};
      G4ThreeVector primaryPos{0, 0, 0};
      if (event && event->GetPrimaryVertex()) {
        primaryPos = event->GetPrimaryVertex()->GetPosition();
        const auto* particle = event->GetPrimaryVertex()->GetPrimary();
        if (particle) {
          primaryDir = particle->GetMomentumDirection();
        }
      }
      const G4double depth = (pos - primaryPos).dot(primaryDir) / mm;
      fEventAction->AddLongitudinalDeposit(depth, edep / MeV);
    }
  }
}

}  // namespace B1
