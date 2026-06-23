#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;

namespace B1
{

class RunAction;

class EventAction : public G4UserEventAction
{
  public:
    explicit EventAction(RunAction* runAction);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void AddEdep(G4double edep) { fEdep += edep; }
    void UpdateShowerBounds(const G4ThreeVector& pos);
    void AddLongitudinalDeposit(G4double depth_mm, G4double edep);
    void SetPrimaryEnergy(G4double e) { fPrimaryEnergy = e; }
    void SetPrimaryDirection(const G4ThreeVector& dir) { fPrimaryDir = dir; }
    void SetPrimaryPosition(const G4ThreeVector& pos) { fPrimaryPos = pos; }

  private:
    G4double ComputeL90() const;

    RunAction* fRunAction = nullptr;
    G4double fEdep = 0.;
    G4double fPrimaryEnergy = 0.;
    G4ThreeVector fPrimaryDir{0, 0, 1};
    G4ThreeVector fPrimaryPos{0, 0, 0};
    G4double fShowerZMin = 1e9;
    G4double fShowerZMax = -1e9;

    static constexpr G4int kProfileBins = 240;
    static constexpr G4double kProfileMaxMm = 1200.;
    G4double fLongProfile[kProfileBins]{};
};

}  // namespace B1

#endif
