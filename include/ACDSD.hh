#ifndef ACDSD_HH
#define ACDSD_HH

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include <map>

class ACDSD : public G4VSensitiveDetector
{
  public:
    explicit ACDSD(const G4String& name);
    ~ACDSD() override = default;

    void Initialize(G4HCofThisEvent* /*hce*/) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

    G4int GetNTilesHit() const;
    G4double GetTotalEdepMeV() const { return fTotalEdep / CLHEP::MeV; }
    G4bool GetVeto() const { return GetNTilesHit() > 0; }

    static constexpr G4double kTileThresholdMeV = 0.15;

  private:
    std::map<G4int, G4double> fTileEdep;
    G4double fTotalEdep = 0.;
};

#endif
