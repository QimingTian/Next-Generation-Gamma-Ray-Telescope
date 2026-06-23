#ifndef SIPMSD_HH
#define SIPMSD_HH

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <map>
#include <vector>

class SiPMSD : public G4VSensitiveDetector
{
  public:
    struct HitData
    {
      G4int sipmID = -1;
      G4ThreeVector position;
      G4double time = 0.;
      G4double wavelength_nm = 0.;
      G4String processName;
    };

    explicit SiPMSD(const G4String& name);
    ~SiPMSD() override = default;

    void Initialize(G4HCofThisEvent* hce) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;

    const std::vector<HitData>& GetHits() const { return fHits; }
    const std::map<G4int, G4int>& GetSiPMCounts() const { return fSiPMCounts; }
    G4int GetNCherenkov() const { return fNCherenkov; }
    G4int GetNScint() const { return fNScint; }

  private:
    std::vector<HitData> fHits;
    std::map<G4int, G4int> fSiPMCounts;
    G4int fNCherenkov = 0;
    G4int fNScint = 0;
    G4bool fStoreHits = false;
    G4bool fApplyFilter = false;
};

#endif
