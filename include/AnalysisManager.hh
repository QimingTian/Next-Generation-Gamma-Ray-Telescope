#ifndef B1AnalysisManager_h
#define B1AnalysisManager_h 1

#include "globals.hh"

namespace B1
{

/// Thread-safe output via G4AnalysisManager (CSV ntuples).
class AnalysisManager
{
  public:
    static AnalysisManager* Instance();
    static void DeleteInstance();

    void Book();
    G4bool IsBooked() const { return fBooked; }
    void Open(const G4String& fileName);
    void Close();

    void FillPhoton(G4int eventID, G4int sipmID, G4double x, G4double y, G4double z,
                    G4double t, G4double wavelength_nm, const G4String& process);
    void FillEvent(G4int eventID, G4double ePrimary, G4double dirX, G4double dirY,
                   G4double dirZ, G4int nCherenkov, G4int nScint, G4double edep,
                   G4double showerLength, G4double l90_mm, G4int truncated,
                   G4int acdNTiles = 0, G4double acdEdepMeV = 0., G4int acdVeto = 0);
    void FillSiPMSummary(G4int eventID, G4int sipmID, G4int photonCount);

  private:
    AnalysisManager() = default;
    ~AnalysisManager() = default;

    static AnalysisManager* fInstance;

    G4int fPhotonNtuple = -1;
    G4int fEventNtuple = -1;
    G4int fSiPMNtuple = -1;
    G4bool fBooked = false;
};

}  // namespace B1

#endif
