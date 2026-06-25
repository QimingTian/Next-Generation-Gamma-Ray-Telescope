#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

namespace B1
{

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    static constexpr G4int kSiPMGridSize = 18;
    static constexpr G4int kSiPMsPerFace = kSiPMGridSize * kSiPMGridSize;
    static constexpr G4int kTotalSiPMs = kSiPMsPerFace * 6;

    // Fermi-LAT-style ACD: 16×16 scintillator tiles per face, 1 cm thick.
    static constexpr G4int kACDGridSize = 16;
    static constexpr G4int kACDTilesPerFace = kACDGridSize * kACDGridSize;
    static constexpr G4int kTotalACDTiles = kACDTilesPerFace * 6;

  private:
    void SetupOpticalSurfaces(G4Material* shell_mat, G4LogicalVolume* logicShell,
                              G4LogicalVolume* logicSiPM);
    static G4bool ACDEnabled();

    G4LogicalVolume* fScoringVolume = nullptr;
    G4LogicalVolume* fLogicSiPM = nullptr;
    G4LogicalVolume* fLogicACDTile = nullptr;
};

}  // namespace B1

#endif
