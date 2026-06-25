#include "DetectorConstruction.hh"

#include "FilterTransmission.hh"
#include "OpticalMaterials.hh"
#include "RuntimeConfig.hh"
#include "SiPMSD.hh"
#include "ACDSD.hh"

#include <cstring>

#include "G4Box.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

namespace B1
{

namespace
{
void PlaceSiPMArray(G4LogicalVolume* parent, G4ThreeVector faceCenter, G4ThreeVector uDir,
                    G4ThreeVector vDir, int nRows, int nCols, G4LogicalVolume* logicSiPM,
                    G4RotationMatrix* rotation, G4int& copyNoBase)
{
  const G4double spacing = 50.0 * mm;
  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      const G4double offsetU = (i - (nRows - 1) / 2.0) * spacing;
      const G4double offsetV = (j - (nCols - 1) / 2.0) * spacing;
      const G4ThreeVector pos = faceCenter + offsetU * uDir + offsetV * vDir;
      new G4PVPlacement(rotation, pos, logicSiPM, "SiPM", parent, false, copyNoBase++, false);
    }
  }
}

G4bool AcdEnabledFromEnv()
{
  const char* env = std::getenv("G4ACD");
  if (!env) return true;
  return !(std::strcmp(env, "0") == 0 || std::strcmp(env, "off") == 0 ||
           std::strcmp(env, "OFF") == 0 || std::strcmp(env, "false") == 0);
}

void PlaceACDFace(G4LogicalVolume* parent, G4ThreeVector faceCenter, G4ThreeVector uDir,
                  G4ThreeVector vDir, G4LogicalVolume* logicTile, G4RotationMatrix* rotation,
                  G4int& copyNoBase, G4int gridSize)
{
  const G4double spacing = 60.0 * mm;
  for (int i = 0; i < gridSize; ++i) {
    for (int j = 0; j < gridSize; ++j) {
      const G4double offsetU = (i - (gridSize - 1) / 2.0) * spacing;
      const G4double offsetV = (j - (gridSize - 1) / 2.0) * spacing;
      const G4ThreeVector pos = faceCenter + offsetU * uDir + offsetV * vDir;
      new G4PVPlacement(rotation, pos, logicTile, "ACDTile", parent, false, copyNoBase++, false);
    }
  }
}
}  // namespace

G4bool DetectorConstruction::ACDEnabled() { return AcdEnabledFromEnv(); }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto* nist = G4NistManager::Instance();
  const G4bool checkOverlaps = false;

  const G4double world_size = 1.2 * m;
  auto* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  auto* solidWorld = new G4Box("World", 0.5 * world_size, 0.5 * world_size, 0.5 * world_size);
  auto* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  auto* physWorld =
    new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, checkOverlaps);

  const G4double telescope_size = 1.02 * m;
  auto* shell_mat = nist->FindOrBuildMaterial("G4_Al");
  auto* solidShell =
    new G4Box("TelescopeShell", 0.5 * telescope_size, 0.5 * telescope_size, 0.5 * telescope_size);
  auto* logicShell = new G4LogicalVolume(solidShell, shell_mat, "TelescopeShell");
  new G4PVPlacement(nullptr, G4ThreeVector(), logicShell, "TelescopeShell", logicWorld, false, 0,
                  checkOverlaps);

  const G4double lxe_size = 1.0 * m;
  auto* lxe_mat = new G4Material("LXe", 54.0, 131.29 * g / mole, 2.953 * g / cm3, kStateLiquid);
  OpticalMaterials::SetLXeProperties(lxe_mat);

  auto* solidLXe = new G4Box("LXe", 0.5 * lxe_size, 0.5 * lxe_size, 0.5 * lxe_size);
  auto* logicLXe = new G4LogicalVolume(solidLXe, lxe_mat, "LXe");
  new G4PVPlacement(nullptr, G4ThreeVector(), logicLXe, "LXe", logicShell, false, 0, checkOverlaps);

  const G4double sipm_sizeXY = 10.0 * mm;
  const G4double sipm_thickness = 1.0 * mm;
  auto* sipm_mat = nist->FindOrBuildMaterial("G4_Si");
  auto* solidSiPM =
    new G4Box("SiPM", 0.5 * sipm_sizeXY, 0.5 * sipm_sizeXY, 0.5 * sipm_thickness);
  fLogicSiPM = new G4LogicalVolume(solidSiPM, sipm_mat, "SiPM");

  SetupOpticalSurfaces(shell_mat, logicShell, fLogicSiPM);

  auto* rotPosZ = new G4RotationMatrix();
  rotPosZ->rotateX(CLHEP::pi);
  auto* rotPosX = new G4RotationMatrix();
  rotPosX->rotateY(CLHEP::halfpi);
  auto* rotNegX = new G4RotationMatrix();
  rotNegX->rotateY(-CLHEP::halfpi);
  auto* rotPosY = new G4RotationMatrix();
  rotPosY->rotateZ(CLHEP::pi);
  rotPosY->rotateX(CLHEP::halfpi);
  auto* rotNegY = new G4RotationMatrix();
  rotNegY->rotateX(-CLHEP::halfpi);

  const G4double half = 0.5 * lxe_size;
  const G4double offset = 0.5 * sipm_thickness + 0.01 * mm;
  G4int copyNo = 0;

  PlaceSiPMArray(logicLXe, {0, 0, +half - offset}, {1, 0, 0}, {0, 1, 0}, kSiPMGridSize,
                 kSiPMGridSize, fLogicSiPM, rotPosZ, copyNo);
  PlaceSiPMArray(logicLXe, {0, 0, -half + offset}, {1, 0, 0}, {0, -1, 0}, kSiPMGridSize,
                 kSiPMGridSize, fLogicSiPM, nullptr, copyNo);
  PlaceSiPMArray(logicLXe, {+half - offset, 0, 0}, {0, 1, 0}, {0, 0, -1}, kSiPMGridSize,
                 kSiPMGridSize, fLogicSiPM, rotPosX, copyNo);
  PlaceSiPMArray(logicLXe, {-half + offset, 0, 0}, {0, 1, 0}, {0, 0, 1}, kSiPMGridSize,
                 kSiPMGridSize, fLogicSiPM, rotNegX, copyNo);
  PlaceSiPMArray(logicLXe, {0, +half - offset, 0}, {1, 0, 0}, {0, 0, -1}, kSiPMGridSize,
                 kSiPMGridSize, fLogicSiPM, rotPosY, copyNo);
  PlaceSiPMArray(logicLXe, {0, -half + offset, 0}, {1, 0, 0}, {0, 0, 1}, kSiPMGridSize,
                 kSiPMGridSize, fLogicSiPM, rotNegY, copyNo);

  if (ACDEnabled()) {
    auto* acd_mat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    const G4double acd_thickness = 10.0 * mm;
    const G4double acd_tile = 58.0 * mm;
    const G4double shell_half = 0.5 * telescope_size;
    const G4double acd_gap = 2.0 * mm;
    const G4double acd_offset = shell_half + acd_gap + 0.5 * acd_thickness;

    auto* solidACD =
      new G4Box("ACDTileSolid", 0.5 * acd_tile, 0.5 * acd_tile, 0.5 * acd_thickness);
    fLogicACDTile = new G4LogicalVolume(solidACD, acd_mat, "ACDTile");
    fLogicACDTile->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.35)));

    G4int acdCopy = 0;
    PlaceACDFace(logicWorld, {0, 0, +acd_offset}, {1, 0, 0}, {0, 1, 0}, fLogicACDTile, nullptr,
                 acdCopy, kACDGridSize);
    PlaceACDFace(logicWorld, {0, 0, -acd_offset}, {1, 0, 0}, {0, -1, 0}, fLogicACDTile, nullptr,
                 acdCopy, kACDGridSize);
    PlaceACDFace(logicWorld, {+acd_offset, 0, 0}, {0, 1, 0}, {0, 0, -1}, fLogicACDTile, rotPosX,
                 acdCopy, kACDGridSize);
    PlaceACDFace(logicWorld, {-acd_offset, 0, 0}, {0, 1, 0}, {0, 0, 1}, fLogicACDTile, rotNegX,
                 acdCopy, kACDGridSize);
    PlaceACDFace(logicWorld, {0, +acd_offset, 0}, {1, 0, 0}, {0, 0, -1}, fLogicACDTile, rotPosY,
                 acdCopy, kACDGridSize);
    PlaceACDFace(logicWorld, {0, -acd_offset, 0}, {1, 0, 0}, {0, 0, 1}, fLogicACDTile, rotNegY,
                 acdCopy, kACDGridSize);
  }

  logicShell->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.2)));
  logicLXe->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3)));
  fLogicSiPM->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)));

  fScoringVolume = logicLXe;
  return physWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  if (fLogicSiPM) {
    auto* sipmSD = new SiPMSD("SiPMSD");
    G4SDManager::GetSDMpointer()->AddNewDetector(sipmSD);
    fLogicSiPM->SetSensitiveDetector(sipmSD);
  }
  if (fLogicACDTile) {
    auto* acdSD = new ACDSD("ACDSD");
    G4SDManager::GetSDMpointer()->AddNewDetector(acdSD);
    fLogicACDTile->SetSensitiveDetector(acdSD);
  }
}

void DetectorConstruction::SetupOpticalSurfaces(G4Material* /*shell_mat*/,
                                              G4LogicalVolume* logicShell,
                                              G4LogicalVolume* logicSiPM)
{
  const char* filterEnv = std::getenv("G4FILTER");
  const G4bool applyFilter =
    !(filterEnv && (std::strcmp(filterEnv, "0") == 0 || std::strcmp(filterEnv, "off") == 0 ||
                    std::strcmp(filterEnv, "OFF") == 0 || std::strcmp(filterEnv, "false") == 0));
  if (applyFilter) {
    const char* csv = std::getenv("G4FILTER_CSV");
    FilterTransmission::Instance().Initialize(
      csv ? G4String(csv) : G4String("../analysis/output/filter_transmission.csv"));
  }

  const G4int nE = 11;
  G4double energies[nE];
  G4double reflectivity[nE];
  G4double efficiency[nE];
  const G4double wavelengths[] = {600, 500, 400, 300, 250, 200, 190, 178, 175, 170, 160};
  constexpr G4double kPDE = 0.25;
  for (G4int i = 0; i < nE; ++i) {
    energies[i] = (1240.0 / wavelengths[i]) * eV;
    reflectivity[i] = 0.05;
    const G4double filterT =
      applyFilter ? FilterTransmission::Instance().Transmission(wavelengths[i]) : 1.0;
    efficiency[i] = kPDE * filterT;
  }

  auto* shellSurface = new G4OpticalSurface("AlShellSurface");
  shellSurface->SetType(dielectric_metal);
  shellSurface->SetModel(unified);
  shellSurface->SetFinish(polished);
  auto* shellMPT = new G4MaterialPropertiesTable();
  shellMPT->AddProperty("REFLECTIVITY", energies, reflectivity, nE);
  shellSurface->SetMaterialPropertiesTable(shellMPT);
  new G4LogicalSkinSurface("AlShellSkin", logicShell, shellSurface);

  auto* sipmSurface = new G4OpticalSurface("SiPMSurface");
  sipmSurface->SetType(dielectric_metal);
  sipmSurface->SetModel(unified);
  sipmSurface->SetFinish(polished);
  auto* sipmMPT = new G4MaterialPropertiesTable();
  sipmMPT->AddProperty("EFFICIENCY", energies, efficiency, nE);
  sipmSurface->SetMaterialPropertiesTable(sipmMPT);
  new G4LogicalSkinSurface("SiPMSkin", logicSiPM, sipmSurface);
}

}  // namespace B1
