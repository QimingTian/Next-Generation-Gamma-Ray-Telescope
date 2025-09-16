//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

//
// DetectorConstruction.cc
//

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "SiPMSD.hh" 

namespace B1
{

void PlaceSiPMArray(G4LogicalVolume* parent, G4ThreeVector faceCenter, G4ThreeVector uDir,
                    G4ThreeVector vDir, int nRows, int nCols, G4LogicalVolume* logicSiPM,
                    G4RotationMatrix* rotation)
{
  G4double spacing = 50.0 * mm;
  int idx = 0;

  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      G4double offsetU = (i - (nRows - 1) / 2.0) * spacing;
      G4double offsetV = (j - (nCols - 1) / 2.0) * spacing;
      G4ThreeVector pos = faceCenter + offsetU * uDir + offsetV * vDir;

      new G4PVPlacement(rotation, pos, logicSiPM, "SiPM", parent, false, idx++, true);
    }
  }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  // -------- World --------
  G4double world_size = 1.2 * m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld = new G4Box("World", 0.5 * world_size, 0.5 * world_size, 0.5 * world_size);
  auto logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");

  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false,
                                     0, checkOverlaps);

  // -------- Shell --------
  G4double telescope_size = 1.0 * m;
  G4Material* shell_mat = nist->FindOrBuildMaterial("G4_Al");

  auto solidShell =
    new G4Box("TelescopeShell", 0.5 * telescope_size, 0.5 * telescope_size, 0.5 * telescope_size);
  auto logicShell = new G4LogicalVolume(solidShell, shell_mat, "TelescopeShell");

  // ---- Add optical surface for inner shell reflectivity ----
  auto shellSurface = new G4OpticalSurface("ShellSurface");
  shellSurface->SetType(dielectric_metal);
  shellSurface->SetModel(unified);
  shellSurface->SetFinish(polished);

  auto shellMPT = new G4MaterialPropertiesTable();
  G4double photonEnergy[2] = {1.5 * eV, 6.2 * eV};
  G4double reflectivity[2] = {0.005, 0.005}; // 0.5% reflectivity
  shellMPT->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, 2);
  shellSurface->SetMaterialPropertiesTable(shellMPT);

  new G4LogicalSkinSurface("ShellSkinSurface", logicShell, shellSurface);

  new G4PVPlacement(nullptr, G4ThreeVector(), logicShell, "TelescopeShell", logicWorld, false, 0,
                    checkOverlaps);

  // -------- Liquid Xenon --------
  G4double gap = 1.0 * cm;
  G4double lxe_size = telescope_size - 2 * gap;
  G4Material* lxe_mat = nist->FindOrBuildMaterial("G4_lXe");

  G4MaterialPropertiesTable* LXeMPT = new G4MaterialPropertiesTable();
  // Set wavelength-dependent refractive index for LXe
  const G4int NUM_RINDEX = 11;
  G4double wavelength_nm[NUM_RINDEX] = {160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210};
  G4double photonEnergy[NUM_RINDEX];
  for (int i = 0; i < NUM_RINDEX; ++i) {
    photonEnergy[i] = 1239.841984 / wavelength_nm[i]; // E[eV] = 1239.841984 / lambda[nm]
  }
  G4double refractiveIndex[NUM_RINDEX] = {2.10, 1.90, 1.82, 1.74, 1.70, 1.66, 1.62, 1.60, 1.58, 1.56, 1.54};
  LXeMPT->AddProperty("RINDEX", photonEnergy, refractiveIndex, NUM_RINDEX);
  LXeMPT->AddConstProperty("SCINTILLATIONYIELD", 25000. / MeV);
  LXeMPT->AddConstProperty("FASTTIMECONSTANT", 4.3 * ns, true);
  LXeMPT->AddConstProperty("SLOWTIMECONSTANT", 22 * ns, true);
  lxe_mat->SetMaterialPropertiesTable(LXeMPT);

  auto solidLXe = new G4Box("LXe", 0.5 * lxe_size, 0.5 * lxe_size, 0.5 * lxe_size);
  auto logicLXe = new G4LogicalVolume(solidLXe, lxe_mat, "LXe");

  new G4PVPlacement(nullptr, G4ThreeVector(), logicLXe, "LXe", logicShell, false, 0, checkOverlaps);

  // -------- SiPM logical volume --------
  G4double sipm_sizeXY = 10.0 * mm;  // 1cm �� 1cm
  G4double sipm_thickness = 1.0 * mm;

  G4Material* sipm_mat = nist->FindOrBuildMaterial("G4_Si");
  auto solidSiPM = new G4Box("SiPM", 0.5 * sipm_sizeXY, 0.5 * sipm_sizeXY, 0.5 * sipm_thickness);
  auto logicSiPM = new G4LogicalVolume(solidSiPM, sipm_mat, "SiPM");

  auto sipmSD = new SiPMSD("SiPMSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(sipmSD);
  logicSiPM->SetSensitiveDetector(sipmSD);

  G4RotationMatrix* rotPosZ = new G4RotationMatrix();
  rotPosZ->rotateX(CLHEP::pi);  

  G4RotationMatrix* rotNegZ = nullptr;

  G4RotationMatrix* rotPosX = new G4RotationMatrix();
  rotPosX->rotateY(CLHEP::halfpi);  //

  G4RotationMatrix* rotNegX = new G4RotationMatrix();
  rotNegX->rotateY(-CLHEP::halfpi);  // 

  G4RotationMatrix* rotPosY = new G4RotationMatrix();
  rotPosY->rotateZ(CLHEP::pi);
  rotPosY->rotateX(CLHEP::halfpi);

  G4RotationMatrix* rotNegY = new G4RotationMatrix();
  rotNegY->rotateX(-CLHEP::halfpi);

  G4double half = 0.5 * lxe_size;
  G4double offset = 0.5 * sipm_thickness + 0.01 * mm;

  PlaceSiPMArray(logicLXe, {0, 0, +half - offset}, {1, 0, 0}, {0, 1, 0}, 18, 18, logicSiPM,
                 rotPosZ);  // +Z
  PlaceSiPMArray(logicLXe, {0, 0, -half + offset}, {1, 0, 0}, {0, -1, 0}, 18, 18, logicSiPM,
                 rotNegZ);  // -Z
  PlaceSiPMArray(logicLXe, {+half - offset, 0, 0}, {0, 1, 0}, {0, 0, -1}, 18, 18, logicSiPM,
                 rotPosX);  // +X
  PlaceSiPMArray(logicLXe, {-half + offset, 0, 0}, {0, 1, 0}, {0, 0, 1}, 18, 18, logicSiPM,
                 rotNegX);  // -X
  PlaceSiPMArray(logicLXe, {0, +half - offset, 0}, {1, 0, 0}, {0, 0, -1}, 18, 18, logicSiPM,
                 rotPosY);  // +Y
  PlaceSiPMArray(logicLXe, {0, -half + offset, 0}, {1, 0, 0}, {0, 0, 1}, 18, 18, logicSiPM,
                 rotNegY);  // -Y

  //
  auto shellVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.2));
  shellVis->SetForceSolid(true);
  logicShell->SetVisAttributes(shellVis);

  auto lxeVis = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));
  lxeVis->SetForceSolid(true);
  logicLXe->SetVisAttributes(lxeVis);

  auto sipmVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  sipmVis->SetForceSolid(true);
  logicSiPM->SetVisAttributes(sipmVis);

  fScoringVolume = logicLXe;

  return physWorld;
}

}  // namespace B1
