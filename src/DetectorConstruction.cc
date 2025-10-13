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

#include "SiPMSD.hh"  // ����SiPM SDͷ�ļ�

namespace B1
{

// ====== �޸İ棺PlaceSiPMArray �������ټ�����ת����ֻ������� ======
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

      new G4PVPlacement(rotation, pos, logicSiPM, "SiPM", parent, false, idx++, false);
    }
  }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = false;  // Disable overlap checking for faster startup

  // -------- World (Space/Vacuum Environment) --------
  G4double world_size = 1.2 * m;
  // Use vacuum instead of air to simulate space environment
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  auto solidWorld = new G4Box("World", 0.5 * world_size, 0.5 * world_size, 0.5 * world_size);
  auto logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");

  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false,
                                     0, checkOverlaps);

  // -------- Shell --------
  G4double telescope_size = 1.02 * m;
  G4Material* shell_mat = nist->FindOrBuildMaterial("G4_Al");

  // Set up aluminum shell optical properties
  G4MaterialPropertiesTable* AlMPT = new G4MaterialPropertiesTable();
  // Use same energy points as LXe for consistency
  const G4int AlNUM = 13;
  G4double alPhotonEnergy[AlNUM] = {
    1240.0/600.0 * eV,  // 600 nm -> 2.07 eV
    1240.0/400.0 * eV,  // 400 nm -> 3.10 eV
    1240.0/210.0 * eV,  // 210 nm -> 5.90 eV
    1240.0/205.0 * eV,  // 205 nm -> 6.05 eV
    1240.0/200.0 * eV,  // 200 nm -> 6.20 eV
    1240.0/195.0 * eV,  // 195 nm -> 6.36 eV
    1240.0/190.0 * eV,  // 190 nm -> 6.53 eV
    1240.0/185.0 * eV,  // 185 nm -> 6.70 eV
    1240.0/180.0 * eV,  // 180 nm -> 6.89 eV
    1240.0/175.0 * eV,  // 175 nm -> 7.09 eV
    1240.0/170.0 * eV,  // 170 nm -> 7.29 eV
    1240.0/165.0 * eV,  // 165 nm -> 7.52 eV
    1240.0/160.0 * eV   // 160 nm -> 7.75 eV
  };
  // 5% reflectivity for all wavelengths
  G4double alReflectivity[AlNUM] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  AlMPT->AddProperty("REFLECTIVITY", alPhotonEnergy, alReflectivity, AlNUM);
  shell_mat->SetMaterialPropertiesTable(AlMPT);

  auto solidShell =
    new G4Box("TelescopeShell", 0.5 * telescope_size, 0.5 * telescope_size, 0.5 * telescope_size);
  auto logicShell = new G4LogicalVolume(solidShell, shell_mat, "TelescopeShell");

  new G4PVPlacement(nullptr, G4ThreeVector(), logicShell, "TelescopeShell", logicWorld, false, 0,
                    checkOverlaps);

  // -------- Liquid Xenon --------
  G4double gap = 1.0 * cm;
  G4double lxe_size = 1.0 * m;  // Fixed at 100cm
  G4Material* lxe_mat = nist->FindOrBuildMaterial("G4_lXe");

  // ���ù�ѧ����
  G4MaterialPropertiesTable* LXeMPT = new G4MaterialPropertiesTable();
  const G4int NUM = 10;
  // Wavelength-dependent refractive index and absorption length data
  // Wavelengths: 160, 170, 175, 178, 200, 250, 300, 400, 500, 600 nm
  // Convert wavelengths to photon energies: E = hc/λ (where hc = 1240 eV·nm)
  // Energies must be in increasing order
  G4double photonEnergy[NUM] = {
    1240.0/600.0 * eV,  // 600 nm -> 2.07 eV
    1240.0/500.0 * eV,  // 500 nm -> 2.48 eV
    1240.0/400.0 * eV,  // 400 nm -> 3.10 eV
    1240.0/300.0 * eV,  // 300 nm -> 4.13 eV
    1240.0/250.0 * eV,  // 250 nm -> 4.96 eV
    1240.0/200.0 * eV,  // 200 nm -> 6.20 eV
    1240.0/178.0 * eV,  // 178 nm -> 6.97 eV
    1240.0/175.0 * eV,  // 175 nm -> 7.09 eV
    1240.0/170.0 * eV,  // 170 nm -> 7.29 eV
    1240.0/160.0 * eV   // 160 nm -> 7.75 eV
  };
  G4double refractiveIndex[NUM] = {1.4, 1.4, 1.4, 1.4, 1.4, 1.58, 1.60, 1.70, 1.74, 2.10};
  // Absorption length: wavelength-dependent values (meters) - your exact data
  G4double absorptionLength[NUM] = {
    100.0 * m,  // 600 nm -> 100 m
    100.0 * m,  // 500 nm -> 100 m
    100.0 * m,  // 400 nm -> 100 m
    80.0 * m,   // 300 nm -> 80 m
    70.0 * m,   // 250 nm -> 70 m
    60.0 * m,   // 200 nm -> 60 m
    50.0 * m,   // 178 nm -> 50 m
    50.0 * m,   // 175 nm -> 50 m
    20.0 * m,   // 170 nm -> 20 m
    1.0 * m     // 160 nm -> 1 m
  };
  
  // Rayleigh scattering length: wavelength-dependent values (meters) - based on your exact data
  // Data points: 150nm->0cm, 200nm->150cm, 250nm->570cm, 300nm->1400cm, 350nm->2800cm, 400nm->5000cm, 450nm->9000cm
  // Interpolated for the energy points used
  G4double rayleighLength[NUM] = {
    280.0 * m,  // 600 nm -> extrapolated (very long)
    135.0 * m,   // 500 nm -> extrapolated  
    50.0 * m,   // 400 nm -> 5000 cm = 50 m
    14.0 * m,   // 300 nm -> 1400 cm = 14 m
    5.7 * m,    // 250 nm -> 570 cm = 5.7 m
    1.5 * m,    // 200 nm -> 150 cm = 1.5 m
    0.8 * m,    // 178 nm -> interpolated between 200nm(150cm) and 250nm(570cm)
    0.6 * m,    // 175 nm -> interpolated between 200nm(150cm) and 250nm(570cm)
    0.3 * m,    // 170 nm -> interpolated between 150nm(0cm) and 200nm(150cm)
    0.1 * m     // 160 nm -> interpolated between 150nm(0cm) and 200nm(150cm)
  };
  
  LXeMPT->AddProperty("RINDEX", photonEnergy, refractiveIndex, NUM);
  LXeMPT->AddProperty("ABSLENGTH", photonEnergy, absorptionLength, NUM);
  LXeMPT->AddProperty("RAYLEIGH", photonEnergy, rayleighLength, NUM);
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

  // ����Ϊ����̽����
  auto sipmSD = new SiPMSD("SiPMSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(sipmSD);
  logicSiPM->SetSensitiveDetector(sipmSD);

// -------- �̶���ת������ --------
  // +Z�� SiPM ���� -Z��Ĭ��SiPM����+Z����Ҫ��X��ת180�ȣ�
  G4RotationMatrix* rotPosZ = new G4RotationMatrix();
  rotPosZ->rotateX(CLHEP::pi);  // 180����ת

  // -Z�� SiPM ���� +Z��Ĭ�ϣ�������ת��
  G4RotationMatrix* rotNegZ = nullptr;

  // +X�� SiPM ���� -X����Y��ת +90�ȣ�
  G4RotationMatrix* rotPosX = new G4RotationMatrix();
  rotPosX->rotateY(CLHEP::halfpi);  // +90��

  // -X�� SiPM ���� +X����Y��ת -90�ȣ�
  G4RotationMatrix* rotNegX = new G4RotationMatrix();
  rotNegX->rotateY(-CLHEP::halfpi);  // -90��

  // +Y�� SiPM ���� -Y����X��ת+90�ȣ�����Z��ת180�ȣ�
  G4RotationMatrix* rotPosY = new G4RotationMatrix();
  rotPosY->rotateZ(CLHEP::pi);
  rotPosY->rotateX(CLHEP::halfpi);

  // -Y�� SiPM ���� +Y����X��ת-90�ȣ�
  G4RotationMatrix* rotNegY = new G4RotationMatrix();
  rotNegY->rotateX(-CLHEP::halfpi);

  // -------- ������ 6 ���� --------
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

  // -------- ���ӻ� --------
  auto shellVis = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.2));
  shellVis->SetForceSolid(true);
  logicShell->SetVisAttributes(shellVis);

  auto lxeVis = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3));
  lxeVis->SetForceSolid(true);
  logicLXe->SetVisAttributes(lxeVis);

  auto sipmVis = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  sipmVis->SetForceSolid(true);
  logicSiPM->SetVisAttributes(sipmVis);

  // �趨 scoring volume
  fScoringVolume = logicLXe;

  return physWorld;
}

}  // namespace B1
