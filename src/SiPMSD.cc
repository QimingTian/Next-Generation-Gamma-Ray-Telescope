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
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SiPMSD.hh"

#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

SiPMSD::SiPMSD(const G4String& name) : G4VSensitiveDetector(name)
{
  // ע��̽������ SD ������
  G4SDManager::GetSDMpointer()->AddNewDetector(this);
}

SiPMSD::~SiPMSD() {}

void SiPMSD::Initialize(G4HCofThisEvent* /*hce*/)
{
  fHits.clear();  // ����ϸ��¼��� hits
}

G4bool SiPMSD::ProcessHits(G4Step* step, G4TouchableHistory* /*history*/)
{
  G4cout << "[DEBUG] SiPMSD::ProcessHits called!" << G4endl;

  auto track = step->GetTrack();

  // Only record photons with wavelength > 190 nm (E <= 6.52 eV)
  G4double photonEnergy = track->GetTotalEnergy() / eV; // in eV
  if (photonEnergy > 6.52) return false;

  // �ж��ǲ����� Cerenkov ���̲����Ĺ���
  const G4VProcess* creationProcess = track->GetCreatorProcess();
  G4String processName = creationProcess ? creationProcess->GetProcessName() : "Unknown";

  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4ThreeVector pos = preStepPoint->GetPosition();
  G4double time = preStepPoint->GetGlobalTime();

  fHits.push_back({pos, time, processName});
  return true;
}

void SiPMSD::EndOfEvent(G4HCofThisEvent* /*hce*/)
{
  G4cout << "[SiPMSD] Total Cherenkov hits in this event: " << fHits.size() << G4endl;

  // ��ѡ���ԣ����ÿ�� hit ����Ϣ
  /*
  for (const auto& hit : fHits)
  {
    G4cout << "  Hit at " << hit.position
           << ", time = " << hit.time / ns << " ns" << G4endl;
  }
  */
}