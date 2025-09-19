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

#include "SteppingAction.hh"

#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include <cstdio>  // for std::ifstream
#include <fstream>
#include <iostream>

namespace B1
{

SteppingAction::SteppingAction(EventAction* eventAction)
  : fEventAction(eventAction), fScoringVolume(nullptr)
{}

// ÿ��ִ��
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // �ӳ���scoring volume������ÿ������
  if (!fScoringVolume) {
    auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // �õ�ǰstep���߼���
  G4LogicalVolume* volume =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // ֻ��SiPM�߼����еĹ��Ӵ���
  if (volume->GetName() != "SiPM") return;

  // ֻ�ܹ���
  G4Track* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;

  // �¼�ID
  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  // ��ȡ��Ϣ
  G4StepPoint* postPoint = step->GetPostStepPoint();
  G4ThreeVector pos = postPoint->GetPosition();
  G4double time = postPoint->GetGlobalTime();

  // ===== Remove all file writing to cherenkov_photons.csv =====
  // (No file output here)

  // ͳ
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
}

}  // namespace B1
