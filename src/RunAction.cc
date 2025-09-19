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
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

// RunAction.cc
#include "RunAction.hh"

#include "G4AccumulableManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SiPMSD.hh"

#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <mutex>

namespace B1
{

RunAction::RunAction() : G4UserRunAction(), fEdep(0.), fEdep2(0.)
{
  // �Զ��嵥λ
  const G4double milligray = 1.e-3 * gray;
  const G4double microgray = 1.e-6 * gray;
  const G4double nanogray = 1.e-9 * gray;
  const G4double picogray = 1.e-12 * gray;

  new G4UnitDefinition("milligray", "milliGy", "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy", "Dose", microgray);
  new G4UnitDefinition("nanogray", "nanoGy", "Dose", nanogray);
  new G4UnitDefinition("picogray", "picoGy", "Dose", picogray);

  // ע��accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  G4AccumulableManager::Instance()->Reset();
  // No file writing here; master will write at end of run
  G4cout << "[RunAction] Ready to collect event summaries." << G4endl;
}

void RunAction::AddEdep(G4double edep)
{
  fEdep += edep;
  fEdep2 += edep * edep;
}

void RunAction::AddEventSummary(int eventID, const std::map<std::string, int>& processCounts)
{
  std::lock_guard<std::mutex> lock(fSummaryMutex);
  fEventSummaries.emplace_back(eventID, processCounts);
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4AccumulableManager::Instance()->Merge();

  G4double edep = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();

  G4double rms = edep2 - edep * edep / nofEvents;
  if (rms > 0.)
    rms = std::sqrt(rms);
  else
    rms = 0.;

  const auto detConstruction = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep / mass;
  G4double rmsDose = rms / mass;

  const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction) {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy, "Energy");
  }

  if (IsMaster()) {
    G4cout << G4endl << "--------------------End of Global Run-----------------------";
    // Write the CSV file from all collected event summaries
    std::ofstream outFile("cherenkov_photons.csv", std::ios::trunc);
    outFile << "EventID,Process,Count\n";
    for (const auto& entry : fEventSummaries) {
      int eventID = entry.first;
      const auto& processCounts = entry.second;
      for (const auto& proc : processCounts) {
        outFile << eventID << "," << proc.first << "," << proc.second << "\n";
      }
    }
    outFile.close();
    G4cout << "[RunAction] Wrote cherenkov_photons.csv with all event summaries." << G4endl;
  }
  else {
    G4cout << G4endl << "--------------------End of Local Run------------------------";
  }

  G4cout << G4endl << " The run consists of " << nofEvents << " " << runCondition << G4endl
         << " Cumulated dose per run, in scoring volume : " << G4BestUnit(dose, "Dose")
         << " rms = " << G4BestUnit(rmsDose, "Dose") << G4endl
         << "------------------------------------------------------------" << G4endl;
}

}  // namespace B1
