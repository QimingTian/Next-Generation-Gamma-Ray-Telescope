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
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

// EventAction.cc

#include "EventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"
#include "SiPMSD.hh"

#include <fstream>
#include <iostream>
#include <map> // Added for process counting

namespace B1
{

EventAction::EventAction(RunAction* runAction)
  : fRunAction(runAction), fEdep(0.), fIsFirstEvent(true)
{}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
  fProcessCounts.clear();
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  fRunAction->AddEdep(fEdep);

  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  auto sipmSD = (SiPMSD*)sdManager->FindSensitiveDetector("SiPMSD");

  if (!sipmSD) {
    G4cerr << "[EventAction] Error: Can't find SiPMSD!" << G4endl;
    return;
  }

  const auto& hits = sipmSD->GetHits();

  // Count number of hits for each process type
  for (const auto& hit : hits) {
    fProcessCounts[hit.processName]++;
  }

  // Store this event's summary in RunAction (thread-safe)
  fRunAction->AddEventSummary(event->GetEventID(), fProcessCounts);
}

}  // namespace B1
