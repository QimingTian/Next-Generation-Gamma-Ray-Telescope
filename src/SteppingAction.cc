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

// 每步执行
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // 延迟拿scoring volume，避免每步都查
  if (!fScoringVolume) {
    auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // 拿当前step的逻辑体
  G4LogicalVolume* volume =
    step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  // 只对SiPM逻辑体中的光子处理
  if (volume->GetName() != "SiPM") return;

  // 只管光子
  G4Track* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) return;

  // 事件ID
  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  // 读取信息
  G4StepPoint* postPoint = step->GetPostStepPoint();
  G4ThreeVector pos = postPoint->GetPosition();
  G4double time = postPoint->GetGlobalTime();

  // ===== 写文件逻辑：先判断文件是否存在 =====
  std::ifstream checkFile("cherenkov_photons.csv");
  bool fileExists = checkFile.good();
  checkFile.close();

  std::ofstream outFile("cherenkov_photons.csv", std::ios::app);
  if (!outFile.is_open()) {
    G4cerr << "[SteppingAction] ERROR: Can't open cherenkov_photons.csv" << G4endl;
    return;
  }

  // 如果文件不存在，写入表头
  if (!fileExists) {
    outFile << "EventID,PosX_mm,PosY_mm,PosZ_mm,Time_ns\n";
  }

  // 写入数据
  outFile << eventID << "," << pos.x() / mm << "," << pos.y() / mm << "," << pos.z() / mm << ","
          << time / ns << "\n";
  outFile.close();

  // 能量沉积统计
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
}

}  // namespace B1