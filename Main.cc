#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "ActionInitialization.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace B1;

int main(int argc, char** argv)
{
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  G4SteppingVerbose::UseBestUnit(4);

  const char* seedEnv = std::getenv("G4SEED");
  long seed = seedEnv ? std::atol(seedEnv) : 123456789L;
  CLHEP::HepRandom::setTheSeed(seed);

  auto runType = ui ? G4RunManagerType::Default : G4RunManagerType::Serial;
  auto* runManager = G4RunManagerFactory::CreateRunManager(runType);

  runManager->SetUserInitialization(new DetectorConstruction());
  auto* physicsList = new PhysicsList();
  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
  auto* actionInit = new ActionInitialization();
  const char* outputTag = std::getenv("G4OUTPUT_TAG");
  if (outputTag) {
    actionInit->SetOutputTag(outputTag);
  }
  runManager->SetUserInitialization(actionInit);

  G4VisExecutive* visManager = nullptr;
  if (ui) {
    visManager = new G4VisExecutive(argc, argv);
    visManager->Initialize();
  }

  auto* UImanager = G4UImanager::GetUIpointer();

  if (!ui) {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  } else {
    ui->SessionStart();
    delete ui;
  }

  AnalysisManager::DeleteInstance();
  delete visManager;
  delete runManager;
  return 0;
}
