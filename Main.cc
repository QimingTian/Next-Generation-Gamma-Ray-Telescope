#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4Threading.hh"
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
#include <thread>

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

  // Batch macros: MT enabled unless G4MT=0. Optical+MT can segfault on some builds;
  // use parallel serial jobs (run_paper_quality_parallel.sh) as the default fast path.
  G4RunManagerType runType = G4RunManagerType::Serial;
  if (ui) {
    runType = G4RunManagerType::Default;
  } else {
    const char* mtEnv = std::getenv("G4MT");
    if (mtEnv && std::string(mtEnv) != "0") {
      runType = G4RunManagerType::MT;
    }
  }

  auto* runManager = G4RunManagerFactory::CreateRunManager(runType);

  if (runType == G4RunManagerType::MT) {
    G4int nThreads = static_cast<G4int>(std::thread::hardware_concurrency());
    if (nThreads < 1) {
      nThreads = 1;
    }
    if (const char* threadsEnv = std::getenv("G4NUM_THREADS")) {
      nThreads = static_cast<G4int>(std::atol(threadsEnv));
    }
    if (nThreads < 1) {
      nThreads = 1;
    }
    runManager->SetNumberOfThreads(nThreads);
    G4cout << "### Multithreaded mode: " << nThreads << " threads" << G4endl;
  }

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
