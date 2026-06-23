#include "RunAction.hh"

#include "AnalysisManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <cstdio>

namespace
{
void EnsureDirectory(const G4String& path)
{
  mkdir(path.c_str(), 0755);
}

G4String GetGitCommit()
{
  FILE* pipe = popen("git rev-parse --short HEAD 2>/dev/null", "r");
  if (!pipe) return "unknown";
  char buffer[128];
  G4String result = "unknown";
  if (fgets(buffer, sizeof(buffer), pipe)) {
    result = buffer;
    while (!result.empty() && (result.back() == '\n' || result.back() == '\r')) {
      result.pop_back();
    }
  }
  pclose(pipe);
  return result;
}

void WriteRunMeta(const G4String& prefix, const G4Run* run)
{
  std::ostringstream metaPath;
  metaPath << prefix << "_meta.json";
  std::ofstream out(metaPath.str());
  if (!out) return;

  const char* seedEnv = std::getenv("G4SEED");
  const char* macroEnv = std::getenv("G4MACRO");

  out << "{\n"
      << "  \"run_id\": " << run->GetRunID() << ",\n"
      << "  \"events\": " << run->GetNumberOfEvent() << ",\n"
      << "  \"geant4_version\": \"11.3.2\",\n"
      << "  \"random_seed\": \"" << (seedEnv ? seedEnv : "123456789") << "\",\n"
      << "  \"macro\": \"" << (macroEnv ? macroEnv : "") << "\",\n"
      << "  \"git_commit\": \"" << GetGitCommit() << "\",\n"
      << "  \"timestamp\": \"" << std::time(nullptr) << "\"\n"
      << "}\n";
}
}  // namespace

namespace B1
{

RunAction::RunAction() = default;

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  EnsureDirectory("data");

  if (!AnalysisManager::Instance()->IsBooked()) {
    AnalysisManager::Instance()->Book();
  }

  std::ostringstream oss;
  oss << "data/" << fOutputTag << "_run" << run->GetRunID();
  fOutputPrefix = oss.str();
  AnalysisManager::Instance()->Open(fOutputPrefix);

  G4cout << "### Run " << run->GetRunID() << " start, output: " << fOutputPrefix << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  AnalysisManager::Instance()->Close();
  WriteRunMeta(fOutputPrefix, run);

  G4cout << "### Run " << run->GetRunID() << " end. Events: " << run->GetNumberOfEvent()
         << G4endl;
}

}  // namespace B1
