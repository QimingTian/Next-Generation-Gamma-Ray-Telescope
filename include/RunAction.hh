#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

namespace B1
{

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    void SetOutputTag(const G4String& tag) { fOutputTag = tag; }

  private:
    G4String fOutputTag = "default";
    G4String fOutputPrefix;
};

}  // namespace B1

#endif
