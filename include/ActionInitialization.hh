#ifndef B1ActionInitialization_h
#define B1ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

namespace B1
{

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization() = default;
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

    void SetOutputTag(const G4String& tag) { fOutputTag = tag; }

  private:
    G4String fOutputTag = "default";
};

}  // namespace B1

#endif
