#ifndef B1PhysicsList_h
#define B1PhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList() override = default;

    void SetPhysicsHadronic(G4bool enable) { fHadronic = enable; }

  private:
    G4bool fHadronic = false;
};

#endif
