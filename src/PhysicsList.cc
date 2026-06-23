#include "PhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "G4SystemOfUnits.hh"

#ifdef WITH_HADRONIC
#include "G4HadronPhysicsFTFP_BERT.hh"
#endif

PhysicsList::PhysicsList()
{
  RegisterPhysics(new G4EmStandardPhysics_option4());

#ifdef WITH_HADRONIC
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
#endif

  auto* optical = new G4OpticalPhysics();
  RegisterPhysics(optical);

  auto* opParams = G4OpticalParameters::Instance();
  opParams->SetScintTrackSecondariesFirst(true);
  opParams->SetCerenkovTrackSecondariesFirst(true);
  opParams->SetCerenkovMaxPhotonsPerStep(300);
  opParams->SetScintByParticleType(true);

  RegisterPhysics(new G4DecayPhysics());
}
