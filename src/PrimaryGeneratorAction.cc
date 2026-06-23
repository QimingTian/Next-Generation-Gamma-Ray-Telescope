#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

namespace B1
{

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fGPS = new G4GeneralParticleSource();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  fGPS->GeneratePrimaryVertex(event);
}

}  // namespace B1
