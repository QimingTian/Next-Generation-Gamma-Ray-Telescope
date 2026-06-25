#include "EventAction.hh"

#include "AnalysisManager.hh"
#include "RunAction.hh"
#include "RuntimeConfig.hh"
#include "SiPMSD.hh"
#include "ACDSD.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

EventAction::EventAction(RunAction* runAction) : fRunAction(runAction) {}

void EventAction::BeginOfEventAction(const G4Event* event)
{
  fEdep = 0.;
  fShowerZMin = 1e9;
  fShowerZMax = -1e9;
  fPrimaryEnergy = 0.;
  fPrimaryDir = G4ThreeVector(0, 0, 1);
  fPrimaryPos = G4ThreeVector(0, 0, 0);
  for (G4int i = 0; i < kProfileBins; ++i) {
    fLongProfile[i] = 0.;
  }

  const auto* vertex = event->GetPrimaryVertex();
  if (vertex) {
    fPrimaryPos = vertex->GetPosition();
    const auto* particle = vertex->GetPrimary();
    if (particle) {
      fPrimaryEnergy = particle->GetKineticEnergy();
      fPrimaryDir = particle->GetMomentumDirection();
    }
  }
}

void EventAction::UpdateShowerBounds(const G4ThreeVector& pos)
{
  if (pos.z() < fShowerZMin) fShowerZMin = pos.z();
  if (pos.z() > fShowerZMax) fShowerZMax = pos.z();
}

void EventAction::AddLongitudinalDeposit(G4double depth_mm, G4double edep)
{
  if (depth_mm < 0. || depth_mm >= kProfileMaxMm || edep <= 0.) return;
  const G4int bin = static_cast<G4int>(depth_mm / kProfileMaxMm * kProfileBins);
  if (bin >= 0 && bin < kProfileBins) {
    fLongProfile[bin] += edep;
  }
}

G4double EventAction::ComputeL90() const
{
  G4double total = 0.;
  for (G4int i = 0; i < kProfileBins; ++i) {
    total += fLongProfile[i];
  }
  if (total <= 0.) return 0.;

  G4double cumulative = 0.;
  for (G4int i = 0; i < kProfileBins; ++i) {
    cumulative += fLongProfile[i];
    if (cumulative >= 0.9 * total) {
      return (static_cast<G4double>(i) + 0.5) * kProfileMaxMm / kProfileBins;
    }
  }
  return kProfileMaxMm;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  const G4int eventID = event->GetEventID();
  G4int nCherenkov = 0;
  G4int nScint = 0;

  auto* sd = dynamic_cast<SiPMSD*>(G4SDManager::GetSDMpointer()->FindSensitiveDetector("SiPMSD"));
  if (sd) {
    nCherenkov = sd->GetNCherenkov();
    nScint = sd->GetNScint();

    if (RuntimeConfig::Instance().WritePhotons()) {
      for (const auto& hit : sd->GetHits()) {
        AnalysisManager::Instance()->FillPhoton(eventID, hit.sipmID, hit.position.x(),
                                                hit.position.y(), hit.position.z(), hit.time,
                                                hit.wavelength_nm, hit.processName);
      }
    }
    for (const auto& entry : sd->GetSiPMCounts()) {
      AnalysisManager::Instance()->FillSiPMSummary(eventID, entry.first, entry.second);
    }
  }

  const G4double showerLength =
    (fShowerZMax > fShowerZMin) ? (fShowerZMax - fShowerZMin) : 0.;
  const G4double l90 = ComputeL90();

  const G4double halfBox = 500. * mm;
  const G4bool truncated =
    (fShowerZMin < -halfBox + 5. * mm) || (fShowerZMax > halfBox - 5. * mm) ||
    (l90 > 0.95 * kProfileMaxMm);

  G4int acdNTiles = 0;
  G4double acdEdepMeV = 0.;
  G4int acdVeto = 0;

  auto* acdSD = dynamic_cast<ACDSD*>(G4SDManager::GetSDMpointer()->FindSensitiveDetector("ACDSD"));
  if (acdSD) {
    acdNTiles = acdSD->GetNTilesHit();
    acdEdepMeV = acdSD->GetTotalEdepMeV();
    acdVeto = acdSD->GetVeto() ? 1 : 0;
  }

  AnalysisManager::Instance()->FillEvent(eventID, fPrimaryEnergy, fPrimaryDir.x(),
                                         fPrimaryDir.y(), fPrimaryDir.z(), nCherenkov, nScint,
                                         fEdep, showerLength, l90, truncated ? 1 : 0, acdNTiles,
                                         acdEdepMeV, acdVeto);
}

}  // namespace B1
