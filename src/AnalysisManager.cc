#include "AnalysisManager.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

AnalysisManager* AnalysisManager::fInstance = nullptr;

AnalysisManager* AnalysisManager::Instance()
{
  if (!fInstance) {
    fInstance = new AnalysisManager();
  }
  return fInstance;
}

void AnalysisManager::DeleteInstance()
{
  delete fInstance;
  fInstance = nullptr;
}

void AnalysisManager::Book()
{
  if (fBooked) return;

  auto analysis = G4AnalysisManager::Instance();
  analysis->SetDefaultFileType("csv");
  analysis->SetVerboseLevel(1);
  analysis->SetNtupleMerging(false);

  fPhotonNtuple = analysis->CreateNtuple("photons", "Optical photon hits");
  analysis->CreateNtupleIColumn(fPhotonNtuple, "EventID");
  analysis->CreateNtupleIColumn(fPhotonNtuple, "SiPMID");
  analysis->CreateNtupleDColumn(fPhotonNtuple, "X_mm");
  analysis->CreateNtupleDColumn(fPhotonNtuple, "Y_mm");
  analysis->CreateNtupleDColumn(fPhotonNtuple, "Z_mm");
  analysis->CreateNtupleDColumn(fPhotonNtuple, "Time_ns");
  analysis->CreateNtupleDColumn(fPhotonNtuple, "Wavelength_nm");
  analysis->CreateNtupleSColumn(fPhotonNtuple, "Process");
  analysis->FinishNtuple();

  fEventNtuple = analysis->CreateNtuple("events", "Per-event summary");
  analysis->CreateNtupleIColumn(fEventNtuple, "EventID");
  analysis->CreateNtupleDColumn(fEventNtuple, "E_primary_GeV");
  analysis->CreateNtupleDColumn(fEventNtuple, "DirX");
  analysis->CreateNtupleDColumn(fEventNtuple, "DirY");
  analysis->CreateNtupleDColumn(fEventNtuple, "DirZ");
  analysis->CreateNtupleIColumn(fEventNtuple, "N_Cherenkov");
  analysis->CreateNtupleIColumn(fEventNtuple, "N_Scint");
  analysis->CreateNtupleDColumn(fEventNtuple, "Edep_MeV");
  analysis->CreateNtupleDColumn(fEventNtuple, "ShowerLength_mm");
  analysis->CreateNtupleDColumn(fEventNtuple, "L90_mm");
  analysis->CreateNtupleIColumn(fEventNtuple, "L90_truncated");
  analysis->FinishNtuple();

  fSiPMNtuple = analysis->CreateNtuple("sipm_summary", "Per-SiPM photon counts");
  analysis->CreateNtupleIColumn(fSiPMNtuple, "EventID");
  analysis->CreateNtupleIColumn(fSiPMNtuple, "SiPMID");
  analysis->CreateNtupleIColumn(fSiPMNtuple, "PhotonCount");
  analysis->FinishNtuple();

  fBooked = true;
}

void AnalysisManager::Open(const G4String& fileName)
{
  G4AnalysisManager::Instance()->OpenFile(fileName);
}

void AnalysisManager::Close()
{
  G4AnalysisManager::Instance()->Write();
  G4AnalysisManager::Instance()->CloseFile();
}

void AnalysisManager::FillPhoton(G4int eventID, G4int sipmID, G4double x, G4double y,
                                 G4double z, G4double t, G4double wavelength_nm,
                                 const G4String& process)
{
  auto analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(fPhotonNtuple, 0, eventID);
  analysis->FillNtupleIColumn(fPhotonNtuple, 1, sipmID);
  analysis->FillNtupleDColumn(fPhotonNtuple, 2, x / mm);
  analysis->FillNtupleDColumn(fPhotonNtuple, 3, y / mm);
  analysis->FillNtupleDColumn(fPhotonNtuple, 4, z / mm);
  analysis->FillNtupleDColumn(fPhotonNtuple, 5, t / ns);
  analysis->FillNtupleDColumn(fPhotonNtuple, 6, wavelength_nm);
  analysis->FillNtupleSColumn(fPhotonNtuple, 7, process);
  analysis->AddNtupleRow(fPhotonNtuple);
}

void AnalysisManager::FillEvent(G4int eventID, G4double ePrimary, G4double dirX,
                                G4double dirY, G4double dirZ, G4int nCherenkov,
                                G4int nScint, G4double edep, G4double showerLength,
                                G4double l90_mm, G4int truncated)
{
  auto analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(fEventNtuple, 0, eventID);
  analysis->FillNtupleDColumn(fEventNtuple, 1, ePrimary / GeV);
  analysis->FillNtupleDColumn(fEventNtuple, 2, dirX);
  analysis->FillNtupleDColumn(fEventNtuple, 3, dirY);
  analysis->FillNtupleDColumn(fEventNtuple, 4, dirZ);
  analysis->FillNtupleIColumn(fEventNtuple, 5, nCherenkov);
  analysis->FillNtupleIColumn(fEventNtuple, 6, nScint);
  analysis->FillNtupleDColumn(fEventNtuple, 7, edep / MeV);
  analysis->FillNtupleDColumn(fEventNtuple, 8, showerLength / mm);
  analysis->FillNtupleDColumn(fEventNtuple, 9, l90_mm / mm);
  analysis->FillNtupleIColumn(fEventNtuple, 10, truncated);
  analysis->AddNtupleRow(fEventNtuple);
}

void AnalysisManager::FillSiPMSummary(G4int eventID, G4int sipmID, G4int photonCount)
{
  if (photonCount <= 0) return;
  auto analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(fSiPMNtuple, 0, eventID);
  analysis->FillNtupleIColumn(fSiPMNtuple, 1, sipmID);
  analysis->FillNtupleIColumn(fSiPMNtuple, 2, photonCount);
  analysis->AddNtupleRow(fSiPMNtuple);
}

}  // namespace B1
