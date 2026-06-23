#include "OpticalMaterials.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

namespace B1
{

namespace
{
G4double WavelengthToEnergy(G4double wavelength_nm)
{
  return (1240.0 / wavelength_nm) * eV;
}

// Segel et al. (2002) / Aprile et al. (2012) VUV dispersion for LXe.
G4double LXeRefractiveIndex(G4double wavelength_nm)
{
  const G4double lambda_um = wavelength_nm * 1e-3;
  if (lambda_um < 0.14) return 2.10;
  if (lambda_um > 0.60) return 1.40;
  // Empirical fit consistent with n~1.62 at 190 nm (paper Cherenkov angle).
  return 1.35 + 0.75 * std::exp(-(lambda_um - 0.14) / 0.04);
}

// Seidel et al. (2002) Rayleigh scattering: L ~ lambda^4
G4double LXeRayleighLength(G4double wavelength_nm)
{
  const G4double lambda_m = wavelength_nm * 1e-9;
  const G4double lambda_ref = 200e-9;
  const G4double L_ref = 1.5;  // 150 cm at 200 nm (Grace et al. 2017)
  return L_ref * std::pow(lambda_m / lambda_ref, 4.0) * m;
}

G4double LXeAbsorptionLength(G4double wavelength_nm)
{
  if (wavelength_nm >= 178) return 100.0 * m;
  if (wavelength_nm >= 170) return 20.0 * m;
  return 1.0 * m;
}
}  // namespace

void OpticalMaterials::SetLXeProperties(G4Material* lxe)
{
  // Energies must be in increasing order (Geant4 requirement).
  const G4double wavelengths_nm[] = {600, 500, 400, 300, 250, 200, 190, 178, 175, 170, 160};
  const G4int nPoints = 11;

  G4double energies[nPoints];
  G4double rindex[nPoints];
  G4double abslen[nPoints];
  G4double rayleigh[nPoints];

  for (G4int i = 0; i < nPoints; ++i) {
    energies[i] = WavelengthToEnergy(wavelengths_nm[i]);
    rindex[i] = LXeRefractiveIndex(wavelengths_nm[i]);
    abslen[i] = LXeAbsorptionLength(wavelengths_nm[i]);
    rayleigh[i] = LXeRayleighLength(wavelengths_nm[i]);
  }

  auto* mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", energies, rindex, nPoints);
  mpt->AddProperty("ABSLENGTH", energies, abslen, nPoints);
  mpt->AddProperty("RAYLEIGH", energies, rayleigh, nPoints);

  // Scintillation: 25000 photons/MeV, dual exponential (Aprile et al. 2012)
  mpt->AddConstProperty("SCINTILLATIONYIELD", 25000. / MeV);
  mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
  mpt->AddConstProperty("FASTTIMECONSTANT", 4.3 * ns, true);
  mpt->AddConstProperty("SLOWTIMECONSTANT", 22.0 * ns, true);
  // FAST/SLOW components configured by G4OpticalPhysics via G4OpticalParameters

  lxe->SetMaterialPropertiesTable(mpt);
}

G4MaterialPropertiesTable* OpticalMaterials::BuildScintillationSpectrum()
{
  return nullptr;
}

}  // namespace B1
