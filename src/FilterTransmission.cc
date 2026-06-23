#include "FilterTransmission.hh"

#include <algorithm>
#include <fstream>
#include <sstream>

namespace B1
{

namespace
{
constexpr G4double kBuiltinWavelengths[] = {600, 500, 400, 300, 250, 200, 190, 178, 175, 170, 160};
constexpr G4double kBuiltinTransmission[] = {0.95, 0.94, 0.93, 0.91, 0.90, 0.89, 0.88, 0.45, 0.01, 0.001, 0.0001};
}  // namespace

FilterTransmission& FilterTransmission::Instance()
{
  static FilterTransmission instance;
  return instance;
}

FilterTransmission::FilterTransmission()
{
  UseBuiltinTable();
}

void FilterTransmission::Initialize(const G4String& csvPath)
{
  LoadFromCsv(csvPath);
}

void FilterTransmission::UseBuiltinTable()
{
  fN = static_cast<G4int>(sizeof(kBuiltinWavelengths) / sizeof(kBuiltinWavelengths[0]));
  for (G4int i = 0; i < fN; ++i) {
    fWavelengths[i] = kBuiltinWavelengths[i];
    fTransmission[i] = kBuiltinTransmission[i];
  }
}

void FilterTransmission::LoadFromCsv(const G4String& path)
{
  std::ifstream in(path.c_str());
  if (!in) {
    return;
  }

  G4String line;
  std::getline(in, line);  // header
  G4int n = 0;
  while (std::getline(in, line) && n < kMaxPoints) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    G4String token;
    G4double wl = 0.;
    G4double t = 0.;
    if (!std::getline(ss, token, ',')) continue;
    wl = std::stod(token);
    if (!std::getline(ss, token, ',')) continue;
    t = std::stod(token);
    fWavelengths[n] = wl;
    fTransmission[n] = t;
    ++n;
  }
  if (n > 0) {
    fN = n;
  }
}

G4double FilterTransmission::Transmission(G4double wavelength_nm) const
{
  if (fN <= 0) return 1.;
  if (wavelength_nm <= fWavelengths[0]) return fTransmission[0];
  if (wavelength_nm >= fWavelengths[fN - 1]) return fTransmission[fN - 1];

  for (G4int i = 0; i < fN - 1; ++i) {
    if (wavelength_nm >= fWavelengths[i] && wavelength_nm <= fWavelengths[i + 1]) {
      const G4double t =
        (wavelength_nm - fWavelengths[i]) / (fWavelengths[i + 1] - fWavelengths[i]);
      return fTransmission[i] + t * (fTransmission[i + 1] - fTransmission[i]);
    }
  }
  return fTransmission[fN - 1];
}

}  // namespace B1
