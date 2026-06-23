#ifndef B1FilterTransmission_h
#define B1FilterTransmission_h 1

#include "globals.hh"

namespace B1
{

/// Load filter T(lambda) from CSV or use built-in table.
class FilterTransmission
{
  public:
    static FilterTransmission& Instance();
    void Initialize(const G4String& csvPath);
    G4double Transmission(G4double wavelength_nm) const;

  private:
    FilterTransmission();
    void LoadFromCsv(const G4String& path);
    void UseBuiltinTable();

    static constexpr G4int kMaxPoints = 64;
    G4int fN = 0;
    G4double fWavelengths[kMaxPoints]{};
    G4double fTransmission[kMaxPoints]{};
};

}  // namespace B1

#endif
