#ifndef B1RuntimeConfig_h
#define B1RuntimeConfig_h 1

#include "globals.hh"

namespace B1
{

/// Runtime flags from environment variables (read once at startup).
class RuntimeConfig
{
  public:
    static RuntimeConfig& Instance();

    G4bool WritePhotons() const { return fWritePhotons; }
    G4bool ApplyFilter() const { return fApplyFilter; }
    G4String FilterCsvPath() const { return fFilterCsvPath; }

  private:
    RuntimeConfig();

    G4bool fWritePhotons = false;
    G4bool fApplyFilter = true;
    G4String fFilterCsvPath;
};

}  // namespace B1

#endif
