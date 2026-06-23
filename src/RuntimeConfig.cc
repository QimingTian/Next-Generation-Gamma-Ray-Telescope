#include "RuntimeConfig.hh"

#include <cstdlib>
#include <cstring>

namespace B1
{

namespace
{
G4bool EnvTruthy(const char* value, G4bool defaultValue)
{
  if (!value) return defaultValue;
  if (std::strcmp(value, "0") == 0 || std::strcmp(value, "false") == 0 ||
      std::strcmp(value, "FALSE") == 0 || std::strcmp(value, "off") == 0 ||
      std::strcmp(value, "OFF") == 0) {
    return false;
  }
  return true;
}
}  // namespace

RuntimeConfig& RuntimeConfig::Instance()
{
  static RuntimeConfig instance;
  return instance;
}

RuntimeConfig::RuntimeConfig()
{
  if (const char* v = std::getenv("G4WRITE_PHOTONS")) {
    fWritePhotons = EnvTruthy(v, false);
  }
  if (const char* v = std::getenv("G4FILTER")) {
    fApplyFilter = EnvTruthy(v, true);
  }
  if (const char* path = std::getenv("G4FILTER_CSV")) {
    fFilterCsvPath = path;
  } else {
    fFilterCsvPath = "../analysis/output/filter_transmission.csv";
  }
}

}  // namespace B1
