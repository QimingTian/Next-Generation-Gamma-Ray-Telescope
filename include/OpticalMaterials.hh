#ifndef B1OpticalMaterials_h
#define B1OpticalMaterials_h 1

class G4Material;
class G4MaterialPropertiesTable;

namespace B1
{

/// Literature-based optical properties for liquid xenon and detector surfaces.
/// References: Segel et al. (2002), Aprile et al. (2012), Grace et al. (2017).
class OpticalMaterials
{
  public:
    static void SetLXeProperties(G4Material* lxe);
    static G4MaterialPropertiesTable* BuildScintillationSpectrum();
};

}  // namespace B1

#endif
