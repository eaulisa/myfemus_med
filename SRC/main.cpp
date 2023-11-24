#include "fem_ttu.hpp"
#include "FemusMedCoupling.hpp"

using namespace femus;

int main(int argc, char **args) {

  unsigned numberOfUniformRefinements = 4;
  MyFEMuSMED myfemus(argc, args, "../MESH/square_quad.neu", numberOfUniformRefinements);

  FemusMedCoupling F2M(myfemus._mlSol.GetSolutionLevel(numberOfUniformRefinements - 1));


  unsigned mshType = 2;
  MEDCouplingMesh* mesh0 = F2M.Femus2MedMesh("linear");
  MEDCouplingFieldDouble* nodeFieldUV = F2M.Femus2MedNodeField({"P", "U", "V"}, "biquadratic");
  MEDCouplingFieldDouble* cellFieldP = F2M.Femus2MedCellField({"P", "U", "V"}, "biquadratic");

  MEDCouplingMesh* mesh2 = nodeFieldUV->getMesh();

  MEDCouplingFieldDouble* targetfield = F2M.GetFieldProjectionOnMesh(nodeFieldUV, mesh2, "Identity");

  F2M.AddField(targetfield);

  if (F2M.GetProcessId() == 0) {
    targetfield->writeVTK(targetfield->getName());
  }

  F2M.Med2FemusNodeField({"U", "V"}, "biquadratic");

  targetfield = F2M.GetFieldProjectionOnMesh(cellFieldP, mesh2, "P0P0");
  F2M.AddField(targetfield);

  if (F2M.GetProcessId() == 0) {
    targetfield->writeVTK(targetfield->getName());
  }

  myfemus._vtkIO.Write("RESU", "biquadratic", {"All"}, 1);

  F2M.EraseAllFields();

  return 0;
}
