#include "fem_ttu.hpp"
#include "FemusMedCoupling.hpp"

using namespace femus;

int main(int argc, char **args) {

  unsigned numberOfUniformRefinements = 1;
  MyFEMuSMED myfemus(argc, args, "../MESH/square_quad.neu", numberOfUniformRefinements);

  FemusMedCoupling PM(myfemus._mlSol.GetSolutionLevel(numberOfUniformRefinements - 1));


  unsigned mshType = 2;
  PM.Femus2MedMesh("linear");
  PM.Femus2MedNodeField({"P", "U", "V"}, "biquadratic");
  PM.Femus2MedCellField({"P", "U", "V"}, "biquadratic");
  //PM.Femus2MedNodeField({"U", "P", "V"}, "biquadratic");
  //PM.Femus2MedNodeField({"P", "P", "P"}, "biquadratic");


  //PM.BuildProjectionMatrixBetweenMeshes(PM._medMesh[0],PM._medMesh[0]);

  MEDCouplingFieldDouble * targetfield = PM.GetFieldProjectionOnMesh(PM._medField[0], static_cast <MEDCouplingUMesh*>(PM._medMesh[2]), "P1P1");

  PM._medField[0]->decrRef();
  PM._medField[0] = targetfield;
  //PM.PushBackField(targetfield);

  if(PM.GetProcessId() == 0){
    targetfield->writeVTK(targetfield->getName());
  }

  PM.Med2FemusNodeField({"U", "V"}, "biquadratic");

  myfemus._vtkIO.Write("RESU", "biquadratic", {"All"}, 1);


  targetfield = PM.GetFieldProjectionOnMesh(PM._medField[1], static_cast <MEDCouplingUMesh*>(PM._medMesh[0]), "P0P0");
  if(PM.GetProcessId() == 0){
    targetfield->writeVTK(targetfield->getName());
  }
  targetfield->decrRef();


  return 0;
}
