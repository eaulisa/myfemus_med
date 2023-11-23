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
  //PM.Femus2MedNodeField({"U", "P", "V"}, "biquadratic");
  //PM.Femus2MedNodeField({"P", "P", "P"}, "biquadratic");


  //PM.BuildProjectionMatrixBetweenMeshes(PM._medMesh[0],PM._medMesh[0]);

  MEDCouplingFieldDouble * targetfield = PM.GetFieldProjectionOnMesh(PM._medField[0], static_cast <MEDCouplingUMesh*>(PM._medMesh[0]), "P1P1");
  PM.PushBackField(targetfield);

  if(PM.GetProcessId() == 0){
    targetfield->writeVTK(targetfield->getName());
  }

  PM.Med2FemusNodeField({"U", "V"}, "linear");

  myfemus._vtkIO.Write("RESU", "biquadratic", {"All"}, 1);

  return 0;
}
