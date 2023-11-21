#include "fem_ttu.hpp"
#include "FemusMedCoupling.hpp"

using namespace femus;

int main(int argc, char **args) {

  unsigned numberOfUniformRefinements = 1;
  MyFEMuSMED myfemus(argc, args, "../MESH/square_tri.neu", numberOfUniformRefinements);

  FemusMedCoupling PM(myfemus._mlSol.GetSolutionLevel(numberOfUniformRefinements - 1));


  unsigned mshType = 2;
  PM.Femus2MedMesh("linear");
  PM.Femus2MedNodeField({"P", "U", "V"}, "quadratic");
  //PM.Femus2MedNodeField({"U", "P", "V"}, "biquadratic");
  //PM.Femus2MedNodeField({"P", "P", "P"}, "biquadratic");

  PM.BuildProjectionMatrixBetweenFields(PM._medField[0], static_cast <MEDCouplingUMesh*>(PM._medField[0]->getMesh()));
  //PM.BuildProjectionMatrixBetweenMeshes(PM._medMesh[0],PM._medMesh[0]);


  //unsigned nIterations = 2;

  // for (unsigned it = 1; it <= nIterations; it++) {
  //   std::cout << "ITERATION " << it << "\n";
  //   PM.solve(it);
  //   PM.write(it);
  // }
  return 0;
}
