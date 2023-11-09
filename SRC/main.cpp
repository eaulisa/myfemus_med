#include "fem_ttu.hpp"

using namespace femus;

int main(int argc, char **args) {

  MyFEMuSMED PM(argc, args, "../MESH/quad64.med");
  PM.init();
  PM.femus_to_med_mesh();
  PM.femus_to_med_field("U");

  unsigned nIterations = 2;

  for (unsigned it = 1; it <= nIterations; it++) {
    std::cout << "ITERATION " << it << "\n";
    PM.solve(it);
    PM.write(it);
  }
  return 0;
}
