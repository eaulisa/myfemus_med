
#include "fem_ttu.hpp"
//#include "FemusMedCoupling.hpp"

#include <FemusInit.hpp>
#include <NonLinearImplicitSystem.hpp>

namespace femus {

// =======================================================
  MyFEMuSMED::MyFEMuSMED(int argc, char *args[], const std::string &mesh_file, const unsigned &numberOfUniformLevels)
    : _mpinit(argc, args, MPI_COMM_WORLD),
      _numberOfUniformLevels(numberOfUniformLevels),
      _mlMsh(_numberOfUniformLevels, _numberOfUniformLevels, mesh_file.c_str(), "seventh", 1.0, nullptr),
      _mlSol(&_mlMsh),
      _mlProb(&_mlSol),
      _vtkIO(&_mlSol) {

    _mlMsh.PrintInfo();
    unsigned const dim = _mlMsh.GetDimension();

    _mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
    _mlSol.AddSolution("V", LAGRANGE, FIRST, 2);
    _mlSol.AddSolution("P", DISCONTINUOUS_POLYNOMIAL, ZERO, 2);

    _mlSol.Initialize("U", SetInitialCondition, &_mlProb);
    _mlSol.Initialize("V", SetInitialCondition, &_mlProb);
    _mlSol.Initialize("P", SetInitialCondition, &_mlProb);

    // attach the boundary condition function and generate boundary data
    _mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    _mlSol.GenerateBdc("All");

    std::vector<std::string> variablesToBePrinted = {"All"};

    // Print
    _vtkIO.SetDebugOutput(true);
    _vtkIO.Write("RESU", "biquadratic", variablesToBePrinted, 0);
  }



// =======================================================
  void MyFEMuSMED::init() {
  }

// =======================================================
  void MyFEMuSMED::solve(unsigned const &it) {
  }

// =======================================================
  void MyFEMuSMED::write(unsigned const &it) {
    _vtkIO.Write("RESU", "biquadratic", {"All"}, it);
  }

//==========================================================================================
  bool SetBoundaryCondition(const std::vector<double> &x, const char SolName[],
                            double &value, const int facename,
                            const double time) {
    bool dirichlet = false; // dirichlet

    value = 0;

    return dirichlet;
  }

//==========================================================================================
  double SetInitialCondition(const MultiLevelProblem *ml_prob,
                             const std::vector<double> &x, const char name[]) {

    double value = 0.;

    if (!strcmp(name, "U")) {
      value = -x[1];
    }

    if (!strcmp(name, "V")) {
      value = x[1] * x[0];
    }

     if (!strcmp(name, "P")) {
      value = x[0];
    }

    return value;
  }


} // namespace femus
