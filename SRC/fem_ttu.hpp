#pragma once

#include <FemusInit.hpp>
#include <MultiLevelProblem.hpp>
#include <MultiLevelSolution.hpp>
#include <VTKWriter.hpp>



using std::vector;

// #include <MEDCouplingFieldDouble.hxx>
// #include <MEDCouplingUMesh.hxx>
// #include <MEDLoader.hxx>
//
//
// #include <MEDCouplingNormalizedUnstructuredMesh.txx>
// //#include <MEDCouplingPointSet.hxx>
// #include <Interpolation2D.txx>
// #include <InterpKernelMatrix.hxx>
// #include <MEDCouplingRemapper.hxx>
//#include <InterpolationPlanar.hxx>


//using namespace MEDCoupling;

namespace femus {

  bool SetBoundaryCondition(const std::vector<double> &x, const char SolName[],
                            double &value, const int facename, const double time);
  double SetInitialCondition(const MultiLevelProblem *ml_prob,
                             const std::vector<double> &x, const char name[]);
  inline void AssembleBoussinesqAppoximation_AD(MultiLevelProblem &ml_prob) {};
  inline void AssembleBoussinesqAppoximation(MultiLevelProblem &ml_prob) {};

  class MyFEMuSMED {
    public:
      // Constructor
      MyFEMuSMED(int argc, char *args[], const std::string &mesh_file, const unsigned &numberOfUniformLevels);
      ~MyFEMuSMED() {};

      void init();
      void update() {}
      void assemble() {}
      void solve(const unsigned &it);
      void write(const unsigned &it);

      // private:
      FemusInit _mpinit;
      unsigned _numberOfUniformLevels;
      MultiLevelMesh _mlMsh;
      MultiLevelSolution _mlSol;
      MultiLevelProblem _mlProb;
      VTKWriter _vtkIO;

  };

} // namespace femus
