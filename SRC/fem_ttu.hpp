#pragma once

#include <FemusInit.hpp>
#include <MultiLevelProblem.hpp>
#include <MultiLevelSolution.hpp>
#include <VTKWriter.hpp>

#include <MEDCouplingFieldDouble.hxx>
#include <MEDCouplingUMesh.hxx>
#include <MEDLoader.hxx>

namespace femus {

bool SetBoundaryCondition(const std::vector<double> &x, const char SolName[],
                          double &value, const int facename, const double time);
double SetInitialCondition(const MultiLevelProblem *ml_prob,
                           const std::vector<double> &x, const char name[]);
inline void AssembleBoussinesqAppoximation_AD(MultiLevelProblem &ml_prob){};
inline void AssembleBoussinesqAppoximation(MultiLevelProblem &ml_prob){};

class MyFEMuSMED {
public:
  // Constructor
  MyFEMuSMED(int argc, char *args[], std::string const mesh_file);
  void femus_to_med_mesh();
  void femus_to_med_field(std::string const fieldName);
  void init();
  void update() {}
  void assemble() {}
  void solve(unsigned const it);
  void write(unsigned const it);

  // private:
  FemusInit _mpinit;
  unsigned int const _numberOfUniformLevels = 1;
  MultiLevelMesh _mlMsh;
  MultiLevelSolution _mlSol;
  MultiLevelProblem _mlProb;
  VTKWriter _vtkIO;

  std::vector<MEDCoupling::MEDCouplingUMesh*> _mesh_med;
};

} // namespace femus
