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

  // class FemusMedCoupling {
  //   public:
  //     // Constructor
  //     FemusMedCoupling(Solution *sol) {
  //       _sol = sol;
  //       _msh = _sol->GetMesh();
  //     };
  //
  //     ~FemusMedCoupling() {
  //       for (unsigned k = 0; k < 3; k++) {
  //         if (_medMesh[k] != NULL) {
  //           _medMesh[k]->decrRef();
  //           _medMesh[k] = NULL;
  //         }
  //       }
  //       for (unsigned k = 0; k < _medField.size(); k++) {
  //         _medField[k]->decrRef();
  //       }
  //     }
  //
  //     void Femus2MedMesh(const unsigned &meshType = 2);
  //     void Femus2MedNodeField(const std::vector<std::string> &fieldName, const unsigned &mshType = 2);
  //
  //     void Femus2MedMesh(const std::string &meshType = "biquadratic") {
  //       if (meshType == "linear") {
  //         Femus2MedMesh(0);
  //       }
  //       else if (meshType == "quadratic") {
  //         Femus2MedMesh(1);
  //       }
  //       else if (meshType == "biquadratic") {
  //         Femus2MedMesh(2);
  //       }
  //       else {
  //         std::cout << "Error in void MyFEMuSMED::femus_to_med_mesh (const std::string &meshType)!!!\n";
  //         std::cout << "Only linear \"0\", quadratic \"1\", or biquadratic \"2\" meshTypes are available\n";
  //         abort();
  //       }
  //     }
  //
  //     void Femus2MedNodeField(const std::vector<std::string> &fieldName, const std::string &meshType = "biquadratic") {
  //       if (meshType == "linear") {
  //         Femus2MedNodeField(fieldName, 0);
  //       }
  //       else if (meshType == "quadratic") {
  //         Femus2MedNodeField(fieldName, 1);
  //       }
  //       else if (meshType == "biquadratic") {
  //         Femus2MedNodeField(fieldName, 2);
  //       }
  //       else {
  //         std::cout << "Error in void MyFEMuSMED::femus_to_med_field (const std::vector<std::string> &fieldName, const std::string &meshType)!!!\n";
  //         std::cout << "Only linear \"0\", quadratic \"1\", or biquadratic \"2\" meshTypes are available\n";
  //         abort();
  //       }
  //     }
  //
  //     void BuildProjectionMatrixBetweenMeshes(const MEDCouplingUMesh *med_source_mesh, MEDCouplingUMesh *med_target_mesh);
  //     void BuildProjectionMatrixBetweenFields(MEDCouplingFieldDouble *sourceField, MEDCouplingUMesh *med_target_mesh);
  //
  //     Solution *_sol;
  //     Mesh *_msh;
  //
  //     MEDCoupling::MEDCouplingUMesh* _medMesh[3] = {NULL, NULL, NULL};
  //     std::vector <MEDCoupling::MEDCouplingFieldDouble *> _medField;
  //     std::vector <std::vector<std::string>> _medFieldNames;
  // };



} // namespace femus
