#pragma once

// #include <FemusInit.hpp>
// #include <MultiLevelProblem.hpp>
// #include <MultiLevelSolution.hpp>
// #include <VTKWriter.hpp>
//


using std::vector;

#include <MEDCouplingFieldDouble.hxx>
#include <MEDCouplingUMesh.hxx>
#include <MEDLoader.hxx>


#include <MEDCouplingNormalizedUnstructuredMesh.txx>
//#include <MEDCouplingPointSet.hxx>
#include <Interpolation2D.txx>
#include <InterpKernelMatrix.hxx>
#include <MEDCouplingRemapper.hxx>
//#include <InterpolationPlanar.hxx>


using namespace MEDCoupling;

namespace femus {


  class FemusMedCoupling {
    public:
      // Constructor
      FemusMedCoupling(Solution *sol) {
        _sol = sol;
        _msh = _sol->GetMesh();
        _iproc = _msh->processor_id();
      };

      ~FemusMedCoupling() {
        for (unsigned k = 0; k < 3; k++) {
          if (_medMesh[k] != NULL) {
            _medMesh[k]->decrRef();
            _medMesh[k] = NULL;
          }
        }
        for (unsigned k = 0; k < _medField.size(); k++) {
          _medField[k]->decrRef();
        }
      }

      void Femus2MedMesh(const unsigned &meshType);
      void Femus2MedMesh(const std::string &stringMeshType = "biquadratic") { Femus2MedMesh(MapMeshType(stringMeshType)); }
      void Femus2MedNodeField(const std::vector<std::string> &fieldName, const unsigned &mshType);
      void Femus2MedNodeField(const std::vector<std::string> &fieldName, const std::string &stringMeshType = "biquadratic") {
        Femus2MedNodeField(fieldName, MapMeshType(stringMeshType));
      }
      void Med2FemusNodeField(const std::vector<std::string> &fieldName, const unsigned &mshType);
      void Med2FemusNodeField(const std::vector<std::string> &fieldName, const std::string &stringMeshType = "biquadratic"){
        Med2FemusNodeField(fieldName, MapMeshType(stringMeshType));
      }

      void Femus2MedCellField(const std::vector<std::string> &fieldName, const unsigned &mshType);
      void Femus2MedCellField(const std::vector<std::string> &fieldName, const std::string &stringMeshType = "biquadratic"){
        Femus2MedCellField(fieldName, MapMeshType(stringMeshType));
      }


      unsigned MapMeshType(const std::string stringMeshType) {
        if (stringMeshType == "linear") return 0;
        else if (stringMeshType == "quadratic") return 1;
        else if (stringMeshType == "biquadratic") return 2;
        else {
          std::cout << "Error! " << stringMeshType << " is not a valid option\n";
          std::cout << "Only linear \"0\", quadratic \"1\", or biquadratic \"2\" meshTypes are available\n";
          abort();
        }
      }

      void BuildProjectionMatrixBetweenMeshes(const MEDCouplingUMesh *med_source_mesh, MEDCouplingUMesh *med_target_mesh);
      MEDCouplingFieldDouble * GetFieldProjectionOnMesh(MEDCouplingFieldDouble *sourceField, MEDCouplingUMesh *med_target_mesh, const std::string &projName);

      unsigned GetSubFieldPosition(const std::string &fieldName, const std::string &subFieldName);
      void GetSubFieldNames(const std::string &fieldName, std::vector<std::string > &subFieldName);


      void PushBackField(MEDCoupling::MEDCouplingFieldDouble *medField) {
        _medField.push_back(medField);
      }
      unsigned GetProcessId(){return _iproc;}

      Solution *_sol = NULL;
      Mesh *_msh = NULL;
      unsigned _iproc;

      MEDCoupling::MEDCouplingUMesh* _medMesh[3] = {NULL, NULL, NULL};
      std::vector <MEDCoupling::MEDCouplingFieldDouble *> _medField;
  };






} // namespace femus

