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

      MEDCouplingUMesh* Femus2MedMesh(const unsigned &meshType);
      MEDCouplingUMesh* Femus2MedMesh(const std::string &stringMeshType = "biquadratic") {
        return Femus2MedMesh(MapMeshType(stringMeshType));
      }
      MEDCouplingFieldDouble* Femus2MedNodeField(const std::vector<std::string> &fieldName, const unsigned &mshType);
      MEDCouplingFieldDouble* Femus2MedNodeField(const std::vector<std::string> &fieldName, const std::string &stringMeshType = "biquadratic") {
        return Femus2MedNodeField(fieldName, MapMeshType(stringMeshType));
      }


      void Med2FemusNodeField(const std::vector<std::string> &subFieldName, const std::string &stringMeshType = "biquadratic") {
        Med2FemusNodeField(subFieldName, MapMeshType(stringMeshType));
      }
      void Med2FemusNodeField(const std::vector<std::string> &subFieldName, const unsigned &mshType) {
        MEDCouplingFieldDouble* thisField = GetField(subFieldName, mshType);
        Med2FemusNodeField(thisField);
      }
      void Med2FemusNodeField(const std::string &fieldName) {
        Med2FemusNodeField(GetField(fieldName));
      }
      void Med2FemusNodeField(MEDCouplingFieldDouble* thisField);



      void Med2FemusCellField(const std::vector<std::string> &subFieldName, const std::string &stringMeshType = "biquadratic") {
        Med2FemusCellField(subFieldName, MapMeshType(stringMeshType));
      }
      void Med2FemusCellField(const std::vector<std::string> &subFieldName, const unsigned &mshType) {
        MEDCouplingFieldDouble* thisField = GetField(subFieldName, mshType);
        Med2FemusCellField(thisField);
      }
      void Med2FemusCellField(const std::string &fieldName) {
        Med2FemusCellField(GetField(fieldName));
      }
      void Med2FemusCellField(MEDCouplingFieldDouble* thisField);



      MEDCouplingFieldDouble* Femus2MedCellField(const std::vector<std::string> &fieldName, const unsigned &mshType);
      MEDCouplingFieldDouble* Femus2MedCellField(const std::vector<std::string> &fieldName, const std::string &stringMeshType = "biquadratic") {
        return Femus2MedCellField(fieldName, MapMeshType(stringMeshType));
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

      void BuildProjectionMatrixBetweenMeshes(const MEDCouplingMesh *med_source_mesh, MEDCouplingMesh *med_target_mesh);
      MEDCouplingFieldDouble * GetFieldProjectionOnMesh(MEDCouplingFieldDouble *sourceField, MEDCouplingMesh *med_target_mesh, const std::string &projName);

      unsigned GetSubFieldPosition(const std::string &fieldName, const std::string &subFieldName);
      void GetSubFieldNames(const std::string &fieldName, std::vector<std::string > &subFieldName);


      void AddField(MEDCouplingFieldDouble *field) {
        std::string fieldName = field->getName();
        for (unsigned i = 0; i < _medField.size(); i++) {
          if (_medField[i]->getName() == fieldName) {
            std::cout << "Warning a field named " << fieldName << " already exists and it is being replaced" << std::endl;
            _medField[i]->decrRef();
            _medField[i] = field;
            return;
          }
        }
        _medField.push_back(field);
      }

      unsigned GetNumberOfFields() {
        return _medField.size();
      }



      void EraseAllFields() {
        for (unsigned i = 0; i < _medField.size(); i++) {
          _medField[i]->decrRef();
        }
        _medField.resize(0);
      }

      void EraseField(const std::string &fieldName) {
        for (unsigned i = 0; i < _medField.size(); i++) {
          if (_medField[i]->getName() == fieldName) {
            _medField[i]->decrRef();
            _medField.erase(_medField.begin() + i);
            return;
          }
        }
        std::cout << "Warning a field named " << fieldName << " cannot be removed, because it does not exist" << std::endl;
      }

      void EraseField(const unsigned &i) {
        if (i < _medField.size()) {
          _medField[i]->decrRef();
          _medField.erase(_medField.begin() + i);
        }
        else {
          std::cout << "Warning there is no field in the " << i << " position to be erased " << std::endl;
        }
        return;
      }


      unsigned GetProcessId() {
        return _iproc;
      }

      std::string GetFieldName(const std::vector<std::string> &subFieldName, const unsigned &mshType);
      MEDCouplingFieldDouble* GetField(const std::vector<std::string> &subFieldName, const unsigned &mshType);
      MEDCouplingFieldDouble* GetField(const std::string & fieldName);
      MEDCouplingFieldDouble* GetField(const unsigned &i) {
        return (i < _medField.size()) ? _medField[i] : NULL;
      }


  private:
      Solution *_sol = NULL;
      Mesh *_msh = NULL;
      unsigned _iproc;

      MEDCouplingUMesh* _medMesh[3] = {NULL, NULL, NULL};
      std::vector <MEDCouplingFieldDouble *> _medField;
  };






} // namespace femus


