
#include "fem_ttu.hpp"

#include <FemusInit.hpp>
#include <NonLinearImplicitSystem.hpp>

namespace femus {

// =======================================================
MyFEMuSMED::MyFEMuSMED(int argc, char *args[], std::string const mesh_file)
    : _mpinit{argc, args, MPI_COMM_WORLD},
      _mlMsh{_numberOfUniformLevels, _numberOfUniformLevels, mesh_file.c_str(), "seventh", 1.0, nullptr},
      _mlSol{&_mlMsh},
      _mlProb{&_mlSol},
      _vtkIO(&_mlSol) {

  _mlMsh.PrintInfo();
  unsigned const dim = _mlMsh.GetDimension();

  _mlSol.AddSolution("U", LAGRANGE, SECOND, 2);

  std::vector<unsigned> solVIndex(dim);
  solVIndex[0] = _mlSol.GetIndex("U");
  _mlSol.Initialize("U", SetInitialCondition, &_mlProb);

  // attach the boundary condition function and generate boundary data
  _mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  _mlSol.GenerateBdc("U");

  std::vector<std::string> variablesToBePrinted = {"All"};

  // Print
  _vtkIO.SetDebugOutput(true);
  _vtkIO.Write("RESU", "biquadratic", variablesToBePrinted, 0);
}

static const std::vector<std::vector<INTERP_KERNEL::NormalizedCellType>> myFemusToMEDMapping = {
    {INTERP_KERNEL::NORM_HEXA8, INTERP_KERNEL::NORM_HEXA20, INTERP_KERNEL::NORM_HEXA27},
    {INTERP_KERNEL::NORM_TETRA4, INTERP_KERNEL::NORM_TETRA10, INTERP_KERNEL::NORM_TETRA10},
    {INTERP_KERNEL::NORM_PENTA6, INTERP_KERNEL::NORM_PENTA15, INTERP_KERNEL::NORM_PENTA15},
    {INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_QUAD8, INTERP_KERNEL::NORM_QUAD9},
    {INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_TRI7},
    {INTERP_KERNEL::NORM_SEG2, INTERP_KERNEL::NORM_SEG3, INTERP_KERNEL::NORM_SEG3},
};

// =======================================================
void MyFEMuSMED::femus_to_med_mesh()
{


    Solution* sol = _mlSol.GetSolutionLevel(_numberOfUniformLevels - 1);
    Mesh* msh = sol->GetMesh();    // pointer to the mesh (level) object
    elem* el = msh->el;  // pointer to the elem object in msh (level)

    const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
    unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

    //solution variable
    //unsigned soluIndex = _mlSol.GetIndex("U");    // get the position of "u" in the ml_sol object
    //unsigned soluType = _mlSol.GetSolutionType(soluIndex);    // get the finite element type for "u"
    //std::vector < double >  solu; // local solution

    std::vector < std::vector < double > > x(dim);    // local coordinates
    unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

    std::vector < int > xDof; // local to global pdeSys dofs

    MEDCoupling::MEDCouplingUMesh* mesh = MEDCoupling::MEDCouplingUMesh::New("mesh", dim);

    mesh->allocateCells( msh->_elementOffset[iproc + 1] - msh->_elementOffset[iproc]);

    for(unsigned iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom = msh->GetElementType(iel);
        unsigned nDofs  = msh->GetElementDofNumber(iel, xType);    // number of solution element dofs

        // resize local arrays
        xDof.resize(nDofs);
        //solu.resize(nDofs);

        for(unsigned k = 0; k < dim; k++) {
          x[k].resize(nDofs);
        }

        // local storage of global mapping
        for(unsigned i = 0; i < nDofs; i++) {
          //solution
          //unsigned solDof = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between solution node and solution dof
          //solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
          //sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof

          // coordinates
          unsigned iDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
          xDof[i] = iDof;
          for(unsigned k = 0; k < dim; k++) {
            x[k][i] = (*msh->_topology->_Sol[k])(iDof);      // global extraction and local storage for the element coordinates
          }
        }






        INTERP_KERNEL::NormalizedCellType cell = myFemusToMEDMapping[ielGeom][xType];
        mesh->insertNextCell(cell, nDofs, xDof.data());


    }

    mesh->finishInsertingCells();

    unsigned nvt = msh->_dofOffset[xType][iproc+1] - msh->_dofOffset[xType][iproc];

    MEDCoupling::DataArrayDouble* coordarr = MEDCoupling::DataArrayDouble::New();
    std::vector<double> xcoord(nvt*dim);
    for(unsigned i =msh->_dofOffset[xType][iproc]; i< msh->_dofOffset[xType][iproc+1];i++){
        for(unsigned k =0;k<dim;k++){
            xcoord[i*dim+k] = (*msh->_topology->_Sol[k])(i);
        }
    }
    coordarr->alloc(nvt, dim);
    std::copy(xcoord.data(), xcoord.data() + xcoord.size(), coordarr->getPointer());
    mesh->setCoords(coordarr);
    mesh->zipCoords();

    _mesh_med.push_back(mesh);
}

// =======================================================
void MyFEMuSMED::femus_to_med_field(std::string const fieldName)
{
    Solution* sol = _mlSol.GetSolutionLevel(_numberOfUniformLevels - 1);
    Mesh* msh = sol->GetMesh();    // pointer to the mesh (level) object
    unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

    MEDCoupling::MEDCouplingFieldDouble *field = MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES, MEDCoupling::ONE_TIME);
    field->setMesh(_mesh_med[0]);
    field->setName(fieldName);

    MEDCoupling::DataArrayDouble* fieldarr = MEDCoupling::DataArrayDouble::New();


    //solution variable
    unsigned fieldIndex = _mlSol.GetIndex(fieldName.c_str());    // get the position of "u" in the ml_sol object
    unsigned fieldType = _mlSol.GetSolutionType(fieldIndex);    // get the finite element type for "u"
    unsigned nvt = msh->_dofOffset[fieldType][iproc+1] - msh->_dofOffset[fieldType][iproc];

    std::vector<double> fielddata(nvt);
    for(unsigned i =msh->_dofOffset[fieldType][iproc]; i< msh->_dofOffset[fieldType][iproc+1];i++){

            fielddata[i] = (*sol->_Sol[fieldIndex])(i);

    }
    fieldarr->alloc(nvt);
    std::copy(fielddata.data(), fielddata.data() + fielddata.size(), fieldarr->getPointer());
    field->setArray(fieldarr);
    field->checkConsistencyLight();
    field->setTime(0.0, 0, 0);

    field->writeVTK("field");
}

// =======================================================
void MyFEMuSMED::init() {
}

// =======================================================
void MyFEMuSMED::solve(unsigned const it) {
}

// =======================================================
void MyFEMuSMED::write(unsigned const it) {
  _vtkIO.Write("RESU", "biquadratic", {"All"}, it);
}

//==========================================================================================
bool SetBoundaryCondition(const std::vector<double> &x, const char SolName[],
                          double &value, const int facename,
                          const double time) {
  bool dirichlet = false; // dirichlet

  if (!strcmp(SolName,
              "U")) { // strcmp compares two string in lexiographic sense.
    value = 0;
  }
  return dirichlet;
}

//==========================================================================================
double SetInitialCondition(const MultiLevelProblem *ml_prob,
                           const std::vector<double> &x, const char name[]) {

  double value = 0.;

  if (!strcmp(name, "U")) {
    value = -x[1];
  }

  return value;
}

} // namespace femus
