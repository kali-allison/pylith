// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "ElasticityImplicitCUDA.hh" // implementation of class methods

#include <cuda_runtime.h>
extern cudaError_t setupKernel(const int dim, const int numBasisFuncs, const int numCells, float *K, float **geometry, float **elemMat, float **analytic_gpu, float **geometry_gpu, float **elemMat_gpu);
extern cudaError_t launchKernel(const int spaceDim, const int numBasis, const int elementBatchSize, const int numConcurrentElements,
                                const int N, float *analytic_gpu, float *geometry_gpu, float *elemMat_gpu);
extern cudaError_t cleanupKernel(float *geometry, float *elemMat, float *analytic_gpu, float *geometry_gpu, float *elemMat_gpu);

#include "Quadrature.hh" // USES Quadrature
#include "CellGeometry.hh" // USES CellGeometry

#include "pylith/materials/ElasticMaterial.hh" // USES ElasticMaterial
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/Jacobian.hh" // USES Jacobian

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/array.hh" // USES double_array
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/utils/lapack.h" // USES LAPACKdgesvd

#include "petscmat.h" // USES PetscMat
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimendional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/petscerror.h" // USES CHECK_PETSC_ERROR
#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;

typedef pylith::topology::Field<pylith::topology::Mesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::Mesh>::UpdateAddVisitor UpdateAddVisitor;
typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveMesh::order_type,PetscInt> IndicesVisitor;

// ----------------------------------------------------------------------
// Constructor
pylith::feassemble::ElasticityImplicitCUDA::ElasticityImplicitCUDA(void) :
  _dtm1(-1.0), geometry(NULL), elemMat(NULL), analytic_gpu(NULL), geometry_gpu(NULL), elemMat_gpu(NULL)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::feassemble::ElasticityImplicitCUDA::~ElasticityImplicitCUDA(void)
{ // destructor
  deallocate();
} // destructor
  
// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::ElasticityImplicitCUDA::deallocate(void)
{ // deallocate
  IntegratorElasticity::deallocate();
} // deallocate
  
// ----------------------------------------------------------------------
// Set time step for advancing from time t to time t+dt.
void
pylith::feassemble::ElasticityImplicitCUDA::timeStep(const double dt)
{ // timeStep
  if (_dt != -1.0)
    _dtm1 = _dt;
  else
    _dtm1 = dt;
  _dt = dt;
  if (0 != _material)
    _material->timeStep(_dt);
} // timeStep

// ----------------------------------------------------------------------
// Get stable time step for advancing from time t to time t+dt.
double
pylith::feassemble::ElasticityImplicitCUDA::stableTimeStep(const topology::Mesh& mesh) const
{ // stableTimeStep
  assert(0 != _material);
  return _material->stableTimeStepImplicit(mesh);
} // stableTimeStep

// ----------------------------------------------------------------------
// Set flag for setting constraints for total field solution or
// incremental field solution.
void
pylith::feassemble::ElasticityImplicitCUDA::useSolnIncr(const bool flag)
{ // useSolnIncr
  assert(0 != _material);
  _useSolnIncr = flag;
  _material->useElasticBehavior(!_useSolnIncr);
} // useSolnIncr

// ----------------------------------------------------------------------
// Integrate constributions to residual term (r) for operator.
void
pylith::feassemble::ElasticityImplicitCUDA::integrateResidual(
			  const topology::Field<topology::Mesh>& residual,
			  const double t,
			  topology::SolutionFields* const fields)
{ // integrateResidual
  /// Member prototype for _elasticityResidualXD()
  typedef void (pylith::feassemble::ElasticityImplicitCUDA::*elasticityResidual_fn_type)
    (const double_array&);
  
  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != _logger);
  assert(0 != fields);

  const int setupEvent = _logger->eventId("ElIR setup");
  const int geometryEvent = _logger->eventId("ElIR geometry");
  const int computeEvent = _logger->eventId("ElIR compute");
  const int restrictEvent = _logger->eventId("ElIR restrict");
  const int stateVarsEvent = _logger->eventId("ElIR stateVars");
  const int stressEvent = _logger->eventId("ElIR stress");
  const int updateEvent = _logger->eventId("ElIR update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
  if (cellDim != spaceDim)
    throw std::logic_error("Integration for cells with spatial dimensions "
			   "different than the spatial dimension of the "
			   "domain not implemented yet.");

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  elasticityResidual_fn_type elasticityResidualFn;
  if (1 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicitCUDA::_elasticityResidual1D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicitCUDA::_elasticityResidual2D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    elasticityResidualFn = 
      &pylith::feassemble::ElasticityImplicitCUDA::_elasticityResidual3D;
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vectors for cell values.
  double_array dispTpdtCell(numBasis*spaceDim);
  double_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;
  double_array gravVec(spaceDim);
  double_array quadPtsGlobal(numQuadPts*spaceDim);

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  double_array dispTCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispTSection = 
    fields->get("disp(t)").section();
  assert(!dispTSection.isNull());
  RestrictVisitor dispTVisitor(*dispTSection, dispTCell.size(), &dispTCell[0]);

  double_array dispTIncrCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection = 
    fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());
  RestrictVisitor dispTIncrVisitor(*dispTIncrSection,
				   dispTIncrCell.size(), &dispTIncrCell[0]);

  const ALE::Obj<RealSection>& residualSection = residual.section();
  UpdateAddVisitor residualVisitor(*residualSection, &_cellVector[0]);

  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  assert(0 != _normalizer);
  const double lengthScale = _normalizer->lengthScale();
  const double gravityScale = 
    _normalizer->pressureScale() / (_normalizer->lengthScale() *
				    _normalizer->densityScale());

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(computeEvent);

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry information for current cell
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(*c_iter);
#else
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);
#endif

    // Get state variables for cell.
    _material->retrievePropsAndVars(*c_iter);

    // Reset element vector to zero
    _resetCellVector();

    // Restrict input fields to cell
    dispTVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTVisitor);
    dispTIncrVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTIncrVisitor);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();
    const double_array& quadPtsNondim = _quadrature->quadPts();

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    dispTpdtCell = dispTCell + dispTIncrCell;

    // Compute body force vector if gravity is being used.
    if (0 != _gravityField) {
      const spatialdata::geocoords::CoordSys* cs = fields->mesh().coordsys();
      assert(0 != cs);

      // Get density at quadrature points for this cell
      const double_array& density = _material->calcDensity();

      quadPtsGlobal = quadPtsNondim;
      _normalizer->dimensionalize(&quadPtsGlobal[0], quadPtsGlobal.size(),
          lengthScale);

      // Compute action for element body forces
      for (int iQuad = 0; iQuad < numQuadPts; ++iQuad) {
        const int err = _gravityField->query(&gravVec[0], gravVec.size(),
            &quadPtsGlobal[0], spaceDim, cs);
        if (err)
          throw std::runtime_error("Unable to get gravity vector for point.");
        _normalizer->nondimensionalize(&gravVec[0], gravVec.size(),
            gravityScale);
        const double wt = quadWts[iQuad] * jacobianDet[iQuad] * density[iQuad];
        for (int iBasis = 0, iQ = iQuad * numBasis; iBasis < numBasis; ++iBasis) {
          const double valI = wt * basis[iQ + iBasis];
          for (int iDim = 0; iDim < spaceDim; ++iDim) {
            _cellVector[iBasis * spaceDim + iDim] += valI * gravVec[iDim];
          } // for
        } // for
      } // for
      PetscLogFlops(numQuadPts * (2 + numBasis * (1 + 2 * spaceDim)));
    } // if

    // residualSection->view("After gravity contribution");
    // Compute B(transpose) * sigma, first computing strains
    calcTotalStrainFn(&strainCell, basisDeriv, dispTpdtCell, 
		      numBasis, numQuadPts);
    const double_array& stressCell = _material->calcStress(strainCell, true);

    CALL_MEMBER_FN(*this, elasticityResidualFn)(stressCell);

#if 0 // DEBUGGING
    std::cout << "Updating residual for cell " << *c_iter << std::endl;
    for(int i = 0; i < _quadrature->spaceDim() * _quadrature->numBasis(); ++i) {
      std::cout << "  v["<<i<<"]: " << _cellVector[i] << std::endl;
    }
#endif
    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateClosure(*c_iter, residualVisitor);
  } // for

    _logger->eventEnd(computeEvent);
} // integrateResidual

void
pylith::feassemble::ElasticityImplicitCUDA::calculateGeometry(const int spaceDim, const double coords[], float geometry[])
{
  // Just use identity for now
  PetscMemzero(geometry, spaceDim*spaceDim*sizeof(float));
  for(int d = 0; d < spaceDim; ++d) {
    geometry[d*spaceDim+d] = 1.0;
  }
}

void
pylith::feassemble::ElasticityImplicitCUDA::integrateJacobian(
					topology::Jacobian* jacobian,
					const double t,
					topology::SolutionFields* fields)
{ // integrateJacobian
  /// Member prototype for _elasticityJacobianXD()
  typedef void (pylith::feassemble::ElasticityImplicitCUDA::*elasticityJacobian_fn_type)
    (const double_array&);

  assert(0 != _quadrature);
  assert(0 != _material);
  assert(0 != _logger);
  assert(0 != jacobian);
  assert(0 != fields);

  const int setupEvent = _logger->eventId("ElIJ setup");
  const int geometryEvent = _logger->eventId("ElIJ geometry");
  const int computeEvent = _logger->eventId("ElIJ compute");
  const int restrictEvent = _logger->eventId("ElIJ restrict");
  const int stateVarsEvent = _logger->eventId("ElIJ stateVars");
  const int updateEvent = _logger->eventId("ElIJ update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  const int tensorSize = _material->tensorSize();
  if (cellDim != spaceDim)
    throw std::logic_error("Don't know how to integrate elasticity " \
			   "contribution to Jacobian matrix for cells with " \
			   "different dimensions than the spatial dimension.");

  // Set variables dependent on dimension of cell
  totalStrain_fn_type calcTotalStrainFn;
  if (1 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain1D;
  } else if (2 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain2D;
  } else if (3 == cellDim) {
    calcTotalStrainFn = 
      &pylith::feassemble::IntegratorElasticity::_calcTotalStrain3D;
  } else
    assert(0);

  // Allocate vector for total strain
  double_array dispTpdtCell(numBasis*spaceDim);
  double_array strainCell(numQuadPts*tensorSize);
  strainCell = 0.0;

  // Get cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const int materialId = _material->id();
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", materialId);
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();

  // Get sections
  double_array dispTCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispTSection = 
    fields->get("disp(t)").section();
  assert(!dispTSection.isNull());
  RestrictVisitor dispTVisitor(*dispTSection, dispTCell.size(), &dispTCell[0]);

  double_array dispTIncrCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection = 
    fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());
  RestrictVisitor dispTIncrVisitor(*dispTIncrSection,
				   dispTIncrCell.size(), &dispTIncrCell[0]);

  // Get sparse matrix
  const PetscMat jacobianMat = jacobian->matrix();
  assert(0 != jacobianMat);

  // Get parameters used in integration.
  const double dt = _dt;
  assert(dt > 0);

  const ALE::Obj<SieveMesh::order_type>& globalOrder = 
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default",
					    dispTSection);
  assert(!globalOrder.isNull());

  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    sieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  _logger->eventEnd(setupEvent);
  _logger->eventBegin(computeEvent);

  cudaError_t cerr = setupKernel(spaceDim, numBasis, cells->size(), K, &geometry, &elemMat, &analytic_gpu, &geometry_gpu, &elemMat_gpu);CHECK_PETSC_ERROR_MSG(cerr, "Failed to setup CUDA kernel");

  // Loop over cells
  int c = 0;
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter, ++c) {
    // Compute geometry information for current cell
    coordsVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, coordsVisitor);
#if 1
    calculateGeometry(spaceDim, &coordinatesCell[0], &geometry[c*spaceDim*spaceDim]);
#else
    _quadrature->computeGeometry(coordinatesCell, *c_iter);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Get physical properties and state variables for cell.
    _material->retrievePropsAndVars(*c_iter);

    // Restrict input fields to cell
    dispTVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTVisitor);
    dispTIncrVisitor.clear();
    sieveMesh->restrictClosure(*c_iter, dispTIncrVisitor);

    // Compute current estimate of displacement at time t+dt using
    // solution increment.
    dispTpdtCell = dispTCell + dispTIncrCell;
      
    // Compute strains
    calcTotalStrainFn(&strainCell, basisDeriv, dispTpdtCell, 
		      numBasis, numQuadPts);
      
    // Get "elasticity" matrix at quadrature points for this cell
    const double_array& elasticConsts = 
      _material->calcDerivElastic(strainCell);
#endif
  } // for

  const int elementBatchSize      = 32;
  const int numConcurrentElements = 1;
  launchKernel(spaceDim, numBasis, elementBatchSize, numConcurrentElements, cells->size(), analytic_gpu, geometry_gpu, elemMat_gpu);

  // Update global matrix
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // We would need to request unique points here if we had an interpolated mesh
    IndicesVisitor jacobianVisitor(*dispTSection, *globalOrder,
                                   (int) pow(sieveMesh->getSieve()->getMaxConeSize(), sieveMesh->depth())*spaceDim);
    // Copy values
    for(int i = 0; i < numBasis*spaceDim*numBasis*spaceDim; ++i) {
      _cellMatrix[i] = elemMat[c*numBasis*spaceDim*numBasis*spaceDim+i];
    }
    // Assemble cell contribution into PETSc matrix.
    jacobianVisitor.clear();
    PetscErrorCode err = updateOperator(jacobianMat, *sieveMesh->getSieve(),
					jacobianVisitor, *c_iter,
					&_cellMatrix[0], ADD_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
  }

  cerr = cleanupKernel(geometry, elemMat, analytic_gpu, geometry_gpu, elemMat_gpu);CHECK_PETSC_ERROR_MSG(cerr, "Failed to cleanup CUDA kernel");

  _needNewJacobian = false;
  _material->resetNeedNewJacobian();

  _logger->eventEnd(computeEvent);
} // integrateJacobian
