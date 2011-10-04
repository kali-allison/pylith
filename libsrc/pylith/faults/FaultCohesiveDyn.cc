// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/SubMesh.hh" // USES SubMesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/friction/FrictionModel.hh" // USES FrictionModel
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN
#include "pylith/problems/SolverLinear.hh" // USES SolverLinear

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// Precomputing geometry significantly increases storage but gives a
// slight speed improvement.
//#define PRECOMPUTE_GEOMETRY

// ----------------------------------------------------------------------
typedef pylith::topology::Mesh::SieveMesh SieveMesh;
typedef pylith::topology::Mesh::RealSection RealSection;
typedef pylith::topology::SubMesh::SieveMesh SieveSubMesh;

typedef pylith::topology::Field<pylith::topology::SubMesh>::RestrictVisitor RestrictVisitor;
typedef pylith::topology::Field<pylith::topology::SubMesh>::UpdateAddVisitor UpdateAddVisitor;
typedef ALE::ISieveVisitor::IndicesVisitor<RealSection,SieveSubMesh::order_type,PetscInt> IndicesVisitor;

// ----------------------------------------------------------------------
const double pylith::faults::FaultCohesiveDyn::_zeroTolerance = 1.0e-12;

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void) :
  _dbInitialTract(0),
  _friction(0),
  _jacobian(0),
  _ksp(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveDyn::~FaultCohesiveDyn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FaultCohesiveDyn::deallocate(void)
{ // deallocate
  FaultCohesiveLagrange::deallocate();

  _dbInitialTract = 0; // :TODO: Use shared pointer
  _friction = 0; // :TODO: Use shared pointer

  delete _jacobian; _jacobian = 0;
  if (0 != _ksp) {
    PetscErrorCode err = KSPDestroy(&_ksp); _ksp = 0;
    CHECK_PETSC_ERROR(err);
  } // if
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for the inital tractions
void
pylith::faults::FaultCohesiveDyn::dbInitialTract(spatialdata::spatialdb::SpatialDB* db)
{ // dbInitial
  _dbInitialTract = db;
} // dbInitial

// ----------------------------------------------------------------------
// Get the friction (constitutive) model.  
void
pylith::faults::FaultCohesiveDyn::frictionModel(friction::FrictionModel* const model)
{ // frictionModel
  _friction = model;
} // frictionModel

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDyn::initialize(const topology::Mesh& mesh,
					      const double upDir[3])
{ // initialize
  assert(0 != upDir);
  assert(0 != _quadrature);
  assert(0 != _normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  // Get initial tractions using a spatial database.
  _setupInitialTractions();

  // Setup fault constitutive model.
  assert(0 != _friction);
  assert(0 != _faultMesh);
  assert(0 != _fields);
  _friction->normalizer(*_normalizer);
  _friction->initialize(*_faultMesh, _quadrature);

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(0 != cs);

  // Create field for relative velocity associated with Lagrange vertex k
  _fields->add("relative velocity", "relative_velocity");
  topology::Field<topology::SubMesh>& velRel = 
    _fields->get("relative velocity");
  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  velRel.cloneSection(dispRel);
  velRel.vectorFieldType(topology::FieldBase::VECTOR);

  //logger.stagePop();
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::faults::FaultCohesiveDyn::integrateResidual(
			     const topology::Field<topology::Mesh>& residual,
			     const double t,
			     topology::SolutionFields* const fields)
{ // integrateResidual
  assert(0 != fields);
  assert(0 != _fields);
  assert(0 != _logger);

  // Initial fault tractions have been assembled, so they do not need
  // assembling across processors.

  FaultCohesiveLagrange::integrateResidual(residual, t, fields);

  // No contribution if no initial tractions are specified.
  if (0 == _dbInitialTract)
    return;

  const int setupEvent = _logger->eventId("FaIR setup");
  const int geometryEvent = _logger->eventId("FaIR geometry");
  const int computeEvent = _logger->eventId("FaIR compute");
  const int restrictEvent = _logger->eventId("FaIR restrict");
  const int updateEvent = _logger->eventId("FaIR update");

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int cellDim = _quadrature->cellDim();
  assert(cellDim == spaceDim-1);

  // Get cohesive cell information
  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cells = 
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cells.isNull());
  const SieveMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveMesh::label_sequence::iterator cellsEnd = cells->end();
  const int numCells = cells->size();

  // Get sections associated with cohesive cells
  double_array residualCell(3*numBasis*spaceDim);
  const ALE::Obj<RealSection>& residualSection = residual.section();
  assert(!residualSection.isNull());
  UpdateAddVisitor residualVisitor(*residualSection, &residualCell[0]);

  // Get fault cell information
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& faultCells =
    faultSieveMesh->heightStratum(0);
  assert(!faultCells.isNull());
  assert(faultCells->size() == cells->size());

  // Get sections associated with fault cells
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  double_array dispRelCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dispRelSection = 
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());
  RestrictVisitor dispRelVisitor(*dispRelSection, 
				 dispRelCell.size(), &dispRelCell[0]);

  double_array initialTractionsCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& initialTractionsSection = 
    _fields->get("initial tractions").section();
  assert(!initialTractionsSection.isNull());
  RestrictVisitor initialTractionsVisitor(*initialTractionsSection, 
					  initialTractionsCell.size(),
					  &initialTractionsCell[0]);

  double_array orientationCell(numBasis*spaceDim*spaceDim);
  const ALE::Obj<RealSection>& orientationSection = 
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());
  RestrictVisitor orientationVisitor(*orientationSection, 
				     orientationCell.size(),
				     &orientationCell[0]);

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over cells
  for (SieveMesh::label_sequence::iterator c_iter=cellsBegin;
	 c_iter != cellsEnd;
       ++c_iter) {
    topology::SubMesh::SieveMesh::point_type c_fault = 
      _cohesiveToFault[*c_iter];

    // Compute geometry information for current cell
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(geometryEvent);
#endif
#if defined(PRECOMPUTE_GEOMETRY)
    _quadrature->retrieveGeometry(c_fault);
#else
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, c_fault);
#endif

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(geometryEvent);
    _logger->eventBegin(restrictEvent);
#endif

    // Restrict input fields to cell
    initialTractionsVisitor.clear();
    sieveMesh->restrictClosure(c_fault, initialTractionsVisitor);

    dispRelVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, dispRelVisitor);

    orientationVisitor.clear();
    faultSieveMesh->restrictClosure(c_fault, orientationVisitor);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& basisDeriv = _quadrature->basisDeriv();
    const double_array& jacobianDet = _quadrature->jacobianDet();

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    residualCell = 0.0;

    // Compute action for positive side of fault and Lagrange constraints
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];
      const int iQ = iQuad * numBasis;
      for (int iBasis=0; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];

	// Add entries to residual at
	// iN = DOF at negative node
	// iP = DOF at positive node
	// iL = DOF at constraint node

	// Indices for negative vertex
	const int iBN = 0*numBasis*spaceDim + iBasis*spaceDim;
	
	// Indices for positive vertex
	const int iBP = 1*numBasis*spaceDim + iBasis*spaceDim;
	
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
          const double valIJ = valI * basis[iQ+jBasis];

	  // Indices for fault vertex
	  const int jB = jBasis*spaceDim;

          for (int iDim=0; iDim < spaceDim; ++iDim) {
	    // negative side of the fault
            residualCell[iBN + iDim] -= valIJ * initialTractionsCell[jB + iDim];
	    
	    // positive side of the fault
            residualCell[iBP + iDim] += valIJ * initialTractionsCell[jB + iDim];
	    
#if 0
	    std::cout << "iBasis: " << iBasis
		      << ", jBasis: " << jBasis
		      << ", iDim: " << iDim
		      << ", valIJ: " << valIJ
		      << ", jacobianDet: " << jacobianDet[iQuad]
		      << ", residualN: " << residualCell[iBN + iDim]
		      << ", residualP: " << residualCell[iBP + iDim]
		      << std::endl;
#endif

	  } // for
        } // for
      } // for
    } // for


    // Only apply initial tractions if there is no opening.
    // If there is opening, zero out initial tractions
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      const int iB = iBasis*spaceDim;
      const int iO = iBasis*spaceDim*spaceDim;
      
      double slipNormal = 0.0;
      const int indexN = spaceDim - 1;
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipNormal += 
	  orientationCell[iO + indexN*spaceDim+jDim]*dispRelCell[iB+jDim];
      } // for

      if (slipNormal > _zeroTolerance) {
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  residualCell[iB+iDim] = 0.0;
	} // for
      } // if
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Assemble cell contribution into field
    residualVisitor.clear();
    sieveMesh->updateClosure(*c_iter, residualVisitor);

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numCells*spaceDim*spaceDim*8);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
} // integrateResidual

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::faults::FaultCohesiveDyn::updateStateVars(
				      const double t,
				      topology::SolutionFields* const fields)
{ // updateStateVars
  assert(0 != fields);
  assert(0 != _fields);

  _updateVelRel(*fields);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  double_array tractionTpdtVertex(spaceDim); // Fault coordinate system

  // Get sections
  double_array slipVertex(spaceDim);
  const ALE::Obj<RealSection>& dispRelSection = 
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());

  double_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& velRelSection =
      _fields->get("relative velocity").section();
  assert(!velRelSection.isNull());

  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get relative displacement
    assert(spaceDim == dispRelSection->getFiberDimension(v_fault));
    const double* dispRelVertex = dispRelSection->restrictPoint(v_fault);
    assert(dispRelVertex);

    // Get relative velocity
    assert(spaceDim == velRelSection->getFiberDimension(v_fault));
    const double* velRelVertex = velRelSection->restrictPoint(v_fault);
    assert(velRelVertex);

    // Get orientation
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const double* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(spaceDim == dispTIncrSection->getFiberDimension(v_lagrange));
    const double* lagrangeTIncrVertex = 
      dispTIncrSection->restrictPoint(v_lagrange);

    // Compute slip, slip rate, and fault traction (Lagrange
    // multiplier) at time t+dt in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  dispRelVertex[jDim];
	slipRateVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  velRelVertex[jDim];
	tractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (lagrangeTVertex[jDim]+lagrangeTIncrVertex[jDim]);
      } // for
    } // for

    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    switch (spaceDim) { // switch
    case 1: { // case 1
      const double slipMag = 0.0;
      const double slipRateMag = 0.0;
      const double tractionNormal = tractionTpdtVertex[0];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 1
    case 2: { // case 2
      const double slipMag = fabs(slipVertex[0]);
      const double slipRateMag = fabs(slipRateVertex[0]);
      const double tractionNormal = tractionTpdtVertex[1];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 2
    case 3: { // case 3
      const double slipMag = 
	sqrt(slipVertex[0]*slipVertex[0] + slipVertex[1]*slipVertex[1]);
      const double slipRateMag = 
	sqrt(slipRateVertex[0]*slipRateVertex[0] + 
	     slipRateVertex[1]*slipRateVertex[1]);
      const double tractionNormal = tractionTpdtVertex[2];
      _friction->updateStateVars(slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in "
			     "FaultCohesiveDyn::updateStateVars().");
    } // switch
  } // for
} // updateStateVars

// ----------------------------------------------------------------------
// Constrain solution based on friction.
void
pylith::faults::FaultCohesiveDyn::constrainSolnSpace(
				    topology::SolutionFields* const fields,
				    const double t,
				    const topology::Jacobian& jacobian)
{ // constrainSolnSpace
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (double_array*,
     const double_array&,
     const double_array&,
     const double_array&);

  assert(0 != fields);
  assert(0 != _quadrature);
  assert(0 != _fields);
  assert(0 != _friction);

  _updateVelRel(*fields);
  _sensitivitySetup(jacobian);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  double_array tractionTpdtVertex(spaceDim);

  // Get sections
  double_array dDispRelVertex(spaceDim);
  const ALE::Obj<RealSection>& dispRelSection = 
    _fields->get("relative disp").section();
  assert(!dispRelSection.isNull());

  double_array dSlipVertex(spaceDim);
  double_array slipVertex(spaceDim);
  double_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& velRelSection =
      _fields->get("relative velocity").section();
  assert(!velRelSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  double_array dispTIncrVertexN(spaceDim);
  double_array dispTIncrVertexP(spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  double_array dLagrangeTpdtVertex(spaceDim);
  double_array dLagrangeTpdtVertexGlobal(spaceDim);
  const ALE::Obj<RealSection>& dLagrangeTpdtSection =
      _fields->get("sensitivity dLagrange").section();
  assert(!dLagrangeTpdtSection.isNull());

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 1:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D;
    break;
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDyn::constrainSolnSpace().");
  } // switch


  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    // Get relative displacement
    assert(spaceDim == dispRelSection->getFiberDimension(v_fault));
    const double* dispRelVertex = dispRelSection->restrictPoint(v_fault);
    assert(dispRelVertex);

    // Get relative velocity
    assert(spaceDim == velRelSection->getFiberDimension(v_fault));
    const double* velRelVertex = velRelSection->restrictPoint(v_fault);
    assert(velRelVertex);

    // Get orientation
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const double* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);

    assert(spaceDim == dispTIncrSection->getFiberDimension(v_lagrange));
    const double* lagrangeTIncrVertex = 
      dispTIncrSection->restrictPoint(v_lagrange);

    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipVertex[jDim] += orientationVertex[iDim*spaceDim+jDim] *
	  dispRelVertex[jDim];
	slipRateVertex[jDim] += orientationVertex[iDim*spaceDim+jDim] *
	  velRelVertex[jDim];
	tractionTpdtVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  (lagrangeTVertex[jDim] + lagrangeTIncrVertex[jDim]);
      } // for
    } // for

    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    dLagrangeTpdtVertex = 0.0;
    CALL_MEMBER_FN(*this,
		   constrainSolnSpaceFn)(&dLagrangeTpdtVertex,
					 slipVertex, slipRateVertex,
					 tractionTpdtVertex);

    // Rotate traction back to global coordinate system.
    dLagrangeTpdtVertexGlobal = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	dLagrangeTpdtVertexGlobal[iDim] += 
	  orientationVertex[jDim*spaceDim+iDim] * dLagrangeTpdtVertex[jDim];
      } // for
    } // for

#if 0 // debugging
    std::cout << "slipVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipVertex[iDim];
    std::cout << ",  slipRateVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipRateVertex[iDim];
    std::cout << ",  tractionVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << tractionTpdtVertex[iDim];
    std::cout << ",  dLagrangeTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dLagrangeTpdtVertex[iDim];
    std::cout << ",  dLagrangeTpdtVertexGlobal: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dLagrangeTpdtVertexGlobal[iDim];
    std::cout << std::endl;
#endif
     
    assert(dLagrangeTpdtVertexGlobal.size() ==
        dLagrangeTpdtSection->getFiberDimension(v_fault));
    dLagrangeTpdtSection->updatePoint(v_fault, &dLagrangeTpdtVertexGlobal[0]);
  } // for

  // Solve sensitivity problem for negative side of the fault.
  bool negativeSide = true;
  _sensitivityUpdateJacobian(negativeSide, jacobian, *fields);
  _sensitivityReformResidual(negativeSide);
  _sensitivitySolve();
  _sensitivityUpdateSoln(negativeSide);

  // Solve sensitivity problem for positive side of the fault.
  negativeSide = false;
  _sensitivityUpdateJacobian(negativeSide, jacobian, *fields);
  _sensitivityReformResidual(negativeSide);
  _sensitivitySolve();
  _sensitivityUpdateSoln(negativeSide);

  // Update slip field based on solution of sensitivity problem and
  // increment in Lagrange multipliers.
  const ALE::Obj<RealSection>& sensDispRelSection =
    _fields->get("sensitivity relative disp").section();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get current relative displacement for updating.
    assert(spaceDim == dispRelSection->getFiberDimension(v_fault));
    const double* dispRelVertex = dispRelSection->restrictPoint(v_fault);
    assert(0 != dispRelVertex);

    // Get change in relative displacement from sensitivity solve.
    assert(spaceDim == sensDispRelSection->getFiberDimension(v_fault));
    const double* sensDispRelVertex = 
      sensDispRelSection->restrictPoint(v_fault);
    assert(sensDispRelVertex);

    // Get orientation.
    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    // Get Lagrange multiplier at time t
    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const double* lagrangeTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(lagrangeTVertex);

    // Get Lagrange multiplier increment at time t
    assert(spaceDim == dispTIncrSection->getFiberDimension(v_lagrange));
    const double* lagrangeTIncrVertex = 
      dispTIncrSection->restrictPoint(v_lagrange);
    assert(lagrangeTIncrVertex);

    // Get change in Lagrange multiplier.
    dLagrangeTpdtSection->restrictPoint(v_fault, &dLagrangeTpdtVertex[0],
					dLagrangeTpdtVertex.size());

    // Only update slip if Lagrange multiplier is changing
    double dLagrangeMag = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim)
      dLagrangeMag += dLagrangeTpdtVertex[iDim]*dLagrangeTpdtVertex[iDim];
    if (0.0 == dLagrangeMag)
      continue; // No change, so continue

    // Compute slip and change in slip in fault coordinates.
    dSlipVertex = 0.0;
    slipVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	slipVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  dispRelVertex[jDim];
	dSlipVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] * 
	  sensDispRelVertex[jDim];
      } // for
    } // for

    // Compute normal traction in fault coordinates.
    double tractionNormal = 0.0;
    const int indexN = spaceDim - 1;
    for (int jDim=0; jDim < spaceDim; ++jDim) {
      tractionNormal += orientationVertex[indexN*spaceDim+jDim] *
	(lagrangeTVertex[jDim] + lagrangeTIncrVertex[jDim] + 
	 dLagrangeTpdtVertex[jDim]);
    } // for

    // Do not allow fault interpenetration and set fault opening to
    // zero if fault is under compression.
    if (tractionNormal < -_zeroTolerance || 
	slipVertex[indexN] + dSlipVertex[indexN] < 0.0) {
      dSlipVertex[indexN] = -slipVertex[indexN];
    } // if

    // Compute change relative displacement from change in slip.
    dDispRelVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	dDispRelVertex[iDim] += orientationVertex[jDim*spaceDim+iDim] *
	  dSlipVertex[jDim];
      } // for
    } // for

    // Set change in slip.
    assert(dDispRelVertex.size() ==
        dispRelSection->getFiberDimension(v_fault));
    dispRelSection->updateAddPoint(v_fault, &dDispRelVertex[0]);
    
    // Update Lagrange multiplier increment.
    assert(dLagrangeTpdtVertex.size() ==
	   dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updateAddPoint(v_lagrange, &dLagrangeTpdtVertex[0]);

    // Update displacement field
    dispTIncrVertexN = -0.5*dDispRelVertex;
    assert(dispTIncrVertexN.size() ==
	   dispTIncrSection->getFiberDimension(v_negative));
    dispTIncrSection->updateAddPoint(v_negative, &dispTIncrVertexN[0]);
    
    dispTIncrVertexP = -dispTIncrVertexN;
    assert(dispTIncrVertexP.size() ==
	   dispTIncrSection->getFiberDimension(v_positive));
    dispTIncrSection->updateAddPoint(v_positive, &dispTIncrVertexP[0]);
    
  } // for

#if 0 // DEBUGGING
  //dLagrangeTpdtSection->view("AFTER dLagrange");
  dispTIncrSection->view("AFTER DISP INCR (t->t+dt)");
  //dispRelSection->view("AFTER RELATIVE DISPLACEMENT");
  //velRelSection->view("AFTER RELATIVE VELOCITY");
#endif
} // constrainSolnSpace

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveDyn::adjustSolnLumped(
			 topology::SolutionFields* const fields,
			 const topology::Field<topology::Mesh>& jacobian)
{ // adjustSolnLumped
#if 0
  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (double_array*,
     const double_array&,
     const double_array&,
     const double_array&);

  /// Member prototype for _sensitivitySolveLumpedXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*sensitivitySolveLumped_fn_type)
    (double_array*,
     const double_array&,
     const double_array&,
     const double_array&);

  /// Member prototype for _adjustSolnLumpedXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*adjustSolnLumped_fn_type)
    (double_array*, double_array*, double_array*,
     const double_array&, const double_array&,
     const double_array&, const double_array&,
     const double_array&, const double_array&,
     const double_array&, const double_array&);


  assert(0 != fields);
  assert(0 != _quadrature);

  // Cohesive cells with conventional vertices i and j, and constraint
  // vertex k require three adjustments to the solution:
  //
  //   * DOF k: Compute increment in Lagrange multipliers
  //            dl_k = S^{-1} (-C_ki (A_i^{-1} r_i - C_kj A_j^{-1} r_j + u_i - u_j) - d_k)
  //            S = C_ki (A_i^{-1} + A_j^{-1}) C_ki^T
  //
  //   * Adjust Lagrange multipliers to match friction criterion
  //
  //   * DOF k: Adjust displacement increment (solution) to create slip
  //     consistent with Lagrange multiplier constraints
  //            du_i = +A_i^-1 C_ki^T dlk
  //            du_j = -A_j^-1 C_kj^T dlk

  const int setupEvent = _logger->eventId("FaAS setup");
  const int geometryEvent = _logger->eventId("FaAS geometry");
  const int computeEvent = _logger->eventId("FaAS compute");
  const int restrictEvent = _logger->eventId("FaAS restrict");
  const int updateEvent = _logger->eventId("FaAS update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();
  const int orientationSize = spaceDim * spaceDim;

  // Allocate arrays for vertex values
  double_array tractionTVertex(spaceDim);
  double_array tractionTpdtVertex(spaceDim);
  double_array slipTpdtVertex(spaceDim);
  double_array lagrangeTpdtVertex(spaceDim);
  double_array dLagrangeTpdtVertex(spaceDim);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  // Get section information
  double_array slipVertex(spaceDim);
  const ALE::Obj<RealSection>& slipSection = _fields->get("slip").section();
  assert(!slipSection.isNull());

  double_array slipRateVertex(spaceDim);
  const ALE::Obj<RealSection>& slipRateSection =
      _fields->get("slip rate").section();
  assert(!slipRateSection.isNull());

  double_array orientationVertex(orientationSize);
  const ALE::Obj<RealSection>& orientationSection =
      _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  double_array dispTVertexN(spaceDim);
  double_array dispTVertexP(spaceDim);
  double_array lagrangeTVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTSection = fields->get("disp(t)").section();
  assert(!dispTSection.isNull());

  double_array dispTIncrVertexN(spaceDim);
  double_array dispTIncrVertexP(spaceDim);
  double_array lagrangeTIncrVertex(spaceDim);
  const ALE::Obj<RealSection>& dispTIncrSection =
      fields->get("dispIncr(t->t+dt)").section();
  assert(!dispTIncrSection.isNull());

  const ALE::Obj<RealSection>& dispTIncrAdjSection = fields->get(
    "dispIncr adjust").section();
  assert(!dispTIncrAdjSection.isNull());

  double_array jacobianVertexN(spaceDim);
  double_array jacobianVertexP(spaceDim);
  const ALE::Obj<RealSection>& jacobianSection = jacobian.section();
  assert(!jacobianSection.isNull());

  double_array residualVertexN(spaceDim);
  double_array residualVertexP(spaceDim);
  const ALE::Obj<RealSection>& residualSection =
      fields->get("residual").section();

  const ALE::Obj<SieveMesh>& sieveMesh = fields->mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::order_type>& globalOrder =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", 
					    jacobianSection);
  assert(!globalOrder.isNull());

  adjustSolnLumped_fn_type adjustSolnLumpedFn;
  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  sensitivitySolveLumped_fn_type sensitivitySolveLumpedFn;
  switch (spaceDim) { // switch
  case 1:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped1D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D;
    sensitivitySolveLumpedFn =
      &pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped1D;
    break;
  case 2: 
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped2D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D;
    sensitivitySolveLumpedFn =
      &pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped2D;
    break;
  case 3:
    adjustSolnLumpedFn = 
      &pylith::faults::FaultCohesiveDyn::_adjustSolnLumped3D;
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D;
    sensitivitySolveLumpedFn =
      &pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in "
			   "FaultCohesiveDyn::adjustSolnLumped.");
  } // switch

  _logger->eventEnd(setupEvent);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get slip
    slipSection->restrictPoint(v_fault, &slipVertex[0], slipVertex.size());

    // Get slip rate
    slipRateSection->restrictPoint(v_fault, &slipRateVertex[0],
				   slipRateVertex.size());
    
    // Get fault orientation
    orientationSection->restrictPoint(v_fault, &orientationVertex[0],
				      orientationVertex.size());
    
    // Get Jacobian at vertices on positive and negative sides of the fault.
    jacobianSection->restrictPoint(v_negative, &jacobianVertexN[0],
				   jacobianVertexN.size());
    jacobianSection->restrictPoint(v_positive, &jacobianVertexP[0],
				   jacobianVertexP.size());
    
    // Get residual at cohesive cell's vertices.
    residualSection->restrictPoint(v_negative, &residualVertexN[0], 
				   residualVertexN.size());
    residualSection->restrictPoint(v_positive, &residualVertexP[0], 
				   residualVertexP.size());

    // Get disp(t) at cohesive cell's vertices.
    dispTSection->restrictPoint(v_negative, &dispTVertexN[0], 
				dispTVertexN.size());
    dispTSection->restrictPoint(v_positive, &dispTVertexP[0], 
				dispTVertexP.size());
    
    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    dispTSection->restrictPoint(v_lagrange, &lagrangeTVertex[0],
				lagrangeTVertex.size());
    dispTIncrSection->restrictPoint(v_lagrange, &lagrangeTIncrVertex[0],
				    lagrangeTIncrVertex.size());
    

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    CALL_MEMBER_FN(*this,
		   adjustSolnLumpedFn)(&lagrangeTIncrVertex,
				       &dispTIncrVertexN,
				       &dispTIncrVertexP,
				       slipVertex,
				       orientationVertex,
				       dispTVertexN,
				       dispTVertexP,
				       residualVertexN,
				       residualVertexP,
				       jacobianVertexN,
				       jacobianVertexP);

    
    // Compute Lagrange multiplier at time t+dt
    lagrangeTpdtVertex = lagrangeTVertex + lagrangeTIncrVertex;
    dLagrangeTpdtVertex = 0.0;
    
    // :TODO: Rotate fault tractions to fault coordinate system.
    
    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    CALL_MEMBER_FN(*this,
		   constrainSolnSpaceFn)(&dLagrangeTpdtVertex,
					 slipVertex, slipRateVertex,
					 tractionTpdtVertex);
    CALL_MEMBER_FN(*this,
       sensitivitySolveLumpedFn)(&slipVertex,
           dLagrangeTpdtVertex, jacobianVertexN, jacobianVertexP);

    // :TODO: Rotate fault tractions to global coordinate system.
    
    lagrangeTIncrVertex += dLagrangeTpdtVertex;

    // :TODO: Refactor this into sensitivitySolveLumpedXD().
    switch (spaceDim) { // switch
    case 1: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexP[0] > 0.0);

      const double dlp = lagrangeTIncrVertex[0];

      // Update displacements at negative vertex
      dispTIncrVertexN[0] = +1.0 / jacobianVertexN[0] * dlp;
  
      // Update displacements at positive vertex
      dispTIncrVertexP[0] = -1.0 / jacobianVertexP[0] * dlp;
  
      break;
    } // case 1
    case 2: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      assert(jacobianVertexP[1] > 0.0);

      // Check to make sure Jacobian is same at all DOF for
      // vertices i and j (means S is diagonal with equal enties).
      assert(jacobianVertexN[0] == jacobianVertexN[1]);
      assert(jacobianVertexP[0] == jacobianVertexP[1]);

      const double Cpx = orientationVertex[0];
      const double Cpy = orientationVertex[1];
      const double Cqx = orientationVertex[2];
      const double Cqy = orientationVertex[3];

      const double dlp = lagrangeTIncrVertex[0];
      const double dlq = lagrangeTIncrVertex[1];

      const double dlx = Cpx * dlp + Cqx * dlq;
      const double dly = Cpy * dlp + Cqy * dlq;
  
      // Update displacements at negative vertex.
      dispTIncrVertexN[0] = dlx / jacobianVertexN[0];
      dispTIncrVertexN[1] = dly / jacobianVertexN[0];
  
      // Update displacements at positive vertex.
      dispTIncrVertexP[0] = -dlx / jacobianVertexP[0];
      dispTIncrVertexP[1] = -dly / jacobianVertexP[0];

      break;
    } // case 2
    case 3: {
      assert(jacobianVertexN[0] > 0.0);
      assert(jacobianVertexN[1] > 0.0);
      assert(jacobianVertexN[2] > 0.0);
      assert(jacobianVertexP[0] > 0.0);
      assert(jacobianVertexP[1] > 0.0);
      assert(jacobianVertexP[2] > 0.0);

      // Check to make sure Jacobian is same at all DOF for
      // vertices i and j (means S is diagonal with equal enties).
      assert(jacobianVertexN[0] == jacobianVertexN[1] && 
	     jacobianVertexN[0] == jacobianVertexN[2]);
      assert(jacobianVertexP[0] == jacobianVertexP[1] && 
	     jacobianVertexP[0] == jacobianVertexP[2]);

      const double Cpx = orientationVertex[0];
      const double Cpy = orientationVertex[1];
      const double Cpz = orientationVertex[2];
      const double Cqx = orientationVertex[3];
      const double Cqy = orientationVertex[4];
      const double Cqz = orientationVertex[5];
      const double Crx = orientationVertex[6];
      const double Cry = orientationVertex[7];
      const double Crz = orientationVertex[8];

      const double dlp = lagrangeTIncrVertex[0];
      const double dlq = lagrangeTIncrVertex[1];
      const double dlr = lagrangeTIncrVertex[2];

      const double dlx = Cpx * dlp + Cqx * dlq + Crx * dlr;
      const double dly = Cpy * dlp + Cqy * dlq + Cry * dlr;
      const double dlz = Cpz * dlp + Cqz * dlq + Crz * dlr;

      // Update displacements at negative vertex.
      dispTIncrVertexN[0] = dlx / jacobianVertexN[0];
      dispTIncrVertexN[1] = dly / jacobianVertexN[1];
      dispTIncrVertexN[2] = dlz / jacobianVertexN[2];

      // Update displacements at positive vertex.
      dispTIncrVertexP[0] = -dlx / jacobianVertexP[0];
      dispTIncrVertexP[1] = -dly / jacobianVertexP[1];
      dispTIncrVertexP[2] = -dlz / jacobianVertexP[2];

      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in "
			     "FaultCohesiveDyn::adjustSolnLumped().");
    } // switch

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Compute contribution to adjusting solution only if Lagrange
    // constraint is local (the adjustment is assembled across processors).
    if (globalOrder->isLocal(v_lagrange)) {
      // Adjust displacements to account for Lagrange multiplier values
      // (assumed to be zero in perliminary solve).
      assert(dispTIncrVertexN.size() == 
	     dispTIncrAdjSection->getFiberDimension(v_negative));
      dispTIncrAdjSection->updateAddPoint(v_negative, &dispTIncrVertexN[0]);
      
      assert(dispTIncrVertexP.size() == 
	     dispTIncrAdjSection->getFiberDimension(v_positive));
      dispTIncrAdjSection->updateAddPoint(v_positive, &dispTIncrVertexP[0]);
    } // if

    // The Lagrange multiplier and slip are NOT assembled across processors.

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry of 1.0.
    assert(lagrangeTIncrVertex.size() == 
	   dispTIncrSection->getFiberDimension(v_lagrange));
    dispTIncrSection->updatePoint(v_lagrange, &lagrangeTIncrVertex[0]);

    // Update the slip estimate based on adjustment to the Lagrange
    // multiplier values.
    assert(slipVertex.size() ==
        slipSection->getFiberDimension(v_fault));
    slipSection->updatePoint(v_fault, &slipVertex[0]);
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif
#endif
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field<pylith::topology::SubMesh>&
pylith::faults::FaultCohesiveDyn::vertexField(const char* name,
                                               const topology::SolutionFields* fields)
{ // vertexField
  assert(0 != _faultMesh);
  assert(0 != _quadrature);
  assert(0 != _normalizer);
  assert(0 != _fields);
  assert(0 != _friction);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  double scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field<topology::SubMesh>& dispRel = 
      _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    _globalToFault(&buffer);
    return buffer;

  } else if (0 == strcasecmp("slip_rate", name)) {
    const topology::Field<topology::SubMesh>& velRel = 
      _fields->get("relative velocity");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.copy(velRel);
    buffer.label("slip_rate");
    _globalToFault(&buffer);
    return buffer;

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      0);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.label("strike_dir");
    buffer.scale(1.0);
    buffer.copy(dirSection);
    return buffer;

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      1);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.label("dip_dir");
    buffer.scale(1.0);
    buffer.copy(dirSection);
    return buffer;

  } else if (0 == strcasecmp("normal_dir", name)) {
    const ALE::Obj<RealSection>& orientationSection = _fields->get(
      "orientation").section();
    assert(!orientationSection.isNull());
    const int space = (0 == cohesiveDim) ? 0 : (1 == cohesiveDim) ? 1 : 2;
    const ALE::Obj<RealSection>& dirSection = orientationSection->getFibration(
      space);
    assert(!dirSection.isNull());
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    buffer.label("normal_dir");
    buffer.scale(1.0);
    buffer.copy(dirSection);
    return buffer;

  } else if (0 == strcasecmp("initial_traction", name)) {
    assert(0 != _dbInitialTract);
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    topology::Field<topology::SubMesh>& tractions =
        _fields->get("initial tractions");
    buffer.copy(tractions);
    _globalToFault(&buffer);
    return buffer;

  } else if (0 == strcasecmp("traction", name)) {
    assert(0 != fields);
    const topology::Field<topology::Mesh>& dispT = fields->get("disp(t)");
    _allocateBufferVectorField();
    topology::Field<topology::SubMesh>& buffer =
        _fields->get("buffer (vector)");
    _calcTractions(&buffer, dispT);
    return buffer;

  } else if (_friction->hasPropStateVar(name)) {
    return _friction->getField(name);

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveDyn::vertexField().");

  // Satisfy return values
  assert(0 != _fields);
  const topology::Field<topology::SubMesh>& buffer = _fields->get(
    "buffer (vector)");

  return buffer;
} // vertexField

// ----------------------------------------------------------------------
void
pylith::faults::FaultCohesiveDyn::_setupInitialTractions(void)
{ // _setupInitialTractions
  assert(0 != _normalizer);

  // If no initial tractions specified, leave method
  if (0 == _dbInitialTract)
    return;

  assert(0 != _normalizer);
  const double pressureScale = _normalizer->pressureScale();
  const double lengthScale = _normalizer->lengthScale();

  const int spaceDim = _quadrature->spaceDim();

  // Create section to hold initial tractions.
  _fields->add("initial tractions", "initial_tractions");
  topology::Field<topology::SubMesh>& initialTractions = 
    _fields->get("initial tractions");
  topology::Field<topology::SubMesh>& dispRel = _fields->get("relative disp");
  initialTractions.cloneSection(dispRel);
  initialTractions.scale(pressureScale);

  double_array initialTractionsVertex(spaceDim);
  double_array initialTractionsVertexGlobal(spaceDim);
  const ALE::Obj<RealSection>& initialTractionsSection = 
    initialTractions.section();
  assert(!initialTractionsSection.isNull());

  const ALE::Obj<RealSection>& orientationSection =
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  const spatialdata::geocoords::CoordSys* cs = _faultMesh->coordsys();
  assert(0 != cs);

  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  double_array coordsVertex(spaceDim);
  const ALE::Obj<RealSection>& coordsSection =
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordsSection.isNull());


  assert(0 != _dbInitialTract);
  _dbInitialTract->open();
  switch (spaceDim) { // switch
  case 1: {
    const char* valueNames[] = { "traction-normal" };
    _dbInitialTract->queryVals(valueNames, 1);
    break;
  } // case 1
  case 2: {
    const char* valueNames[] = { "traction-shear", "traction-normal" };
    _dbInitialTract->queryVals(valueNames, 2);
    break;
  } // case 2
  case 3: {
    const char* valueNames[] = { "traction-shear-leftlateral",
				 "traction-shear-updip", "traction-normal" };
    _dbInitialTract->queryVals(valueNames, 3);
    break;
  } // case 3
  default:
    std::cerr << "Bad spatial dimension '" << spaceDim << "'." << std::endl;
    assert(0);
    throw std::logic_error("Bad spatial dimension in Neumann.");
  } // switch

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    coordsSection->restrictPoint(v_fault, &coordsVertex[0], coordsVertex.size());

    assert(spaceDim*spaceDim == orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    _normalizer->dimensionalize(&coordsVertex[0], coordsVertex.size(),
				lengthScale);

    initialTractionsVertex = 0.0;
    int err = _dbInitialTract->query(&initialTractionsVertex[0], 
				     initialTractionsVertex.size(),
				     &coordsVertex[0], coordsVertex.size(), cs);
    if (err) {
      std::ostringstream msg;
      msg << "Could not find parameters for physical properties at \n" << "(";
      for (int i = 0; i < spaceDim; ++i)
	msg << "  " << coordsVertex[i];
      msg << ") in friction model " << label() << "\n"
	  << "using spatial database '" << _dbInitialTract->label() << "'.";
      throw std::runtime_error(msg.str());
    } // if
    _normalizer->nondimensionalize(&initialTractionsVertex[0],
				   initialTractionsVertex.size(), 
				   pressureScale);

    // Rotate tractions from fault coordinate system to global
    // coordinate system
    initialTractionsVertexGlobal = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	initialTractionsVertexGlobal[iDim] += 
	  orientationVertex[jDim*spaceDim+iDim] *
	  initialTractionsVertex[jDim];
      } // for
    } // for

    assert(initialTractionsVertexGlobal.size() ==
	   initialTractionsSection->getFiberDimension(v_fault));
    initialTractionsSection->updatePoint(v_fault, 
					 &initialTractionsVertexGlobal[0]);
  } // for

  // Close properties database
  _dbInitialTract->close();

  initialTractions.complete(); // Assemble contributions

  //initalTractions.view("INITIAL FORCES"); // DEBUGGING
} // _setupInitialTractions

// ----------------------------------------------------------------------
// Compute tractions on fault surface using solution.
void
pylith::faults::FaultCohesiveDyn::_calcTractions(
    topology::Field<topology::SubMesh>* tractions,
    const topology::Field<topology::Mesh>& dispT)
{ // _calcTractions
  assert(0 != tractions);
  assert(0 != _faultMesh);
  assert(0 != _fields);
  assert(0 != _normalizer);

  // Fiber dimension of tractions matches spatial dimension.
  const int spaceDim = _quadrature->spaceDim();
  double_array tractionsVertex(spaceDim);

  // Get sections.
  const ALE::Obj<RealSection>& dispTSection = dispT.section();
  assert(!dispTSection.isNull());

  const ALE::Obj<RealSection>& orientationSection = 
    _fields->get("orientation").section();
  assert(!orientationSection.isNull());

  // Allocate buffer for tractions field (if necessary).
  const ALE::Obj<RealSection>& tractionsSection = tractions->section();
  if (tractionsSection.isNull()) {
    ALE::MemoryLogger& logger = ALE::MemoryLogger::singleton();
    //logger.stagePush("Fault");

    const topology::Field<topology::SubMesh>& dispRel = 
      _fields->get("relative disp");
    tractions->cloneSection(dispRel);

    //logger.stagePop();
  } // if
  const double pressureScale = _normalizer->pressureScale();
  tractions->label("traction");
  tractions->scale(pressureScale);
  tractions->zero();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;

    assert(spaceDim == dispTSection->getFiberDimension(v_lagrange));
    const double* dispTVertex = dispTSection->restrictPoint(v_lagrange);
    assert(dispTVertex);

    assert(spaceDim*spaceDim == 
	   orientationSection->getFiberDimension(v_fault));
    const double* orientationVertex = 
      orientationSection->restrictPoint(v_fault);
    assert(orientationVertex);

    // Rotate tractions to fault coordinate system.
    tractionsVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
	tractionsVertex[iDim] += orientationVertex[iDim*spaceDim+jDim] *
	  dispTVertex[jDim];
      } // for
    } // for

    assert(tractionsVertex.size() == 
	   tractionsSection->getFiberDimension(v_fault));
    tractionsSection->updatePoint(v_fault, &tractionsVertex[0]);
  } // for

  PetscLogFlops(numVertices * (1 + spaceDim) );

#if 0 // DEBUGGING
  tractions->view("TRACTIONS");
#endif

} // _calcTractions

// ----------------------------------------------------------------------
// Update slip rate associated with Lagrange vertex k corresponding
// to diffential velocity between conventional vertices i and j.
void
pylith::faults::FaultCohesiveDyn::_updateVelRel(const topology::SolutionFields& fields)
{ // _updateVelRel
  assert(0 != _fields);

  const int spaceDim = _quadrature->spaceDim();

  // Get section information
  const ALE::Obj<RealSection>& velocitySection =
      fields.get("velocity(t)").section();
  assert(!velocitySection.isNull());

  double_array velRelVertex(spaceDim);
  const ALE::Obj<RealSection>& velRelSection =
      _fields->get("relative velocity").section();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get values
    const double* velocityVertexN = velocitySection->restrictPoint(v_negative);
    assert(0 != velocityVertexN);
    assert(spaceDim == velocitySection->getFiberDimension(v_negative));

    const double* velocityVertexP = velocitySection->restrictPoint(v_positive);
    assert(0 != velocityVertexP);
    assert(spaceDim == velocitySection->getFiberDimension(v_positive));

    // Compute relative velocity
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      const double value = velocityVertexP[iDim] - velocityVertexN[iDim];
      velRelVertex[iDim] = fabs(value) > _zeroTolerance ? value : 0.0;
    } // for

    // Update slip rate field.
    assert(velRelVertex.size() == 
	   velRelSection->getFiberDimension(v_fault));
    velRelSection->updatePoint(v_fault, &velRelVertex[0]);
  } // for

  PetscLogFlops(numVertices*spaceDim*spaceDim*4);
} // _updateVelRel

// ----------------------------------------------------------------------
// Setup sensitivity problem to compute change in slip given change in Lagrange multipliers.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySetup(const topology::Jacobian& jacobian)
{ // _sensitivitySetup
  assert(0 != _fields);
  assert(0 != _quadrature);

  const int spaceDim = _quadrature->spaceDim();

  // Setup fields involved in sensitivity solve.
  if (!_fields->hasField("sensitivity solution")) {
    _fields->add("sensitivity solution", "sensitivity_soln");
    topology::Field<topology::SubMesh>& solution =
        _fields->get("sensitivity solution");
    const topology::Field<topology::SubMesh>& dispRel =
        _fields->get("relative disp");
    solution.cloneSection(dispRel);
    solution.createScatter(solution.mesh());
  } // if
  const topology::Field<topology::SubMesh>& solution =
      _fields->get("sensitivity solution");

  if (!_fields->hasField("sensitivity residual")) {
    _fields->add("sensitivity residual", "sensitivity_residual");
    topology::Field<topology::SubMesh>& residual =
        _fields->get("sensitivity residual");
    residual.cloneSection(solution);
    residual.createScatter(solution.mesh());
  } // if

  if (!_fields->hasField("sensitivity dispRel")) {
    _fields->add("sensitivity relative disp", "sensitivity_relative_disp");
    topology::Field<topology::SubMesh>& dispRel =
        _fields->get("sensitivity relative disp");
    dispRel.cloneSection(solution);
  } // if
  topology::Field<topology::SubMesh>& dispRel =
    _fields->get("sensitivity relative disp");
  dispRel.zero();

  if (!_fields->hasField("sensitivity dLagrange")) {
    _fields->add("sensitivity dLagrange", "sensitivity_dlagrange");
    topology::Field<topology::SubMesh>& dLagrange =
        _fields->get("sensitivity dLagrange");
    dLagrange.cloneSection(solution);
  } // if
  topology::Field<topology::SubMesh>& dLagrange =
    _fields->get("sensitivity dLagrange");
  dLagrange.zero();

  // Setup Jacobian sparse matrix for sensitivity solve.
  if (0 == _jacobian)
    _jacobian = new topology::Jacobian(solution, jacobian.matrixType());
  assert(0 != _jacobian);
  _jacobian->zero();

  // Setup PETSc KSP linear solver.
  if (0 == _ksp) {
    PetscErrorCode err = 0;
    err = KSPCreate(_faultMesh->comm(), &_ksp); CHECK_PETSC_ERROR(err);
    err = KSPSetInitialGuessNonzero(_ksp, PETSC_FALSE); CHECK_PETSC_ERROR(err);
    double rtol = 0.0;
    double atol = 0.0;
    double dtol = 0.0;
    int maxIters = 0;
    err = KSPGetTolerances(_ksp, &rtol, &atol, &dtol, &maxIters); 
    CHECK_PETSC_ERROR(err);
    rtol = _zeroTolerance;
    atol = 0.001*_zeroTolerance;
    err = KSPSetTolerances(_ksp, rtol, atol, dtol, maxIters);
    CHECK_PETSC_ERROR(err);

    PC pc;
    err = KSPGetPC(_ksp, &pc); CHECK_PETSC_ERROR(err);
    err = PCSetType(pc, PCJACOBI); CHECK_PETSC_ERROR(err);
    err = KSPSetType(_ksp, KSPGMRES); CHECK_PETSC_ERROR(err);

    err = KSPAppendOptionsPrefix(_ksp, "friction_");
    err = KSPSetFromOptions(_ksp); CHECK_PETSC_ERROR(err);
  } // if
} // _sensitivitySetup

// ----------------------------------------------------------------------
// Update the Jacobian values for the sensitivity solve.
void
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateJacobian(const bool negativeSide,
                                                             const topology::Jacobian& jacobian,
                                                             const topology::SolutionFields& fields)
{ // _sensitivityUpdateJacobian
  assert(0 != _quadrature);
  assert(0 != _fields);

  const int numBasis = _quadrature->numBasis();
  const int spaceDim = _quadrature->spaceDim();
  const int subnrows = numBasis*spaceDim;
  const int submatrixSize = subnrows * subnrows;

  // Get solution field
  const topology::Field<topology::Mesh>& solutionDomain = fields.solution();
  const ALE::Obj<RealSection>& solutionDomainSection = solutionDomain.section();
  assert(!solutionDomainSection.isNull());

  // Get cohesive cells
  const ALE::Obj<SieveMesh>& sieveMesh = fields.mesh().sieveMesh();
  assert(!sieveMesh.isNull());
  const ALE::Obj<SieveMesh::label_sequence>& cellsCohesive =
    sieveMesh->getLabelStratum("material-id", id());
  assert(!cellsCohesive.isNull());
  const SieveMesh::label_sequence::iterator cellsCohesiveBegin =
    cellsCohesive->begin();
  const SieveMesh::label_sequence::iterator cellsCohesiveEnd =
    cellsCohesive->end();

  // Visitor for Jacobian matrix associated with domain.
  double_array jacobianSubCell(submatrixSize);
  const PetscMat jacobianDomainMatrix = jacobian.matrix();
  assert(0 != jacobianDomainMatrix);
  const ALE::Obj<SieveMesh::order_type>& globalOrderDomain =
    sieveMesh->getFactory()->getGlobalOrder(sieveMesh, "default", solutionDomainSection);
  assert(!globalOrderDomain.isNull());
  const ALE::Obj<SieveMesh::sieve_type>& sieve = sieveMesh->getSieve();
  assert(!sieve.isNull());
  ALE::ISieveVisitor::NConeRetriever<SieveMesh::sieve_type> ncV(*sieve,
      (size_t) pow(sieve->getMaxConeSize(), std::max(0, sieveMesh->depth())));
  int_array indicesGlobal(subnrows);

  // Get fault Sieve mesh
  const ALE::Obj<SieveSubMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());

  // Get sensitivity solution field
  const ALE::Obj<RealSection>& solutionFaultSection =
    _fields->get("sensitivity solution").section();
  assert(!solutionFaultSection.isNull());

  // Visitor for Jacobian matrix associated with fault.
  assert(0 != _jacobian);
  const PetscMat jacobianFaultMatrix = _jacobian->matrix();
  assert(0 != jacobianFaultMatrix);
  const ALE::Obj<SieveSubMesh::order_type>& globalOrderFault =
    faultSieveMesh->getFactory()->getGlobalOrder(faultSieveMesh, "default", solutionFaultSection);
  assert(!globalOrderFault.isNull());
  // We would need to request unique points here if we had an interpolated mesh
  IndicesVisitor jacobianFaultVisitor(*solutionFaultSection,
				      *globalOrderFault,
				      (int) pow(faultSieveMesh->getSieve()->getMaxConeSize(),
						faultSieveMesh->depth())*spaceDim);

  const int iCone = (negativeSide) ? 0 : 1;
  for (SieveMesh::label_sequence::iterator c_iter=cellsCohesiveBegin;
       c_iter != cellsCohesiveEnd;
       ++c_iter) {
    // Get cone for cohesive cell
    ncV.clear();
    ALE::ISieveTraversal<SieveMesh::sieve_type>::orientedClosure(*sieve,
								 *c_iter, ncV);
    const int coneSize = ncV.getSize();
    assert(coneSize == 3*numBasis);
    const SieveMesh::point_type *cohesiveCone = ncV.getPoints();
    assert(0 != cohesiveCone);

    const SieveMesh::point_type c_fault = _cohesiveToFault[*c_iter];
    jacobianSubCell = 0.0;

    // Get indices
    for (int iBasis = 0; iBasis < numBasis; ++iBasis) {
      // negative side of the fault: iCone=0
      // positive side of the fault: iCone=1
      const int v_domain = cohesiveCone[iCone*numBasis+iBasis];
      
      for (int iDim=0, iB=iBasis*spaceDim; iDim < spaceDim; ++iDim) {
	if (globalOrderDomain->isLocal(v_domain))
	  indicesGlobal[iB+iDim] = globalOrderDomain->getIndex(v_domain) + iDim;
	else
	  indicesGlobal[iB+iDim] = -1;

	// Set matrix diagonal entries to 1.0 (used when vertex is not
	// local).  This happens if a vertex is not on the same
	// processor as the cohesive cell.
	jacobianSubCell[(iB+iDim)*numBasis*spaceDim+iB+iDim] = 1.0;
      } // for
    } // for
    
    PetscErrorCode err = MatGetValues(jacobianDomainMatrix, 
				      indicesGlobal.size(), &indicesGlobal[0],
				      indicesGlobal.size(), &indicesGlobal[0],
				      &jacobianSubCell[0]);
    CHECK_PETSC_ERROR_MSG(err, "Restrict from PETSc Mat failed.");

    // Insert cell contribution into PETSc Matrix
    jacobianFaultVisitor.clear();
    err = updateOperator(jacobianFaultMatrix, *faultSieveMesh->getSieve(),
			 jacobianFaultVisitor, c_fault,
			 &jacobianSubCell[0], INSERT_VALUES);
    CHECK_PETSC_ERROR_MSG(err, "Update to PETSc Mat failed.");
  } // for

  _jacobian->assemble("final_assembly");

  //_jacobian->view(); // DEBUGGING
} // _sensitivityUpdateJacobian

// ----------------------------------------------------------------------
// Reform residual for sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivityReformResidual(const bool negativeSide)
{ // _sensitivityReformResidual
  /** Compute residual -L^T dLagrange
   *
   * Note: We need all entries for L, even those on other processors,
   * so we compute L rather than extract entries from the Jacoiab.
   */

  const double signFault = (negativeSide) ? 1.0 : -1.0;

  // Get cell information
  const int numQuadPts = _quadrature->numQuadPts();
  const double_array& quadWts = _quadrature->quadWts();
  assert(quadWts.size() == numQuadPts);
  const int spaceDim = _quadrature->spaceDim();
  const int numBasis = _quadrature->numBasis();


  double_array basisProducts(numBasis*numBasis);

  // Get fault cell information
  const ALE::Obj<SieveMesh>& faultSieveMesh = _faultMesh->sieveMesh();
  assert(!faultSieveMesh.isNull());
  const ALE::Obj<SieveSubMesh::label_sequence>& cells =
    faultSieveMesh->heightStratum(0);
  assert(!cells.isNull());
  const SieveSubMesh::label_sequence::iterator cellsBegin = cells->begin();
  const SieveSubMesh::label_sequence::iterator cellsEnd = cells->end();
  const int numCells = cells->size();

  // Get sections
  double_array coordinatesCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& coordinates = 
    faultSieveMesh->getRealSection("coordinates");
  assert(!coordinates.isNull());
  RestrictVisitor coordsVisitor(*coordinates, 
				coordinatesCell.size(), &coordinatesCell[0]);

  double_array dLagrangeCell(numBasis*spaceDim);
  const ALE::Obj<RealSection>& dLagrangeSection = 
    _fields->get("sensitivity dLagrange").section();
  assert(!dLagrangeSection.isNull());
  RestrictVisitor dLagrangeVisitor(*dLagrangeSection, 
				   dLagrangeCell.size(), &dLagrangeCell[0]);

  double_array residualCell(numBasis*spaceDim);
  topology::Field<topology::SubMesh>& residual =
      _fields->get("sensitivity residual");
  const ALE::Obj<RealSection>& residualSection = residual.section();
  UpdateAddVisitor residualVisitor(*residualSection, &residualCell[0]);

  residual.zero();

  // Loop over cells
  for (SieveSubMesh::label_sequence::iterator c_iter=cellsBegin;
       c_iter != cellsEnd;
       ++c_iter) {
    // Compute geometry
    coordsVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, coordsVisitor);
    _quadrature->computeGeometry(coordinatesCell, *c_iter);

    // Restrict input fields to cell
    dLagrangeVisitor.clear();
    faultSieveMesh->restrictClosure(*c_iter, dLagrangeVisitor);

    // Get cell geometry information that depends on cell
    const double_array& basis = _quadrature->basis();
    const double_array& jacobianDet = _quadrature->jacobianDet();

    // Compute product of basis functions.
    // Want values summed over quadrature points
    basisProducts = 0.0;
    for (int iQuad=0; iQuad < numQuadPts; ++iQuad) {
      const double wt = quadWts[iQuad] * jacobianDet[iQuad];

      for (int iBasis=0, iQ=iQuad*numBasis; iBasis < numBasis; ++iBasis) {
        const double valI = wt*basis[iQ+iBasis];
	
        for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	  
	  basisProducts[iBasis*numBasis+jBasis] += valI*basis[iQ+jBasis];
	} // for
      } // for
    } // for

    residualCell = 0.0;
    
    for (int iBasis=0; iBasis < numBasis; ++iBasis) {
      for (int jBasis=0; jBasis < numBasis; ++jBasis) {
	const double l = signFault * basisProducts[iBasis*numBasis+jBasis];
	for (int iDim=0; iDim < spaceDim; ++iDim) {
	  residualCell[iBasis*spaceDim+iDim] += 
	    l * dLagrangeCell[jBasis*spaceDim+iDim];
	} // for
      } // for
    } // for

    // Assemble cell contribution into field
    residualVisitor.clear();
    faultSieveMesh->updateClosure(*c_iter, residualVisitor);    
  } // for
} // _sensitivityReformResidual

// ----------------------------------------------------------------------
// Solve sensitivity problem.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolve(void)
{ // _sensitivitySolve
  assert(0 != _fields);
  assert(0 != _jacobian);
  assert(0 != _ksp);

  const topology::Field<topology::SubMesh>& residual =
      _fields->get("sensitivity residual");
  const topology::Field<topology::SubMesh>& solution =
      _fields->get("sensitivity solution");

  // Update PetscVector view of field.
  residual.scatterSectionToVector();

  PetscErrorCode err = 0;
  const PetscMat jacobianMat = _jacobian->matrix();
  err = KSPSetOperators(_ksp, jacobianMat, jacobianMat,
    DIFFERENT_NONZERO_PATTERN); CHECK_PETSC_ERROR(err);

  const PetscVec residualVec = residual.vector();
  const PetscVec solutionVec = solution.vector();
  err = KSPSolve(_ksp, residualVec, solutionVec); CHECK_PETSC_ERROR(err);

  // Update section view of field.
  solution.scatterVectorToSection();

#if 0 // DEBUGGING
  residual.view("SENSITIVITY RESIDUAL");
  solution.view("SENSITIVITY SOLUTION");
#endif
} // _sensitivitySolve

// ----------------------------------------------------------------------
// Update the relative displacement field values based on the
// sensitivity solve.
void
pylith::faults::FaultCohesiveDyn::_sensitivityUpdateSoln(const bool negativeSide)
{ // _sensitivityUpdateSoln
  assert(0 != _fields);
  assert(0 != _quadrature);

  const int spaceDim = _quadrature->spaceDim();

  double_array dispVertex(spaceDim);
  const ALE::Obj<RealSection>& solutionSection =
      _fields->get("sensitivity solution").section();
  const ALE::Obj<RealSection>& dispRelSection =
    _fields->get("sensitivity relative disp").section();

  const double sign = (negativeSide) ? -1.0 : 1.0;

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;

    solutionSection->restrictPoint(v_fault, &dispVertex[0], dispVertex.size());

    dispVertex *= sign;

    assert(dispVertex.size() == dispRelSection->getFiberDimension(v_fault));
    dispRelSection->updateAddPoint(v_fault, &dispVertex[0]);
  } // for
} // _sensitivityUpdateSoln

// ----------------------------------------------------------------------
// Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 1-D.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped1D(
                                     double_array* slip,
				     const double_array& dLagrangeTpdt,
				     const double_array& jacobianN,
				     const double_array& jacobianP)
{ // _sensitivitySolveLumped1D
  assert(0 != slip);

  // Sensitivity of slip to changes in the Lagrange multipliers
  const double Spp = 1.0 / jacobianN[0] + 1.0
    / jacobianP[0];

  const double dlp = dLagrangeTpdt[0];
  (*slip)[0] -= Spp * dlp;

  PetscLogFlops(2);
} // _sensitivitySolveLumped1D

// ----------------------------------------------------------------------
// Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 2-D.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped2D(
                                     double_array* slip,
				     const double_array& dLagrangeTpdt,
				     const double_array& jacobianN,
				     const double_array& jacobianP)
{ // _sensitivitySolveLumped2D
  assert(0 != slip);

  // Sensitivity of slip to changes in the Lagrange multipliers
  // Spp = Sqq = 1.0/Aix + 1.0/Ajx
  assert(jacobianN[0] > 0.0);
  assert(jacobianP[0] > 0.0);

  // Check to make sure Jacobian is same at all DOF for
  // vertices i and j (means S is diagonal with equal enties).
  assert(jacobianN[0] == jacobianN[1]);
  assert(jacobianP[0] == jacobianP[1]);
  
  const double Spp = 1.0 / jacobianN[0] + 1.0 / jacobianP[0];
  const double Sqq = Spp;

  const double dlp = dLagrangeTpdt[0];
  const double dlq = dLagrangeTpdt[1];
  (*slip)[0] -= Spp * dlp;
  (*slip)[1] -= Sqq * dlq;

  PetscLogFlops(7);
} // _sensitivitySolveLumped2D

// ----------------------------------------------------------------------
// Solve slip/Lagrange multiplier sensitivity problem for case of lumped Jacobian in 3-D.
void
pylith::faults::FaultCohesiveDyn::_sensitivitySolveLumped3D(
                                     double_array* slip,
				     const double_array& dLagrangeTpdt,
				     const double_array& jacobianN,
				     const double_array& jacobianP)
{ // _sensitivitySolveLumped3D
  assert(0 != slip);

  // Sensitivity of slip to changes in the Lagrange multipliers
  // Aixjx = 1.0/Aix + 1.0/Ajx
  assert(jacobianN[0] > 0.0);
  assert(jacobianP[0] > 0.0);

  // Check to make sure Jacobian is same at all DOF for
  // vertices i and j (means S is diagonal with equal enties).
  assert(jacobianN[0] == jacobianN[1] && 
	 jacobianN[0] == jacobianN[2]);
  assert(jacobianP[0] == jacobianP[1] && 
	 jacobianP[0] == jacobianP[2]);


  const double Spp = 1.0 / jacobianN[0] + 1.0 / jacobianP[0];
  const double Sqq = Spp;
  const double Srr = Spp;

  const double dlp = dLagrangeTpdt[0];
  const double dlq = dLagrangeTpdt[1];
  const double dlr = dLagrangeTpdt[2];
  (*slip)[0] -= Spp * dlp;
  (*slip)[1] -= Sqq * dlq;
  (*slip)[2] -= Srr * dlr;

  PetscLogFlops(9);
} // _sensitivitySolveLumped3D

// ----------------------------------------------------------------------
// Constrain solution space in 1-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace1D(double_array* dLagrangeTpdt,
         const double_array& slip,
         const double_array& sliprate,
         const double_array& tractionTpdt)
{ // _constrainSolnSpace1D
  assert(0 != dLagrangeTpdt);

    if (tractionTpdt[0] < 0) {
      // if compression, then no changes to solution
    } else {
      // if tension, then traction is zero.

#if 0 // :TODO: FIX THIS
      const double dlp = -tractionTpdt[0] * area;
#else
      const double dlp = -tractionTpdt[0];
#endif // 
      (*dLagrangeTpdt)[0] = dlp;
    } // else

    PetscLogFlops(2);
} // _constrainSolnSpace1D

// ----------------------------------------------------------------------
// Constrain solution space in 2-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D(double_array* dLagrangeTpdt,
         const double_array& slip,
         const double_array& slipRate,
         const double_array& tractionTpdt)
{ // _constrainSolnSpace2D
  assert(0 != dLagrangeTpdt);

  const double slipMag = fabs(slip[0]);
  const double slipRateMag = fabs(slipRate[0]);

  const double tractionNormal = tractionTpdt[1];
  const double tractionShearMag = fabs(tractionTpdt[0]);

  if (tractionNormal < 0 && 0.0 == slip[1]) {
    // if in compression and no opening
    const double frictionStress = _friction->calcFriction(slipMag, slipRateMag,
							  tractionNormal);
    if (tractionShearMag > frictionStress || slipRateMag > 0.0) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip

      if (tractionShearMag > 0.0) {
	// Update traction increment based on value required to stick
	// versus friction
	const double dlp = -(tractionShearMag - frictionStress) *
	  tractionTpdt[0] / tractionShearMag;
	(*dLagrangeTpdt)[0] = dlp;
	(*dLagrangeTpdt)[1] = 0.0;
      } else {
	(*dLagrangeTpdt)[0] = -(*dLagrangeTpdt)[0];
	(*dLagrangeTpdt)[1] = 0.0;
      } // if/else
    } else {
      // friction exceeds value necessary to stick
      // no changes to solution
      assert(0.0 == slipRateMag);
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dLagrangeTpdt)[0] = -tractionTpdt[0];
    (*dLagrangeTpdt)[1] = -tractionTpdt[1];
  } // else

  PetscLogFlops(8);
} // _constrainSolnSpace2D

// ----------------------------------------------------------------------
// Constrain solution space in 3-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D(double_array* dLagrangeTpdt,
         const double_array& slip,
         const double_array& slipRate,
         const double_array& tractionTpdt)
{ // _constrainSolnSpace3D
  assert(0 != dLagrangeTpdt);

  const double slipShearMag = sqrt(slip[0] * slip[0] +
             slip[1] * slip[1]);
  double slipRateMag = sqrt(slipRate[0]*slipRate[0] + 
            slipRate[1]*slipRate[1]);
  
  const double tractionNormal = tractionTpdt[2];
  const double tractionShearMag = 
    sqrt(tractionTpdt[0] * tractionTpdt[0] +
	 tractionTpdt[1] * tractionTpdt[1]);
  
  if (tractionNormal < 0.0 && 0.0 == slip[2]) {
    // if in compression and no opening
    const double frictionStress = 
      _friction->calcFriction(slipShearMag, slipRateMag, tractionNormal);
    if (tractionShearMag > frictionStress || slipRateMag > 0.0) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip
      
      if (tractionShearMag > 0.0) {
	// Update traction increment based on value required to stick
	// versus friction
	const double dlp = -(tractionShearMag - frictionStress) * 
	  tractionTpdt[0] / tractionShearMag;
	const double dlq = -(tractionShearMag - frictionStress) * 
	  tractionTpdt[1] / tractionShearMag;
	
	(*dLagrangeTpdt)[0] = dlp;
	(*dLagrangeTpdt)[1] = dlq;
	(*dLagrangeTpdt)[2] = 0.0;
      } else {
	(*dLagrangeTpdt)[0] = -(*dLagrangeTpdt)[0];
	(*dLagrangeTpdt)[0] = -(*dLagrangeTpdt)[0];
	(*dLagrangeTpdt)[2] = 0.0;
      } // if/else	
      
    } else {
      // else friction exceeds value necessary, so stick
      // no changes to solution
      assert(0.0 == slipRateMag);
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dLagrangeTpdt)[0] = -tractionTpdt[0];
    (*dLagrangeTpdt)[1] = -tractionTpdt[1];
    (*dLagrangeTpdt)[2] = -tractionTpdt[2];
  } // else

  PetscLogFlops(22);
} // _constrainSolnSpace3D


// End of file 
