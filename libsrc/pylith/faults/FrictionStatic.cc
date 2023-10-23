// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/faults/FrictionStatic.hh" // implementation of object methods

#include "pylith/faults/AuxiliaryFactoryRheology.hh" // USES AuxiliaryFactoryRheology
#include "pylith/fekernels/FrictionStatic.hh" // USES FrictionStatic kernels
#include "pylith/utils/error.hh"    // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_DEBUG

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FrictionStatic::FrictionStatic(void) :
    _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryRheology) {
    pylith::utils::PyreComponent::setName("frictionstatic");
} // constructor

// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FrictionStatic::~FrictionStatic(void) {
    deallocate();
} // destructor

// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FrictionStatic::deallocate(void)  {
    FaultFriction::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate

// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::faults::AuxiliaryFactoryRheology*
pylith::faults::FrictionStatic::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory

// ------------------------------------------------------------------------------------------------
// Set database for auxiliary fields.
void
pylith::faults::FrictionStatic::setAuxiliaryFieldDB(spatialdata::spatialdb::SpatialDB* value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setAuxiliaryFieldDB(value="<<typeid(value).name()<<")");

    assert(_auxiliaryFactory);
    _auxiliaryFactory->setQueryDB(value);

    PYLITH_METHOD_END;
} // setAuxiliaryFieldDB

// ------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::faults::FrictionStatic::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).
    //FaultFriction::addAuxiliarySubfields(_auxiliaryFactory);
    _auxiliaryFactory->addCohesion();
    _auxiliaryFactory->addStaticCoefficient();
    
    PYLITH_METHOD_END;
} // addAuxiliarySubfields

// ------------------------------------------------------------------------------------------------
// return residual kernal f0u
PetscBdPointFunc
pylith::faults::FrictionStatic::getF0uKernel(void) const {
    return pylith::fekernels::FrictionStatic::f0u;
}

// ------------------------------------------------------------------------------------------------
// return Jacobian kernal Jf0uu
PetscBdPointJac
pylith::faults::FrictionStatic::getJf0uuKernel(void) const {
    return pylith::fekernels::FrictionStatic::Jf0uu;
}


// End of file
