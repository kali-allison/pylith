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
// Add rheology subfields to auxiliary field.
void
pylith::faults::FrictionStatic::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).
    _auxiliaryFactory->addStaticCoefficient();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields

// ------------------------------------------------------------------------------------------------
// return residual kernal f0u_neg
PetscBdPointFunc
pylith::faults::FrictionStatic::getF0uNegKernel(void) const {
    PYLITH_METHOD_BEGIN;
    const PetscBdPointFunc f0u_neg = pylith::fekernels::FrictionStatic::f0u_neg;
    return f0u_neg;
}

// ------------------------------------------------------------------------------------------------
// return residual kernal f0u_pos
PetscBdPointFunc 
pylith::faults::FrictionStatic::getF0uPosKernel(void) const {
    PYLITH_METHOD_BEGIN;
    const PetscBdPointFunc f0u_pos = pylith::fekernels::FrictionStatic::f0u_pos;
    return f0u_pos;
}

// ------------------------------------------------------------------------------------------------
// return Jacobian kernal Jf0uu_neg
PetscBdPointJac
pylith::faults::FrictionStatic::getJf0uuNegKernel(void) const {
    PYLITH_METHOD_BEGIN;
    const PetscBdPointJac Jf0uu_neg = pylith::fekernels::FrictionStatic::Jf0uu_neg;
    return Jf0uu_neg;
}

// ------------------------------------------------------------------------------------------------
// return residual kernal Jf0uu_pos
PetscBdPointJac 
pylith::faults::FrictionStatic::getJf0uuPosKernel(void) const {
    PYLITH_METHOD_BEGIN;
    const PetscBdPointJac Jf0uu_pos = pylith::fekernels::FrictionStatic::Jf0uu_pos;
    return Jf0uu_pos;
}

// End of file
