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

#include "pylith/faults/FrictionSlipWeakening.hh" // implementation of object methods

#include "pylith/faults/AuxiliaryFactoryRheology.hh" // USES AuxiliaryFactoryRheology
#include "pylith/fekernels/FrictionSlipWeakening.hh" // USES FrictionSlipWeakening kernels
#include "pylith/utils/error.hh"    // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_DEBUG

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FrictionSlipWeakening::FrictionSlipWeakening(void) :
    _auxiliaryFactory(new pylith::faults::AuxiliaryFactoryRheology) {
    pylith::utils::PyreComponent::setName("frictionslipweakening");
} // constructor

// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FrictionSlipWeakening::~FrictionSlipWeakening(void) {
    deallocate();
} // destructor

// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FrictionSlipWeakening::deallocate(void)  {
    FaultFriction::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate

// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::faults::AuxiliaryFactoryRheology*
pylith::faults::FrictionSlipWeakening::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory

// ------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::faults::FrictionSlipWeakening::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).
    _auxiliaryFactory->addStaticCoefficient();
    _auxiliaryFactory->addDynamicCoefficient();
    _auxiliaryFactory->addSlipWeakeningParameter();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields

// ------------------------------------------------------------------------------------------------
// return residual kernal f0u
PetscBdPointFunc 
pylith::faults::FrictionSlipWeakening::getF0uKernel(void) const {
    PYLITH_METHOD_BEGIN;
    const PetscBdPointFunc fu0 = pylith::fekernels::FrictionSlipWeakening::fu0;
    return fu0;
}

// ------------------------------------------------------------------------------------------------
// return residual kernal Jf0uu
PetscBdPointJac 
pylith::faults::FrictionSlipWeakening::getJf0uuKernel(void) const {
    PYLITH_METHOD_BEGIN;
    const PetscBdPointJac Jf0uu = pylith::fekernels::FrictionSlipWeakening::Jf0uu;
    return Jf0uu;
}

// End of file
