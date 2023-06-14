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

#include "pylith/faults/FaultFriction.hh" // implementation of object methods

#include "pylith/feassemble/Integrator.hh" // USES NEW_JACOBIAN_NEVER

#include "pylith/utils/error.hh"    // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_DEBUG

#include "pylith/feassemble/Integrator.hh" // USES NEW_JACOBIAN_NEVER

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultFriction::FaultFriction(void) {}

// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::FaultFriction::~FaultFriction(void)
{
    deallocate();
} // destructor

// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FaultFriction::deallocate(void) {}

// End of file
