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

/** @file libsrc/faults/FaultFriction.hh
 *
 * @brief C++ abstract base class for friction constitutive models.
 */

#if !defined(pylith_faults_frictionstatic_hh)
#define pylith_faults_frictionstatic_hh

#include "faultsfwd.hh"                   // forward declarations
#include "pylith/faults/FaultFriction.hh" // ISA FaultFriction

class pylith::faults::FrictionStatic : public pylith::faults::FaultFriction
{
    friend class TestStaticFriction; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    FrictionStatic(void);

    /// Destructor.
    ~FrictionStatic(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::faults::AuxiliaryFactoryRheology *getAuxiliaryFactory(void);

    /// Add rheology subfields to auxiliary field.
    void addAuxiliarySubfields(void);

    /// Get f0u for negative side of fault.
    PetscBdPointFunc getF0uNegKernel(void) const;

    /// Get f0u for positive side of fault.
    PetscBdPointFunc getF0uPosKernel(void) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::faults::AuxiliaryFactoryRheology* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.


    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:
    FrictionStatic(const FrictionStatic &);                  ///< Not implemented.
    const FrictionStatic &operator=(const FrictionStatic &); /// Not implemented.

}; // class FrictionStatic

#endif // pylith_faults_frictionstatic_hh

// End of file
