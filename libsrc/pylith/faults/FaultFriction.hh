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

#if !defined(pylith_faults_faultfriction_hh)
#define pylith_faults_faultfriction_hh

#include "faultsfwd.hh"                   // forward declarations
#include "pylith/faults/FaultRheology.hh" // ISA FaultRheology

class pylith::faults::FaultFriction : public pylith::faults::FaultRheology
{
    friend class TestStaticFriction; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    FaultFriction(void);

    /// Destructor.
    virtual ~FaultFriction(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual pylith::faults::AuxiliaryFactoryRheology *getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual void addAuxiliarySubfields(pylith::faults::AuxiliaryFactoryRheology*);

    /// Get f0u for negative side of fault.
    virtual
    PetscBdPointFunc getF0uNegKernel(void) const = 0;

    /// Get f0u for positive side of fault.
    virtual
    PetscBdPointFunc getF0uPosKernel(void) const = 0;

    /// Get Jf0uu for negative side of fault.
    virtual
    PetscBdPointJac getJf0uuNegKernel(void) const = 0;

    /// Get Jf0uu for positive side of fault.
    virtual
    PetscBdPointJac getJf0uuPosKernel(void) const = 0;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::faults::AuxiliaryFactoryRheology* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:
    FaultFriction(const FaultFriction &);                  ///< Not implemented.
    const FaultFriction &operator=(const FaultFriction &); /// Not implemented.

}; // class FaultFriction

#endif // pylith_faults_faultfriction_hh

// End of file
