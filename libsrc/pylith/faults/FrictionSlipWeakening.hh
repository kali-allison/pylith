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

#if !defined(pylith_faults_frictionslipweakening_hh)
#define pylith_faults_frictionslipweakeningc_hh

#include "faultsfwd.hh"                   // forward declarations
#include "pylith/faults/FaultFriction.hh" // ISA FaultFriction

class pylith::faults::FrictionSlipWeakening : public pylith::faults::FaultFriction
{
    friend class TestStaticFriction; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    FrictionSlipWeakening(void);

    /// Destructor.
    ~FrictionSlipWeakening(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::faults::AuxiliaryFactoryRheology *getAuxiliaryFactory(void);

    /// Add rheology subfields to auxiliary field.
    void addAuxiliarySubfields(void);

    /// Get f0u for positive side of fault.
    PetscBdPointFunc getF0uKernel(void) const;

    /// Get Jf0uu for positive side of fault.
    PetscBdPointJac getJf0uuKernel(void) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    pylith::faults::AuxiliaryFactoryRheology* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.


    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:
    FrictionSlipWeakening(const FrictionSlipWeakening &);                  ///< Not implemented.
    const FrictionSlipWeakening &operator=(const FrictionSlipWeakening &); /// Not implemented.

}; // class FrictionSlipWeakening

#endif // pylith_faults_frictionslipweakening_hh

// End of file
