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

/** @file modulesrc/faults/FrictionStatic.i
 *
 * @brief Python interface to C++ FrictionStatic object.
 */

namespace pylith {
    namespace faults {
        class FrictionStatic :  public pylith::faults::FaultFriction {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            FrictionStatic(void);

            /// Destructor.
            virtual ~FrictionStatic(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual pylith::faults::AuxiliaryFactoryRheology *getAuxiliaryFactory(void) = 0;

            /// Add rheology subfields to auxiliary field.
            virtual void addAuxiliarySubfields(void);

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

        }; // class FrictionStatic

    } // faults
} // pylith

// End of file
