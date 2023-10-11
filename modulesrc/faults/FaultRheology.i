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

/** @file modulesrc/faults/FaultRheology.i
 *
 * @brief Python interface to C++ FaultRheology object.
 */

namespace pylith {
    namespace faults {
        class FaultRheology : public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            FaultRheology(void);

            /// Destructor.
            virtual ~FaultRheology(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

             /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual
            pylith::faults::AuxiliaryFactoryRheology* getAuxiliaryFactory(void) = 0;

            /// Add rheology subfields to auxiliary field.
            virtual
            void addAuxiliarySubfields(void) = 0;

            /** Get triggers for needing to compute the elastic constants for the RHS Jacobian.
             *
             * @returns Triggers for needing to recompute the RHS Jacobian.
             */
            int getLHSJacobianTriggers(void) const;

            /** Add kernels for updating state variables.
             *
             * @param[inout] kernels Array of kernels for updating state variables.
             * @param[in] coordsys Coordinate system.
             */
            virtual
            void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                        const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Update kernel constants.
             *
             * @param[inout] kernelConstants Array of constants used in integration kernels.
             * @param[in] dt Current time step.
             */
            virtual
            void updateKernelConstants(pylith::real_array* kernelConstants,
                                    const PylithReal dt) const;
            
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


        }; // class FaultRheology

    } // faults
} // pylith

// End of file
