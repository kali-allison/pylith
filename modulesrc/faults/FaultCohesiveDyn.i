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

/** @file modulesrc/faults/FaultCohesiveDyn.i
 *
 * @brief Python interface to C++ FaultCohesiveDyn object.
 */

namespace pylith {
    namespace faults {
        class FaultCohesiveDyn : public pylith::faults::FaultCohesive {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            FaultCohesiveDyn(void);

            /// Destructor.
            virtual ~FaultCohesiveDyn(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set fault rheology.
             *
             * @param[in] rheology Fault rheology.
             */
            void setFaultRheology(pylith::faults::FaultRheology* const rheology);

            /** Get fault rheology.
             *
             * @returns Fault rheology.
             */
            pylith::faults::FaultRheology* getFaultRheology(void) const;

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Set fault traction perturbations.
             *
             * @param names Array of fault traction perturbation names.
             * @param numNames Number of fault traction perturbation names.
             * @param ruptures Array of fault traction perturbations.
             * @param numRuptures Number of fault traction perturbations.
             */
            void setTractionPerturbations(const char* const* names,
                                        const int numNames,
                                        pylith::faults::TractionPerturbation** ruptures,
                                        const int numPerturbations);

            /** Create integrator and set kernels.
             *
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             * @returns Integrator if applicable, otherwise NULL.
             */
            pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution,
                                                             const std::vector<pylith::materials::Material*>& materials);

            /** Create auxiliary field.
             *
             * @param[in] solution Solution field.
             * @param[in\ domainMesh Finite-element mesh associated with integration domain.
             *
             * @returns Auxiliary field if applicable, otherwise NULL.
             */
            pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                          const pylith::topology::Mesh& domainMesh);

            /** Update auxiliary subfields at beginning of time step.
             *
             * @param[out] auxiliaryField Auxiliary field.
             * @param[in] t Current time.
             */
            void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                      const double t);

            // PROTECTED METHODS //////////////////////////////////////////////////////////////////
protected:

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

            /** Update traction perturbation in auxiliary field at beginning of time step.
             *
             * @param[out] auxiliaryField Auxiliary field.
             * @param[in] t Current time.
             */
            void _updateTractionPerturbation(pylith::topology::Field* auxiliaryField,
                                            const double t);

            /** Set kernels for residual.
             *
             * @param[out] integrator Integrator for material.
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             */
            void _setKernelsResidual(pylith::feassemble::IntegratorInterface* integrator,
                                     const pylith::topology::Field& solution,
                                     const std::vector<pylith::materials::Material*>& materials) const;

            /** Set kernels for Jacobian.
             *
             * @param[out] integrator Integrator for material.
             * @param[in] solution Solution field.
             * @param[in] materials Materials in problem.
             */
            void _setKernelsJacobian(pylith::feassemble::IntegratorInterface* integrator,
                                     const pylith::topology::Field& solution,
                                     const std::vector<pylith::materials::Material*>& materials) const;
        
            /** Set kernels for computing updated state variables in auxiliary field.
             *
             * @param[out] integrator Integrator for material.
             * @param[in] solution Solution field.
             */
            void _setKernelsUpdateStateVars(pylith::feassemble::IntegratorInterface* integrator,
                                            const pylith::topology::Field& solution) const;

        }; // class FaultCohesiveDyn

    } // faults
} // pylith

// End of file
