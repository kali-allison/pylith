/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2022 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/FrictionSlipWeakening`.hh
 *
 * Kernels for static friction on a fault.
 *
 * Solution fields: [disp(dim), vel(dim, optional), lagrange(dim)]
 *
 * Auxiliary fields: [...] (not used)
 *
 * LHS Residual
 *
 *
 * LHS Jacobian
 *
 * Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] n Outward normal vector for boundary.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

#if !defined(pylith_fekernels_frictionslipweakening_hh)
#define pylith_fekernels_frictionslipweakening_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/FaultFriction.hh" // USES FaultFriction

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::FrictionSlipWeakening {
    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal staticCoefficient;
        PylithReal dynamicCoefficient;
        PylithReal slipWeakeningCoefficient;
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    static inline
    void setContext(Context* context,
                    const PylithInt dim,
                    const PylithInt numS,
                    const PylithInt numA,
                    const PylithInt sOff[],
                    const PylithInt sOff_x[],
                    const PylithScalar s[],
                    const PylithScalar s_t[],
                    const PylithScalar s_x[],
                    const PylithInt aOff[],
                    const PylithInt aOff_x[],
                    const PylithScalar a[],
                    const PylithScalar a_t[],
                    const PylithScalar a_x[],
                    const PylithReal t,
                    const PylithScalar x[],
                    const PylithScalar n[],
                    const PylithInt numConstants,
                    const PylithScalar constants[]) {
        assert(context);

        const PylithInt i_static = numA-3;
        const PylithInt i_dynamic = numA-2;
        const PylithInt i_slipWeakening = numA-1;

        assert(a);
        assert(aOff);
        assert(aOff[i_static] >= 0);
        assert(aOff[i_dynamic] >= 0);
        assert(aOff[i_slipWeakening] >= 0);

        context->staticCoefficient = a[aOff[i_static]];assert(context->staticCoefficient >= 0.0);
        context->dynamicCoefficient = a[aOff[i_dynamic]];assert(context->dynamicCoefficient >= 0.0);
        context->slipWeakeningCoefficient = a[aOff[i_slipWeakening]];assert(context->slipWeakeningCoefficient >= 0.0);
    } // setContext

    // --------------------------------------------------------------------------------------------
    /** Compute coefficient of friction.
     *
     */
    static inline
    void frictionCoefficient(PetscReal t,
                             PetscReal slip,
                             PetscReal slipVel,
                             void* rheologyContext,
                             PylithReal* coefficient) {
        assert(rheologyContext);
        assert(coefficient);

        Context* context = (Context*)(rheologyContext);
        PetscReal absSlip = abs(slip);
        if (absSlip < context->slipWeakeningCoefficient) {
          *coefficient = context->staticCoefficient - (context->staticCoefficient - context->dynamicCoefficient) * abs(slip)/context->slipWeakeningCoefficient;
        }
        else {
            *coefficient = context->dynamicCoefficient;
        }
    } // friction_coefficient

// --------------------------------------------------------------------------------------------
    // residual fu0
    static inline
    void fu0(const PylithInt dim,
                const PylithInt numS,
                const PylithInt numA,
                const PylithInt sOff[],
                const PylithInt sOff_x[],
                const PylithScalar s[],
                const PylithScalar s_t[],
                const PylithScalar s_x[],
                const PylithInt aOff[],
                const PylithInt aOff_x[],
                const PylithScalar a[],
                const PylithScalar a_t[],
                const PylithScalar a_x[],
                const PylithReal t,
                const PylithScalar x[],
                const PylithScalar n[],
                const PylithInt numConstants,
                const PylithScalar constants[],
                PylithScalar f0[]) {

        // create FaultFriction context
        pylith::fekernels::FaultFriction::Context frictionContext;
        pylith::fekernels::FaultFriction::setContext(&frictionContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants);

        // create FrictionSlipWeakening context
        pylith::fekernels::FrictionSlipWeakening::Context rheologyContext;
        pylith::fekernels::FrictionSlipWeakening::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants);
        
        pylith::fekernels::FaultFriction::f0u(frictionContext,&rheologyContext,frictionCoefficient,f0);
    } // fu0

// --------------------------------------------------------------------------------------------
    // Jacobian Jf0uu on positive side of fault: Jf0uu = 0
    static inline
    void Jf0uu(const PylithInt dim,
                   const PylithInt numS,
                   const PylithInt numA,
                   const PylithInt sOff[],
                   const PylithInt sOff_x[],
                   const PylithScalar s[],
                   const PylithScalar s_t[],
                   const PylithScalar s_x[],
                   const PylithInt aOff[],
                   const PylithInt aOff_x[],
                   const PylithScalar a[],
                   const PylithScalar a_t[],
                   const PylithScalar a_x[],
                   const PylithReal t,
                   const PylithReal s_tshift,
                   const PylithScalar x[],
                   const PylithReal n[],
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar Jf0[]) {
        assert(numS >= 2);
        assert(Jf0);
        assert(sOff);
        assert(n);

        // create FaultFriction context
        pylith::fekernels::FaultFriction::Context frictionContext;
        pylith::fekernels::FaultFriction::setContext(&frictionContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants);

        // create FrictionSlipWeakening context
        pylith::fekernels::FrictionSlipWeakening::Context rheologyContext;
        pylith::fekernels::FrictionSlipWeakening::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants);


        const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

        const PylithInt ncols = spaceDim;
        const PylithInt gOffN = 0;
        const PylithInt gOffP = gOffN+spaceDim*spaceDim;

        // compute slip in fault coords
        PylithReal slipFaultCoords[2];
        pylith::fekernels::FaultFriction::computeDiffFaultCoords(slipFaultCoords,spaceDim, frictionContext.dispN, frictionContext.dispP,frictionContext.n,frictionContext.refDir);
        PylithReal slipMagTangent = (2 == spaceDim) ? fabs(slipFaultCoords[0]) : sqrt(pow(slipFaultCoords[0],2) + pow(slipFaultCoords[1],2));

        // compute traction normal to fault
        PylithReal lagrangeMultFaultCoords[2];
        pylith::fekernels::BoundaryDirections::toTN(lagrangeMultFaultCoords,frictionContext.lagrangeMultGlobalCoords,frictionContext.n);
        PylithReal tractionNormalFaultCoords = (2 == spaceDim) ? -lagrangeMultFaultCoords[1] : -lagrangeMultFaultCoords[2];

        PetscReal dtau_dslip = (rheologyContext.staticCoefficient - rheologyContext.dynamicCoefficient) * tractionNormalFaultCoords / rheologyContext.slipWeakeningCoefficient;
        if (slipMagTangent >= rheologyContext.slipWeakeningCoefficient) {
            dtau_dslip = 0.0;
        }
        for (PylithInt i = 0; i < spaceDim; ++i) {
            Jf0[gOffN+i*ncols+i] += dtau_dslip;
            Jf0[gOffP+i*ncols+i] += dtau_dslip;
        } // for
    } // Jf0uu

}; // FrictionSlipWeakening

#endif // pylith_fekernels_frictionslipweakening_hh

/* End of file */
