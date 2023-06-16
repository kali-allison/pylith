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

/** @file libsrc/fekernels/FaultRheology.hh
 *
 * Kernels for fault rheologies.
 *
 * Solution fields: [disp(dim), vel(dim, optional), lagrange(dim)]
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

#if !defined(pylith_fekernels_faultfriction_hh)
#define pylith_fekernels_faultfriction_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::FaultFriction {
    // PUBLIC TYPEDEFS ////////////////////////////////////////////////////////////////////////////
public:

    struct SlipContext {
        PylithReal t;
        PylithReal slip;
        PylithReal slipRate;
    };

    struct FrictionContext {
        PylithInt dim;
        PylithReal cohesion;
        PylithReal tractionNormal;
        PylithReal slipRate[2];
        PylithReal opening;
        bool openFreeSurface;
    };

    typedef void (*frictioncoeffn_type)(const SlipContext&,
                                        void*, // rheologyContext
                                        PylithReal*); // tractionFriction

    typedef void (*frictiondirfn_type)(const FrictionContext&,
                                       const PylithReal, // frictionMag
                                       PylithReal*); // tractionFriction

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

// --------------------------------------------------------------------------------------------
// create slip context
    static inline
    void setSlipContext(SlipContext* context,
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
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);
        
        const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

        const PylithInt i_disp = 0;
        const PylithInt sOffDispN = sOff[i_disp];
        const PylithInt sOffDispP = sOffDispN+spaceDim;
        const PylithInt fOffLagrange = 0;

        // displacement on + and - sides of fault
        const PylithScalar* dispN = &s[sOffDispN];
        const PylithScalar* dispP = &s[sOffDispP];

        // rate on + and - sides of fault
        const PylithScalar* velN = &s_t[sOffDispN];
        const PylithScalar* velP = &s_t[sOffDispP];

        context->t = t;
        context->slip = dispP - dispN;
        context->slipRate = velP - velN;
    } // setSlipContext

    // --------------------------------------------------------------------------------------------
// create friction context
    static inline
    void setFrictionContext(FrictionContext* context,
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
                    const PylithInt numConstants,
                    const PylithScalar constants[],
                    const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        const PylithInt i_cohesion = numA-5;
        const PylithInt i_tractionNormal = numA-4;
        const PylithInt i_slipRate = numA-3;
        const PylithInt i_opening = numA-2;
        const PylithInt i_openFreeSurface = numA-1;

        assert(a);
        assert(aOff);
        assert(aOff[i_cohesion] >= 0);
        assert(aOff[i_tractionNormal] >= 0);
        assert(aOff[i_slipRate] >= 0);
        assert(aOff[i_opening] >= 0);
        assert(aOff[i_openFreeSurface] >= 0);

        context->cohesion = a[aOff[i_cohesion]];assert(context->cohesion >= 0.0);
        context->tractionNormal = a[aOff[i_tractionNormal]];assert(context->tractionNormal >= 0.0);
        PylithReal slipRate0 = a[aOff[i_slipRate]];assert(slipRate0 >= 0.0);
        context->slipRate[0] = slipRate0;
        PylithReal slipRate1 = a[aOff[i_slipRate]+1];assert(slipRate1 >= 0.0);
        context->slipRate[1] = slipRate1;
        context->opening = a[aOff[i_opening]];assert(context->opening >= 0.0);
        context->openFreeSurface = a[aOff[i_openFreeSurface]];assert(context->openFreeSurface >= 0.0);

    } // setFrictionContext

    // --------------------------------------------------------------------------------------------
    /** Compute friction traction.
     */
    static inline
    void friction(const SlipContext& slipContext,
                  const FrictionContext& frictionContext,
                  void* rheologyContext,
                  frictioncoeffn_type frictionCoefFn,
                  frictiondirfn_type frictionDirFn,
                  PylithReal* tractionFriction) {
        assert(rheologyContext);
        assert(frictionCoefFn);
        assert(frictionDirFn);
        assert(tractionFriction);

        PylithReal frictionCoef = 0.0;
        frictionCoefFn(slipContext, rheologyContext, &frictionCoef);

        PylithReal frictionMag;
        if (0.0 == frictionContext.opening) {
            frictionMag = frictionContext.cohesion + frictionCoef * frictionContext.tractionNormal;
        } else if (!frictionContext.openFreeSurface) {
            frictionMag = frictionContext.cohesion;
        } else {
            frictionMag = 0.0;
        } // if/else

        frictionDirFn(frictionContext, frictionMag, tractionFriction);
    } // friction

    // --------------------------------------------------------------------------------------------
    /** Compute friction direction in 2D.
     */
    static inline
    void frictionDir2D(const FrictionContext& frictionContext,
                       const PylithReal frictionMag,
                       PylithReal* tractionFriction) {
        assert(tractionFriction);
        assert(frictionContext.slipRate);

        const PylithReal slipRate = frictionContext.slipRate[0];
        tractionFriction[0] = -slipRate/fabs(slipRate) * frictionMag;
        tractionFriction[1] = 0.0;
    } // frictionDir2D

    // --------------------------------------------------------------------------------------------
    /** Compute friction direction in 3D.
     */
    static inline
    void frictionDir3D(const FrictionContext& frictionContext,
                       const PylithReal frictionMag,
                       PylithReal* tractionFriction) {
        assert(tractionFriction);
        assert(frictionContext.slipRate);

        const PylithReal slipRate1 = frictionContext.slipRate[0];
        const PylithReal slipRate2 = frictionContext.slipRate[2];
        const PylithReal slipRateMag = sqrt(slipRate1*slipRate1 + slipRate2*slipRate2);
        tractionFriction[0] = -slipRate1/slipRateMag * frictionMag;
        tractionFriction[1] = -slipRate2/slipRateMag * frictionMag;
        tractionFriction[2] = 0.0;
    } // frictionDir3D

// --------------------------------------------------------------------------------------------
    /** f0 function for elasticity equation: f0u = +\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static inline
    void f0u_neg(const PylithInt dim,
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
                 const PylithReal n[],
                 const PylithInt numConstants,
                 const PylithScalar constants[],
                 PylithScalar f0[]) {
        assert(sOff);
        assert(s);
        assert(f0);

        assert(numS >= 2);

        const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1

        const PylithInt fOffN = 0;
        const PylithInt sOffLagrange = sOff[numS-1];
        const PylithScalar* lagrange = &s[sOffLagrange];

        pylith::fekernels::FaultFriction::SlipContext slipContext;
        pylith::fekernels::FaultFriction::setSlipContext(&slipContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::FaultFriction::???);

        pylith::fekernels::FaultFriction::FrictionContext frictionContext;
        pylith::fekernels::FaultFriction::setFrictionContext(&frictionContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, numConstants, constants, pylith::fekernels::FaultFriction::???);

        PylithReal tractionFriction = 0.0;
        friction(slipContext,frictionContext,rheologyContext,frictionCoefFn, frictionDirFn,&tractionFriction);

        for (PylithInt i = 0; i < spaceDim; ++i) {
            f0[fOffN+i] += tractionFriction[i];
        } // for
    } // f0u_neg

}; // FaultFriction

#endif // pylith_fekernels_faultrheology_hh

/* End of file */
