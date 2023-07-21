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

#include "pylith/fekernels/BoundaryDirections.hh" // USES BoundaryDirections

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::FaultFriction {
    // PUBLIC TYPEDEFS ////////////////////////////////////////////////////////////////////////////
public:


    struct Context {
        PylithInt dim;
        PylithScalar t; // time
        PylithScalar dt; // time step
        const PylithScalar* n; // normal to fault
        const PylithScalar* refDir; // reference direction
        const PylithScalar* tractionGlobal; // global coordinates
        const PylithScalar* tractionNormalGlobal; // global coordinates
        PylithScalar* tractionShearFault; // fault coordinates
        PylithScalar tractionNormalFault; // fault coordinates
        PylithScalar tractionShearMagFault; // fault coordinates, magnitude of friction
        const PylithScalar* dispN; // (spaceDim) displacement negative side of fault
        const PylithScalar* dispP; // (spaceDim) displacement positive side of fault
        const PylithScalar* velN; // (spaceDim) velocity negative side of fault
        const PylithScalar* velP; // (spaceDim) velocity positive side of fault
        PylithScalar* slip; // (spaceDim) fault coordinates
        PylithScalar* slipVel; // (spaceDim) fault coordinates
        PylithScalar slipMagTangent; // scalar magnitude of slipRate in fault plane
        PylithScalar slipMagNormal; // scalar magnitude of slipRate in fault plane
        PylithScalar slipVelMagTangent; // scalar magnitude of slipRate in fault plane
        PylithScalar slipVelMagNormal; // scalar magnitude of slipRate in fault plane
        
        PylithReal cohesion;
        PylithReal opening;
        bool openFreeSurface;
    };

    //typedef void (*frictioncoeffn_type)(const Context&,
    //                                    void*, // rheologyContext
    //                                    PylithReal*); // tractionFriction
    
    typedef void (*frictioncoeffn_type)(PetscReal t,
                                        PetscReal slip,
                                        PetscReal slipVel,
                                        void*, // rheologyContext
                                        PylithReal*); // tractionFriction

    typedef void (*frictiondirfn_type)(const Context&,
                                       const PylithReal, // frictionMag
                                       PylithReal*); // tractionFriction

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:


// --------------------------------------------------------------------------------------------
// create context
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
        
        const PylithInt spaceDim = dim + 1; // :KLUDGE: dim passed in is spaceDim-1
        
        context->dim = dim;
        context->t = t;  
        context->n = &n[0];
        context->dt = constants[0]; 
        context->refDir = &constants[1]; 

        // displacement on + and - sides of fault
        const PylithInt i_disp = 0;
        const PylithInt sOffDispN = sOff[i_disp];
        const PylithInt sOffDispP = sOffDispN+spaceDim;
        context->dispN = &s[sOffDispN];
        context->dispP = &s[sOffDispP];

        // velocity on + and - sides of fault
        context->velN = &s_t[sOffDispN];
        context->velP = &s_t[sOffDispP];
        
        // the lagrange multipliers are the tractions on the fault
        const PylithInt sOffLagrange = sOff[numS-1];
        const PylithScalar* lagrangeGlobal = &s[sOffLagrange];
        context->tractionGlobal = &s[sOffLagrange];
        context->tractionNormalGlobal = &s[sOffLagrange+2];

        // dispN, dispP, velN, vel, and lagrange are all in the global coordinate system
        // convert to fault coordinate system

        // compute slip and slipVel (vectors)
        PetscReal slipGlobal[spaceDim];
        PetscReal slipVelGlobal[spaceDim];
        for (PylithInt i = 0; i < spaceDim; ++i) {
            slipGlobal[i] = context->dispP[i] - context->dispN[i];
            slipVelGlobal[i] = context->velP[i] - context->velN[i];
        } // for

        // rotate from global to fault coordinate system
        switch (spaceDim) {
        case 2: {
            pylith::fekernels::BoundaryDirections::toTN(context->slip, slipGlobal, n);
            context->slipMagTangent = context->slip[0];
            context->slipMagNormal = context->slip[1];

            pylith::fekernels::BoundaryDirections::toTN(context->slipVel, slipVelGlobal, n);
            context->slipVelMagTangent = context->slipVel[0];
            context->slipVelMagNormal = context->slipVel[1];

            PetscReal lagrangeFault[2];
            pylith::fekernels::BoundaryDirections::toTN(lagrangeFault, lagrangeGlobal, n);
            context->tractionShearFault = &lagrangeFault[0];
            context->tractionNormalFault = lagrangeFault[1];
            context->tractionShearMagFault = sqrt(lagrangeFault[0]*lagrangeFault[0]);
        } // case 2
        case 3: {
            // get information for rotating vectors from global to fault coordinates
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            PylithScalar tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);
            
            // compute slip in fault coordinates
            pylith::fekernels::BoundaryDirections::toTN3D(context->slip, slipGlobal, refDir1, refDir2, n);
            context->slipMagTangent = sqrt(context->slip[0]*context->slip[0] + context->slip[1]*context->slip[1]);
            context->slipMagNormal = context->slip[2];

            // compute slip velocity in fault coordinates
            pylith::fekernels::BoundaryDirections::toTN3D(context->slipVel, slipVelGlobal, refDir1, refDir2, n);
            context->slipVelMagTangent = sqrt(context->slipVel[0]*context->slipVel[0] + context->slipVel[1]*context->slipVel[1]);
            context->slipVelMagNormal = context->slipVel[2];
            
            // compute tractions
            PylithScalar lagrangeFault[3];
            pylith::fekernels::BoundaryDirections::toTN3D(lagrangeFault, lagrangeGlobal, refDir1, refDir2, n);
            context->tractionShearFault = &lagrangeFault[0];
            context->tractionNormalFault = lagrangeFault[2];
            context->tractionShearMagFault = fabs(lagrangeFault[0]);
        } // case 3
        default:
            assert(0);
        } // switch

        const PylithInt i_cohesion = numA-3;
        const PylithInt i_opening = numA-2;
        const PylithInt i_openFreeSurface = numA-1;

        assert(a);
        assert(aOff);
        assert(aOff[i_cohesion] >= 0);
        assert(aOff[i_opening] >= 0);
        assert(aOff[i_openFreeSurface] >= 0);

        context->cohesion = a[aOff[i_cohesion]];assert(context->cohesion >= 0.0);
        context->opening = a[aOff[i_opening]];assert(context->opening >= 0.0);
        context->openFreeSurface = a[aOff[i_openFreeSurface]];assert(context->openFreeSurface >= 0.0);
    } // setContext

    // --------------------------------------------------------------------------------------------
    /** Compute friction traction.
     */
    static inline
    void friction(const Context& frictionContext,
                  void* rheologyContext,
                  frictioncoeffn_type frictionCoefFn,
                  frictiondirfn_type frictionDirFn,
                  PylithReal* tractionFriction) {
        assert(rheologyContext);
        assert(frictionCoefFn);
        assert(frictionDirFn);
        assert(tractionFriction);

        PylithReal frictionCoef = 0.0;
        frictionCoefFn(frictionContext.t, frictionContext.slipMagTangent, frictionContext.slipVelMagTangent, rheologyContext, &frictionCoef);

        PylithReal frictionMag;
        if (0.0 == frictionContext.opening) {
            frictionMag = frictionContext.cohesion + frictionCoef * frictionContext.tractionNormalFault;
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
    void frictionDir2D(const Context& context,
                       const PylithReal frictionMag,
                       PylithReal* tractionFriction) {
        assert(tractionFriction);
        assert(context.slipVel);

        const PylithReal slipRate = context.slipVel[0];
        tractionFriction[0] = -slipRate/fabs(slipRate) * frictionMag;
        tractionFriction[1] = 0.0;
    } // frictionDir2D

    // --------------------------------------------------------------------------------------------
    /** Compute friction direction in 3D.
     */
    static inline
    void frictionDir3D(const Context& frictionContext,
                       const PylithReal frictionMag,
                       PylithReal* tractionFriction) {
        assert(tractionFriction);
        assert(frictionContext.slipVel);

        tractionFriction[0] = -frictionContext.slipVel[0]/frictionContext.slipVelMagTangent * frictionMag;
        tractionFriction[1] = -frictionContext.slipVel[1]/frictionContext.slipVelMagTangent * frictionMag;
        tractionFriction[2] = 0.0;
    } // frictionDir3D

    // --------------------------------------------------------------------------------------------
    /** f0 function for elasticity equation: f0u = +\lambda (neg side).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static inline
    void f0u(const Context frictionContext,
             void* rheologyContext,
             frictioncoeffn_type frictionCoefficientFn,
             bool faultSidePos,
             PylithScalar f0[]) {
        
        const PylithInt spaceDim = frictionContext.dim + 1; // :KLUDGE: dim passed in is spaceDim-1

        switch (spaceDim) {
        case 2: {
            frictiondirfn_type frictionDirFn = frictionDir2D;
            PylithReal tractionFriction[2];
            friction(frictionContext,rheologyContext,frictionCoefficientFn,frictionDirFn, tractionFriction);
            for (PylithInt i = 0; i < spaceDim; ++i) {
                PylithReal val = tractionFriction[i] - frictionContext.tractionGlobal[i];
                (faultSidePos) ? f0[i] += -val : f0[i] += val;
            }
        }
        case 3: {
            frictiondirfn_type frictionDirFn = frictionDir3D;
            PylithReal tractionFriction[3];
            friction(frictionContext,rheologyContext,frictionCoefficientFn,frictionDirFn, tractionFriction);
            for (PylithInt i = 0; i < spaceDim; ++i) {
                PylithReal val = tractionFriction[i] - frictionContext.tractionGlobal[i];
                (faultSidePos) ? f0[i] += -val : f0[i] += val;
            }
        }
        default:
            assert(0);
        } // switch
    } // f0u

}; // FaultFriction

#endif // pylith_fekernels_faultfriction_hh

/* End of file */
