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
        const PylithScalar* lagrangeMultGlobalCoords; // in global coordinates
        const PylithScalar* dispN; // (spaceDim) displacement negative side of fault
        const PylithScalar* dispP; // (spaceDim) displacement positive side of fault
        const PylithScalar* velN; // (spaceDim) velocity negative side of fault
        const PylithScalar* velP; // (spaceDim) velocity positive side of fault
        
        PylithReal cohesion;
        PylithReal opening;
        bool openFreeSurface;

        Context(void) :
        dim(0),
        t(0.0),
        dt(0.0),
        n(NULL),
        refDir(NULL),
        lagrangeMultGlobalCoords(NULL),
        dispN(NULL),
        dispP(NULL),
        velN(NULL),
        velP(NULL) {}
    };
    
    typedef void (*frictioncoeffn_type)(PetscReal t,
                                        PetscReal slip,
                                        PetscReal slipVel,
                                        void*, // rheologyContext
                                        PylithReal*); // tractionFriction

    typedef void (*frictiondirfn_type)(const Context&,
                                       const PylithReal, // frictionMag
                                       const PylithReal, // slipVelMagTangent
                                       const PylithReal*, // slipVelFaultCoords
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
        context->lagrangeMultGlobalCoords = &s[sOffLagrange];

        const PylithInt i_cohesion = numA-1;

        assert(a);
        assert(aOff);
        assert(aOff[i_cohesion] >= 0);

        context->cohesion = a[aOff[i_cohesion]];assert(context->cohesion >= 0.0);
        context->opening = 0;
        context->openFreeSurface = 1;
        //context->opening = constants[numConstants - 2];assert(context->opening >= 0.0);
        //context->openFreeSurface = constants[numConstants - 1];assert(context->openFreeSurface >= 0.0);
    } // setContext

    // --------------------------------------------------------------------------------------------
    /** Compute slip in fault coords from dispN and dispP in global coords,
     * or compute slipVel in fault coords from velN and velP in global coords.
     *
     * @param[out] diffFaultCoords Values in tagential-normal coordinate system.
     * @param[in] valueXY Values in xy coordinate system.
     * @param[in] normalDir Normal direction unit vector.
     */
    static inline
    void computeDiffFaultCoords(PylithReal diffFaultCoords[],
              const PylithInt spaceDim,
              const PylithReal valuesNGlobalCoords[],
              const PylithReal valuesPGlobalCoords[],
              const PylithReal normalDir[],
              const PylithReal refDir[]) {

        // compute difference between N and P sides of fault
        PetscReal valuesGlobalCoords[spaceDim];
        for (PylithInt i = 0; i < spaceDim; ++i) {
            valuesGlobalCoords[i] = valuesNGlobalCoords[i] - valuesNGlobalCoords[i];
            assert(!isnan(valuesGlobalCoords[i]));
        } // for

        // rotate from global to fault coordinate system
        switch (spaceDim) {
        case 2: {
            pylith::fekernels::BoundaryDirections::toTN(diffFaultCoords, valuesGlobalCoords, normalDir);
            break;
        } // case 2
        case 3: {
            // get information for rotating vectors from global to fault coordinates
            const PylithScalar* refDir1 = &refDir[0];
            const PylithScalar* refDir2 = &refDir[3];
            PylithScalar tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, normalDir);
            pylith::fekernels::BoundaryDirections::toTN3D(diffFaultCoords, valuesGlobalCoords, refDir1, refDir2, normalDir);
            break;
        } // case 3
        default:
            assert(0);
        } // switch

        for (PylithInt i = 0; i < spaceDim; ++i) {
            assert(!isnan(diffFaultCoords[i]));
        } // for

    } // computeDiffFaultCoords



    // --------------------------------------------------------------------------------------------
    /** Compute friction traction.
     */
    static inline
    void friction(const Context& frictionContext,
                  void* rheologyContext,
                  frictioncoeffn_type frictionCoefFn,
                  frictiondirfn_type frictionDirFn,
                  PylithReal* tractionFrictionFaultCoords) {
        assert(rheologyContext);
        assert(frictionCoefFn);
        assert(frictionDirFn);
        assert(tractionFrictionFaultCoords);

        const PylithInt spaceDim = frictionContext.dim + 1; // :KLUDGE: dim passed in is spaceDim-1
        if (2 != spaceDim && 3 != spaceDim ) { assert(0); }

        // compute slip and slipVel in fault coords
        PylithReal slipFaultCoords[2];
        computeDiffFaultCoords(slipFaultCoords,spaceDim, frictionContext.dispN, frictionContext.dispP,frictionContext.n,frictionContext.refDir);
        PylithReal slipMagTangent = (2 == spaceDim) ? fabs(slipFaultCoords[0]) : sqrt(pow(slipFaultCoords[0],2) + pow(slipFaultCoords[1],2));
        PylithReal slipVelFaultCoords[2];
        computeDiffFaultCoords(slipVelFaultCoords,spaceDim, frictionContext.dispN, frictionContext.dispP,frictionContext.n,frictionContext.refDir);
        PylithReal slipVelMagTangent = (2 == spaceDim) ? abs(slipVelFaultCoords[0]) : sqrt(pow(slipVelFaultCoords[0],2) + pow(slipVelFaultCoords[1],2));

        // compute coefficient of friction
        PylithReal frictionCoef = 0.0;
        frictionCoefFn(frictionContext.t, slipMagTangent, slipVelMagTangent, rheologyContext, &frictionCoef);

        // compute magnitude of friction
        PylithReal lagrangeMultFaultCoords[2];
        pylith::fekernels::BoundaryDirections::toTN(lagrangeMultFaultCoords,frictionContext.lagrangeMultGlobalCoords,frictionContext.n);
        PylithReal tractionNormalFaultCoords = (2 == spaceDim) ? -lagrangeMultFaultCoords[1] : -lagrangeMultFaultCoords[2];
        PylithReal frictionMag;
        if (0.0 == frictionContext.opening) {
            frictionMag = frictionContext.cohesion - frictionCoef * tractionNormalFaultCoords;
        } else if (!frictionContext.openFreeSurface) {
            frictionMag = frictionContext.cohesion;
        } else {
            frictionMag = 0.0;
        } // if/else

        frictionDirFn(frictionContext, frictionMag, slipVelMagTangent, slipVelFaultCoords, tractionFrictionFaultCoords);
        for (PylithInt i = 0; i < spaceDim; ++i) {
            assert(!isnan(tractionFrictionFaultCoords[i]));
        }
    } // friction

    // --------------------------------------------------------------------------------------------
    /** Compute friction direction in 2D.
     */
    static inline
    void frictionDir2D(const Context& context,
                       const PylithReal frictionMag,
                       const PylithReal slipVelMagTangent,
                       const PylithReal* slipVelFaultCoords,
                       PylithReal* tractionFrictionFaultCoords) {
        assert(tractionFrictionFaultCoords);

        if (slipVelMagTangent == 0) {
            PylithReal lagrangeMultFaultCoords[2];
            pylith::fekernels::BoundaryDirections::toTN(lagrangeMultFaultCoords,context.lagrangeMultGlobalCoords,context.n);
            PylithReal a = frictionMag - lagrangeMultFaultCoords[0];
            tractionFrictionFaultCoords[0] = a/sqrt(a*a) * frictionMag;
        }
        else {
            tractionFrictionFaultCoords[0] = -slipVelFaultCoords[0]/fabs(slipVelMagTangent) * frictionMag;
        }
        tractionFrictionFaultCoords[1] = 0.0;
    } // frictionDir2D

    // --------------------------------------------------------------------------------------------
    /** Compute friction direction in 3D.
     */
    static inline
    void frictionDir3D(const Context& context,
                       const PylithReal frictionMag,
                       const PylithReal slipVelMagTangent,
                       const PylithReal* slipVelFaultCoords,
                       PylithReal* tractionFrictionFaultCoords) {
        assert(tractionFrictionFaultCoords);

        if (slipVelMagTangent == 0) {
            PylithReal lagrangeMultFaultCoords[3];
            const PylithScalar* refDir1 = &context.refDir[0];
            const PylithScalar* refDir2 = &context.refDir[3];
            pylith::fekernels::BoundaryDirections::toTN3D(lagrangeMultFaultCoords,context.lagrangeMultGlobalCoords,refDir1,refDir2,context.n);
            PylithReal a = frictionMag - lagrangeMultFaultCoords[0];
            PylithReal b = frictionMag - lagrangeMultFaultCoords[1];
            PylithReal temp = sqrt(a*a + b*b);
            tractionFrictionFaultCoords[0] = a/temp * frictionMag;
            tractionFrictionFaultCoords[1] = b/temp * frictionMag;
        }
        else {
            tractionFrictionFaultCoords[0] = -slipVelFaultCoords[0]/slipVelMagTangent * frictionMag;
            tractionFrictionFaultCoords[1] = -slipVelFaultCoords[1]/slipVelMagTangent * frictionMag;
        }
        tractionFrictionFaultCoords[2] = 0.0;
    } // frictionDir3D

    // --------------------------------------------------------------------------------------------
    /** f0 function for elasticity equation: f0u = -(traction_friction - lambda).
     *
     * Solution fields: [disp(dim), ..., lagrange(dim)]
     */
    static inline
    void f0u(const Context frictionContext,
             void* rheologyContext,
             frictioncoeffn_type frictionCoefficientFn,
             PylithScalar f0[]) {
        
        const PylithInt spaceDim = frictionContext.dim + 1; // :KLUDGE: dim passed in is spaceDim-1
        const PylithInt offN = 0;
        const PylithInt offP = offN + spaceDim;

        switch (spaceDim) {
        case 2: {
            frictiondirfn_type frictionDirFn = frictionDir2D;
            PylithReal tractionFrictionFaultCoords[2];
            friction(frictionContext,rheologyContext,frictionCoefficientFn,frictionDirFn, tractionFrictionFaultCoords);
            PylithReal tractionFrictionGlobalCoords[2];
            pylith::fekernels::BoundaryDirections::toXY(tractionFrictionGlobalCoords,tractionFrictionFaultCoords,frictionContext.n);
            for (PylithInt i = 0; i < spaceDim; ++i) {
                PylithReal val = tractionFrictionGlobalCoords[i] - frictionContext.lagrangeMultGlobalCoords[i];
                assert(!isnan(val));
                f0[offN + i] += -val; // was - val
                f0[offP + i] += -val; // was + val
            }
            break;
        }
        case 3: {
            frictiondirfn_type frictionDirFn = frictionDir3D;
            PylithReal tractionFrictionFaultCoords[3];
            friction(frictionContext,rheologyContext,frictionCoefficientFn,frictionDirFn, tractionFrictionFaultCoords);
            PylithReal tractionFrictionGlobalCoords[3];
            const PylithScalar* refDir1 = &frictionContext.refDir[0];
            const PylithScalar* refDir2 = &frictionContext.refDir[3];
            pylith::fekernels::BoundaryDirections::toXYZ(tractionFrictionGlobalCoords,tractionFrictionFaultCoords,refDir1,refDir2,frictionContext.n);
            for (PylithInt i = 0; i < spaceDim; ++i) {
                PylithReal val = tractionFrictionGlobalCoords[i] - frictionContext.lagrangeMultGlobalCoords[i];
                assert(!isnan(val));
                f0[offN + i] += -val;
                f0[offP + i] += +val;
            }
            break;
        }
        default:
            assert(0);
        } // switch
    } // f0u

}; // FaultFriction

#endif // pylith_fekernels_faultfriction_hh

/* End of file */
