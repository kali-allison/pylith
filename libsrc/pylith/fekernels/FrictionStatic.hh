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

/** @file libsrc/fekernels/FrictionStatic`.hh
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

#if !defined(pylith_fekernels_frictionstatic_hh)
#define pylith_fekernels_frictionstatic_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/FaultFriction.hh" // USES FaultFriction

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::FrictionStatic {
    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal staticCoefficient;
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

        const PylithInt i_static = numA-1;

        assert(a);
        assert(aOff);
        assert(aOff[i_static] >= 0);

        context->staticCoefficient = a[aOff[i_static]];assert(context->staticCoefficient >= 0.0);
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
        *coefficient = context->staticCoefficient;
    } // friction_coefficient

    // --------------------------------------------------------------------------------------------
    // residual fu0 on negative side of fault
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
                const PylithScalar n[],
                const PylithInt numConstants,
                const PylithScalar constants[],
                PylithScalar f0[]) {

        // create FaultFriction context
        pylith::fekernels::FaultFriction::Context frictionContext;
        pylith::fekernels::FaultFriction::setContext(&frictionContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants);

        // create FrictionStatic context
        pylith::fekernels::FrictionStatic::Context rheologyContext;
        pylith::fekernels::FrictionStatic::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants);
        
        // call FaultFriction::f0u(), passing in FrictionStatic context
        // Analogous to IsotropicLinearElasticity::f1v_infinitessimalStrain;
        // instead of calling Elasticity::f1v we will call FaultFriction::f0u
        PylithInt faultSidePos = 0;
        pylith::fekernels::FaultFriction::f0u(frictionContext,&rheologyContext,frictionCoefficient,faultSidePos,f0);
    } // f0u_neg

// --------------------------------------------------------------------------------------------
    // residual fu0 on positive side of fault
    static inline
    void f0u_pos(const PylithInt dim,
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

        // create FrictionStatic context
        pylith::fekernels::FrictionStatic::Context rheologyContext;
        pylith::fekernels::FrictionStatic::setContext(&rheologyContext, dim, numS, numA, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x,
            t, x, n, numConstants, constants);
        
        bool faultSidePos = 1;
        pylith::fekernels::FaultFriction::f0u(frictionContext,&rheologyContext,frictionCoefficient,faultSidePos,f0);
    } // f0u_pos


// --------------------------------------------------------------------------------------------
    // Jacobian Jf0uu on negative side of fault: Jf0uu = 0
    static inline
    void Jf0uu_neg(const PylithInt dim,
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

    } // Jf0uu_neg

// --------------------------------------------------------------------------------------------
    // Jacobian Jf0uu on positive side of fault: Jf0uu = 0
    static inline
    void Jf0uu_pos(const PylithInt dim,
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
    } // Jf0uu_pos


}; // FrictionStatic`

#endif // pylith_fekernels_frictionstatic_hh

/* End of file */
