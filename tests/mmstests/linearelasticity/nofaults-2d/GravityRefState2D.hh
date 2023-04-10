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

#include <portinfo>

#include "TestLinearElasticity.hh" // USES TestLinearElasticity_Data

namespace pylith {
    class GravityRefState2D;
}

class pylith::GravityRefState2D {
public:

    // Data factory methods

    static TestLinearElasticity_Data* TriP1(void);

    static TestLinearElasticity_Data* TriP2(void);

    static TestLinearElasticity_Data* TriP3(void);

    static TestLinearElasticity_Data* QuadQ1(void);

    static TestLinearElasticity_Data* QuadQ2(void);

    static TestLinearElasticity_Data* QuadQ3(void);

private:

    GravityRefState2D(void); ///< Not implemented
}; // GravityRefState2D

// End of file
