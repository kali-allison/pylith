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

#include "TestFault.hh" // USES TestFault_Data

namespace pylith {
    class ThreeBlocksStatic;
}

class pylith::ThreeBlocksStatic {
public:

    // Data factory methods
    static TestFault_Data* TriP1(void);

    static TestFault_Data* TriP2(void);

    static TestFault_Data* TriP3(void);

    static TestFault_Data* TriP4(void);

    static TestFault_Data* QuadQ1(void);

    static TestFault_Data* QuadQ1Distorted(void);

    static TestFault_Data* QuadQ2(void);

    static TestFault_Data* QuadQ3(void);

    static TestFault_Data* QuadQ4(void);

private:

    ThreeBlocksStatic(void); ///< Not implemented
}; // ThreeBlocksStatic

// End of file
