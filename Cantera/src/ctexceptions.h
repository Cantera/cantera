/**
 * @file ctexceptions.h
 */

// $Author$
// $Revision$
// $Date$

// Copyright 2001  California Institute of Technology

#ifndef CT_CTEXCEPTIONS_H
#define CT_CTEXCEPTIONS_H

//#include "global.h"
//#include "stringUtils.h"

#include <string>
using namespace std;

namespace Cantera {

    /**
     * Base class for exceptions thrown by Cantera classes
     */
    class CanteraError {
    public:
        CanteraError() {}
        CanteraError(string proc, string msg);
        virtual ~CanteraError(){}
    protected:
    };

    class ArraySizeError : public CanteraError {
    public:
        ArraySizeError(string proc, int sz, int reqd);
    };

    class ElementRangeError : public CanteraError {
    public:
        ElementRangeError(string func, int m, int mmax);
    };


}

#endif
