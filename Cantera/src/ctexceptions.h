/**
 * @file ctexceptions.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_CTEXCEPTIONS_H
#define CT_CTEXCEPTIONS_H

#include <string>
using namespace std;

// See file misc.cpp for implementations of methods/functions declared
// here.

namespace Cantera {

    /**
     * Base class for exceptions thrown by Cantera classes.
     */
    class CanteraError {
    public:
        CanteraError() {}
        CanteraError(string proc, string msg);
        virtual ~CanteraError(){}
    protected:
    };

    /// Array size error.
    class ArraySizeError : public CanteraError {
    public:
        ArraySizeError(string proc, int sz, int reqd);
    };

    /// Exception thrown if an element index is out of range.
    class ElementRangeError : public CanteraError {
    public:
        ElementRangeError(string func, int m, int mmax);
    };

    void deprecatedMethod(string classnm, string oldnm, string newnm);
    void removeAtVersion(string func, string version);
}

#endif
