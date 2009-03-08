/**
 * @file ctexceptions.h
 *
 *   THis contains 
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2005/10/21 21:35:35 $
 *  $Revision: 1.8 $
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
     * @defgroup errorhandling Error Handling
     *
     *  These classes and related functions are used to handle errors
     *  and unknown events within Cantera.
     * 
     *  The general idea is that exceptions are thrown using the common
     *  base class called CanteraError. Derived types of CanteraError
     *  characterize what type of error is thrown. A list of all
     *  of these errors is kept in the Application class. 
     *
     *  Any exceptions which are not caught cause a fatal error exit
     *  from the program.  
     */
    
    /**
     * Base class for exceptions thrown by Cantera classes.
     *
     * @ingroup errorhandling
     */
    class CanteraError {
    public:
        CanteraError() {}
        CanteraError(string proc, string msg);
        virtual ~CanteraError(){}
    protected:
    };

    /// Array size error.
    /// 
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
