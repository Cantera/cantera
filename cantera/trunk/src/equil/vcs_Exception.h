/**
 * @file vcs_Exception.h
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef VCS_EXCEPTION_H
#define VCS_EXCEPTION_H

#include <string>

namespace VCSnonideal
{

class vcsError
{
public:
    vcsError(std::string proc, std::string msg, int errorCode=-1);
    virtual ~vcsError() {}
protected:
    std::string m_proc;
    std::string m_msg;
    int  m_errorCode;
};


//!  Assertion must be true or an error is thrown
/*!
 * Assertion must be true or else a vcsError is thrown. A diagnostic
 * string indicating where the error
 * occured is added to the thrown object.
 *
 * @param expr  Boolean expression that must be true
 * @param proc  Character string or std:string expression indicating the procedure
 *              where the assertion failed
 * @ingroup errorhandling
 */
#define AssertThrowVCS(expr, proc)  ((expr) ? (void) 0 : throw vcsError(proc, std::string("failed Assert: ") + #expr,-1))

#ifdef DEBUG_HKM
#define DebugAssertThrowVCS(expr, proc)  ((expr) ? (void) 0 : throw vcsError(proc, std::string("failed debugAssert: ") + #expr,-1))
#else
#define DebugAssertThrowVCS(expr, proc)
#endif

}

#endif
