/**
 * @file global.h
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_GLOBAL_H
#define CT_GLOBAL_H

#include <string>
#include "ct_defs.h"

namespace Cantera {

    /// Number of errors that have been encountered
    int nErrors();

    /// The last error message
    string lastErrorMessage();

    /// Add an error message
    void setError(string r, string msg);

    /// Print the error messages to f
    void showErrors(ostream& f);

    /// Discard the last error message
    void popError();

    /// Find a file on a search path
    string findInputFile(string name);

    /// Add a directory to the search path
    void addDirectory(string dir);

    /// Write a message
    void write(const string& msg);

    ///  Write a message
    void write(const char* msg);

    string canteraRoot();

    void setTmpDir(string tmp);
    string tmpDir();

    /**
     * Write a diagnostic message to an internal buffer.
     */
    void writelog(const string& msg);
    void writelog(const char* msg);

    
    void getlog(string& s);
    void clearlog();
    doublereal toSI(string unit);
    doublereal actEnergyToSI(string unit);
    //
}

#endif
