/**
 * @file global.h
 *
 * These functions handle various utility functions, and store some parameters in 
 * global storage that are accessible at all times. 
8
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

    class XML_Node;

    /// Number of errors that have been encountered so far
    int nErrors();

    /// The last error message
    string lastErrorMessage();

    /// Add an error message
    void setError(string r, string msg);

    /// Print the error messages to stream f
    void showErrors(ostream& f);

    /// Discard the last error message
    void popError();

    /// Find a file on a search path
    string findInputFile(string name);

    /// Add a directory to the search path
    void addDirectory(string dir);

    /// Delete the application singly defined information
    void appdelete();

    // Write a message (deprecated; use writelog)
    //void write(const string& msg);

    //  Write a message (deprecated; use writelog)
    //void write(const char* msg);

    /// The root directory where Cantera is installed
    string canteraRoot();

    /// Set the temporary file directory. Default: /tmp.
    void setTmpDir(string tmp);

    /// The directory where temporary files may be created
    string tmpDir();


    /**
     * Write a diagnostic message to standard output.
     *
     * There are several versions of function writelog, each designed
     * for a particular environment.  One version is designed for use
     * in C++ programs, or in any environment where messages can be
     * written to the standard output stream. Other versions are
     * written specifically for Matlab and for Python, which call
     * Matlab or Python functions, respectively, to display the
     * message.  This is particularly important for Matlab, since
     * everything written to the standard output simply
     * disappears. For Python, the messages show up, but are not in
     * the right order if Python scripts also print output. Hence the
     * need for separate versions.
     *
     * The C++ version may be linked to an application by linking in
     * the library libctcxx.a (-lctcxx). The Python and Matlab
     * interfaces do not link this library, and instead link to their
     * own versions of writelog.
     */ 

    void writelog(const string& msg);
    void writelog(const char* msg);

    
    //void getlog(string& s);
    //void clearlog();

    /**
     * Return the conversion factor to convert unit string 'unit' to SI units.
     */ 
    doublereal toSI(string unit);
    doublereal actEnergyToSI(string unit);
    //

    XML_Node* get_XML_File(string file);
    void close_XML_File(string file);
}

#endif
