/**
 *  @file cxx/src/writelog.cpp
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.5 $
 * $Date: 2005/01/07 10:26:43 $
 */

#include <string>
#include <iostream>
using namespace std;

namespace Cantera {
    /**
     *
     * writelog():
     *
     * Write a diagnostic message to standard output.
     * -> defined in global.h
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
    void writelog(const string& s) {
        cout << s;
    }

    /**
     * Write an error message and quit.
     */
    void error(const string& msg) {
        cerr << msg << endl;
        exit(-1);
    }

    int userInterface() {
        return 0;
    }

}
