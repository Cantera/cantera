/**
 * @file pycantera.cpp
 *
 * This is the main file for the Python interface module.
 *
 */

// turn off warnings about long names under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


#include "Python.h"
#include "Numeric/arrayobject.h"

#include "ct.h"
#include "ctxml.h"
#include "ctstagn.h" 
#include "ctsurf.h"
#include "ctbdry.h"
#include "ctrpath.h"
#include "ctreactor.h" 
#include "ctfunc.h"
#include "ctonedim.h"

#include <iostream>
using namespace std;

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h"

#include "ctphase_methods.cpp"
#include "ctthermo_methods.cpp"
#include "ctkinetics_methods.cpp"
#include "cttransport_methods.cpp"
#include "ctxml_methods.cpp"
#include "ctflow_methods.cpp"
#include "ctfuncs.cpp"
#include "ctsurf_methods.cpp"
#include "ctbndry_methods.cpp"
#include "ctrpath_methods.cpp"
#include "ctreactor_methods.cpp"
#include "ctfunc_methods.cpp"
#include "ctonedim_methods.cpp"

#include "methods.h"

extern "C" {

    /* Initialization function for the module */

    DL_EXPORT(void) init_cantera(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("_cantera", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}

