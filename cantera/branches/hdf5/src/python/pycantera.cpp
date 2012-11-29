/**
 * @file pycantera.cpp
 *
 * This is the main file for the Python interface module.
 *
 */

#include "cantera/base/config.h"

#include "Python.h"

#ifdef HAS_NUMERIC
#include "Numeric/arrayobject.h"
#endif

#ifdef HAS_NUMARRAY
#include "numarray/arrayobject.h"
#endif

#ifdef HAS_NUMPY
#include "numpy/arrayobject.h"
#endif

#include "clib/ct.h"
#include "clib/ctxml.h"
#include "clib/ctsurf.h"
#include "clib/ctrpath.h"
#include "clib/ctreactor.h"
#include "clib/ctfunc.h"
#include "clib/ctonedim.h"
#include "clib/ctmultiphase.h"

#include <iostream>
using namespace std;

//  constants defined in the module
static PyObject* ErrorObject;

// local includes
#include "pyutils.h"

#include "ctphase_methods.cpp"
#include "ctthermo_methods.cpp"
#include "ctkinetics_methods.cpp"
#include "cttransport_methods.cpp"
#include "ctxml_methods.cpp"
#include "ctfuncs.cpp"
#include "ctsurf_methods.cpp"
#include "ctrpath_methods.cpp"
#include "ctreactor_methods.cpp"
#include "ctfunc_methods.cpp"
#include "ctonedim_methods.cpp"
#include "ctmultiphase_methods.cpp"

static PyObject*
pyct_appdelete(PyObject* self, PyObject* args)
{
    return Py_BuildValue("i",ct_appdelete());
}

#include "methods.h"
#include "pylogger.h"

extern "C" {

    /* Initialization function for the module */

    DL_EXPORT(void) init_cantera(void)
    {
        PyObject* m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */

        /* Create the module and add the functions */
        m = Py_InitModule("_cantera", ct_methods);
        import_array();
        Cantera::Logger* pylog = new Cantera::Py_Logger;
        setLogWriter(pylog);

        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException((char*)"cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
#ifdef HAS_NUMERIC
        PyDict_SetItemString(d, "nummod",PyString_FromString("Numeric"));
#else
#ifdef HAS_NUMARRAY
        PyDict_SetItemString(d, "nummod",PyString_FromString("numarray"));
#else
        PyDict_SetItemString(d, "nummod",PyString_FromString("numpy"));
#endif
#endif
    }

}

