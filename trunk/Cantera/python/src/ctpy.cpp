/**
 * @file ctpy.cpp
 * Cantera Python Interface
 *
 */

// turn off warnings about long names under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Python.h"
#include "ct.h" 
#include <string>
using namespace std;

//  constants defined in the module
static PyObject *ErrorObject;
static PyObject *OneAtmos;
static PyObject *GasCon;

// local includes
#include "pyutils.h" 


static reportCanteraError() {
    PyErr_SetString(ErrorObject, 
        Cantera::Application::errorMessage.back().c_str());
        return NULL;
}


static PyObject *
ct_newPhase(PyObject *self, PyObject *args)
{
    CtPhase *rv=0;
    int job;
    if (!PyArg_ParseTuple(args, "i:ct_newPhase", &job)) 
        return NULL;
    int n = newPhase(job);
    if (n < 0) {
        PyErr_SetString(ErrorObject,"Unknown species thermo manager");
        return NULL;
    }
    return Py_BuildValue("i",n);
}


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"Phase", ct_newPhase,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initct) */

    DL_EXPORT(void) initct(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ct", ct_methods);
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);

        // one atmosphere
        OneAtmos = PyFloat_FromDouble(OneAtm);
        PyDict_SetItemString(d, "OneAtm", OneAtmos);

        // gas constant        
        GasCon = PyFloat_FromDouble(GasConstant);
        PyDict_SetItemString(d, "GasConstant", GasCon);        
    }
}
