
/* Cantera objects */

#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include "Python.h"
#include "Cantera.h"
#include <string>
using namespace std;

static PyObject *ErrorObject;

#include "pyutils.h"

static PyObject *
ct_newCtRxnPath(PyObject *self, PyObject *args)
{
    CtRxnPath *rv;
    PyObject* s;
    if (!PyArg_ParseTuple(args, "O:newCtRxnPath", &s))
        return NULL;
    if (!CtSubstance_check(s)) {
        PyErr_SetString(ErrorObject, "argument must be of type CtSubstance");
    }
    ReactingSubstance* g = s->g;
    rv = newCtRxnPath(g);
    if ( rv == NULL ) return NULL;
    return (PyObject *)rv; 
}


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"Substance", ct_newSubstance,  METH_VARARGS},
    {"CKGas",  ct_newSubstance,  METH_VARARGS},
    {"Substance",  ct_newSubstance,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initcantera) */

    DL_EXPORT(void)
        initcantera(void)
    {
        PyObject *m, *d;
    
        CtSubstance::init();

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
        CtSubstance_Type.ob_type = &PyType_Type;
        CtSubstance_Type.tp_dealloc = (destructor)CtSubstance_dealloc;
        CtSubstance_Type.tp_getattr = (getattrfunc)CtSubstance_getattr;
        CtSubstance_Type.tp_setattr = (setattrfunc)CtSubstance_setattr;

        /* Create the module and add the functions */
        m = Py_InitModule("cantera", ct_methods);
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    } 

}
