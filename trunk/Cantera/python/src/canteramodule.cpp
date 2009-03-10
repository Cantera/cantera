/**
 * @file canteramodule.cpp
 * Cantera Python Interface
 *
 */

//  copyright 2001  David G. Goodwin


// turn off warnings about long names under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Python.h"
#include "Cantera.h"
#include <string>
using namespace std;

//  constants defined in the module
static PyObject *ErrorObject;
static PyObject *OneAtmos;
static PyObject *GasCon;


local includes
#include "pyutils.h" 
#include "CtSubstance.h"
#include "CtRxnPath.h"
#include "CtReactor.h"


/* Function of no arguments returning new CtSubstance object */

static PyObject *
ct_newSubstance(PyObject *self, PyObject *args)
{
    CtSubstance *rv=0;
    char* subtype;
    PyObject* params;
    if (!PyArg_ParseTuple(args, "sO:newSubstance", &subtype, &params)) 
        return NULL;
    if (!PySequence_Check(params)) {
        PyErr_SetString(ErrorObject, 
            "usage: Substance('<substance_type>', [<parameters>])");
        return NULL;
    }
    string stype(subtype);
    if (stype == "chemkin") {
        rv = newChemkinSubstance(params);        
    }
    if ( rv == NULL ) return NULL;
    return (PyObject *)rv;  
}


static PyObject *
ct_newRxnPath(PyObject *self, PyObject *args)
{
    CtRxnPath *rv;
    if (!PyArg_ParseTuple(args, ":newCtRxnPath"))
        return NULL;
    rv = newCtRxnPath();
    if ( rv == NULL ) return NULL;
    return (PyObject *)rv; 
}

static PyObject *
ct_newReactor(PyObject *self, PyObject *args)
{
    CtReactor *rv;
    char* s;
    if (!PyArg_ParseTuple(args, "s:newCtReactor", &s))
        return NULL;
    rv = newCtReactor(s);
    if ( rv == NULL ) return NULL;
    return (PyObject *)rv; 
}



/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"Substance", ct_newSubstance,  METH_VARARGS},
    {"ReactionPath", ct_newRxnPath, METH_VARARGS},
    {"Reactor", ct_newReactor, METH_VARARGS},
    {"ChemEquil", ct_newChemEquil, METH_VARARGS},
    //{"Integrator", ct_newIntegrator, METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initcantera) */

    DL_EXPORT(void)
        initct(void)
    {
        PyObject *m, *d;
    
        CtSubstance::init();
        CtRxnPath::init();
        CtReactor::init();

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */

        CtSubstance_Type.ob_type = &PyType_Type;
        CtSubstance_Type.tp_dealloc = (destructor)CtSubstance_dealloc;
        CtSubstance_Type.tp_getattr = (getattrfunc)CtSubstance_getattr;
        CtSubstance_Type.tp_setattr = (setattrfunc)CtSubstance_setattr;
        CtSubstance_Type.tp_print = (printfunc)CtSubstance_print;

        CtRxnPath_Type.ob_type = &PyType_Type;
        CtRxnPath_Type.tp_dealloc = (destructor)CtRxnPath_dealloc;
        CtRxnPath_Type.tp_getattr = (getattrfunc)CtRxnPath_getattr;
        CtRxnPath_Type.tp_setattr = (setattrfunc)CtRxnPath_setattr;

        CtReactor_Type.ob_type = &PyType_Type;
        CtReactor_Type.tp_dealloc = (destructor)CtReactor_dealloc;
        CtReactor_Type.tp_getattr = (getattrfunc)CtReactor_getattr;
        CtReactor_Type.tp_setattr = (setattrfunc)CtReactor_setattr;
 
        /* Create the module and add the functions */
        m = Py_InitModule("cantera", ct_methods);
 
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
