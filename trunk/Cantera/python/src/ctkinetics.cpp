/**
 * @file ctkinetics.cpp
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
#include <string>
#include <vector>
#include <iostream>
using namespace std;


//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h" 

static PyObject*
kin_newFromXML(PyObject *self, PyObject *args) {
    int mxml, iphase, neighbor1, neighbor2;
    if (!PyArg_ParseTuple(args, "iiii:newFromXML", &mxml, 
            &iphase, &neighbor1, &neighbor2)) 
        return NULL;
    int n = newKineticsFromXML(mxml, iphase, neighbor1, neighbor2);
    if (n < 0) return reportError(n);
    return Py_BuildValue("i",n);
}

static PyObject*
kin_delete(PyObject *self, PyObject *args)
{
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_delete", &kin)) return NULL;    
    delKinetics(kin);
    return Py_BuildValue("i",0);
}

static PyObject*
kin_nspecies(PyObject *self, PyObject *args) {
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_nspecies", &kin)) return NULL;
    return Py_BuildValue("i",kin_nSpecies(kin));        
}

static PyObject*
kin_rstoichcoeff(PyObject *self, PyObject *args) {
    int kin, i, k;
    if (!PyArg_ParseTuple(args, "iii:kin_rstoichcoeff", &kin, &k, &i)) 
        return NULL;
    return Py_BuildValue("d",kin_reactantStoichCoeff(kin, k, i));        
}

static PyObject*
kin_pstoichcoeff(PyObject *self, PyObject *args) {
    int kin, i, k;
    if (!PyArg_ParseTuple(args, "iii:kin_pstoichcoeff", &kin, &k, &i)) 
        return NULL;
    return Py_BuildValue("d",kin_productStoichCoeff(kin, k, i));        
}

static PyObject*
kin_nrxns(PyObject *self, PyObject *args) {
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_nreactions", &kin)) return NULL;
    return Py_BuildValue("i",kin_nReactions(kin));        
}

static PyObject*
kin_isrev(PyObject *self, PyObject *args) {
    int kin, i;
    if (!PyArg_ParseTuple(args, "ii:kin_isrev", &kin, &i)) return NULL;
    return Py_BuildValue("i",kin_isReversible(kin,i));        
}

static PyObject*
kin_rxntype(PyObject *self, PyObject *args) {
    int kin, i;
    if (!PyArg_ParseTuple(args, "ii:kin_rxntype", &kin, &i)) return NULL;
    return Py_BuildValue("i",kin_reactionType(kin,i));        
}

static PyObject*
kin_multiplier(PyObject *self, PyObject *args) {
    int kin, i;
    if (!PyArg_ParseTuple(args, "ii:kin_multiplier", &kin, &i)) return NULL;
    return Py_BuildValue("d",kin_multiplier(kin,i));        
}

static PyObject*
kin_setMultiplier(PyObject *self, PyObject *args) {
    int kin, i;
    double v;
    if (!PyArg_ParseTuple(args, "iid:kin_setMultiplier", &kin, &i, &v)) return NULL;
    return Py_BuildValue("i",kin_setMultiplier(kin,i,v));        
}


static PyObject*
kin_type(PyObject *self, PyObject *args) {
    int kin, i;
    if (!PyArg_ParseTuple(args, "i:kin_type", &kin)) return NULL;
    return Py_BuildValue("i",kin_type(kin));        
}

static PyObject*
kin_start(PyObject *self, PyObject *args) {
    int kin, p;
    if (!PyArg_ParseTuple(args, "ii:kin_start", &kin, &p)) return NULL;
    return Py_BuildValue("i",kin_start(kin,p));        
}

static PyObject*
kin_speciesIndex(PyObject *self, PyObject *args) {
    int kin;
    char *nm, *ph;
    if (!PyArg_ParseTuple(args, "iss:kin_speciesIndex", &kin, &nm, &ph)) 
        return NULL;
    return Py_BuildValue("i",kin_speciesIndex(kin,nm,ph));        
}

static PyObject*
kin_getarray(PyObject *self, PyObject *args)
{
    int kin;
    int job;
    
    if (!PyArg_ParseTuple(args, "ii:kin_getarray", &kin, &job)) 
        return NULL;

    // array attributes
    int iok = -22;
    int nrxns = kin_nReactions(kin);
    int nsp = kin_nSpecies(kin);
    int ix;
    if (job < 45) ix = nrxns; else ix = nsp;
 
    PyArrayObject* x = 
        (PyArrayObject*)PyArray_FromDims(1, &ix, PyArray_DOUBLE);
    double* xd = (double*)x->data;

    switch (job) {
    case 10:
        iok = kin_getFwdRatesOfProgress(kin, nrxns, xd);
        break;
    case 20:
        iok = kin_getRevRatesOfProgress(kin, nrxns, xd);
        break;
    case 30:
        iok = kin_getNetRatesOfProgress(kin, nrxns, xd);
        break;
    case 40:
        iok = kin_getEquilibriumConstants(kin, nrxns, xd);
        break;
    case 50:
        iok = kin_getCreationRates(kin, nsp, xd);
        break;
    case 60:
        iok = kin_getDestructionRates(kin, nsp, xd);
        break;
    case 70:
        iok = kin_getNetProductionRates(kin, nsp, xd);
        break;
    case 80:
        iok = kin_getSourceTerms(kin, nsp, xd);
        break;

        break;
    default:
        ;
    }
    if (iok >= 0) {
        return PyArray_Return(x);
    }
    else 
        return reportError(iok);
}



// string attributes 
static PyObject*
kin_getstring(PyObject *self, PyObject *args)
{
    int kin, job, i, iok = -3;
    int buflen;
    char* output_buf = 0;
    if (!PyArg_ParseTuple(args, "iii:kin_getstring", &kin, &job, &i)) 
        return NULL;
    switch (job) {
    case 1:
        buflen = 80;
        output_buf = new char[buflen];
        iok = kin_getReactionString(kin, i, buflen, output_buf);
        break;
    default:
        iok = -10;
    }
    if (iok >= 0) {
        PyObject* str = Py_BuildValue("s",output_buf);
        delete output_buf;
        return str;
    }
    delete output_buf;
    if (iok == -1) 
        return reportCanteraError();
    else {
        PyErr_SetString(ErrorObject,"Unknown string attribute");
        return NULL;
    }        
}


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"KineticsFromXML", kin_newFromXML,  METH_VARARGS},
    {"delete", kin_delete,  METH_VARARGS},
    {"nspecies", kin_nspecies,  METH_VARARGS},
    {"nreactions", kin_nrxns,  METH_VARARGS},
    {"isreversible", kin_isrev,  METH_VARARGS},
    {"rstoichcoeff", kin_rstoichcoeff,  METH_VARARGS},
    {"pstoichcoeff", kin_pstoichcoeff,  METH_VARARGS},
    {"rxntype", kin_rxntype,  METH_VARARGS},
    {"type", kin_type,  METH_VARARGS},
    {"start", kin_start,  METH_VARARGS},
    {"multiplier", kin_multiplier,  METH_VARARGS},
    {"setMultiplier", kin_setMultiplier,  METH_VARARGS},
    {"speciesIndex", kin_speciesIndex,  METH_VARARGS},
    {"getarray", kin_getarray,  METH_VARARGS},
    {"getstring", kin_getstring,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {
 
    /* Initialization function for the module (*must* be called initctkinetics) */

    DL_EXPORT(void) initctkinetics(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctkinetics", ct_methods);
        import_array();

        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}
