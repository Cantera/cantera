/**
 * @file cttransport.cpp
 * Cantera Python Interface
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
//#include <vector>
#include <iostream>
using namespace std;

//#include "Cantera.h"

//  constants defined in the module
static PyObject *ErrorObject;

// local includes  
#include "pyutils.h" 



/**
 * Create a new Transport object.
 */
static PyObject *
py_transport_new(PyObject *self, PyObject *args) {
    char* model;
    int ph;
    int loglevel;
    if (!PyArg_ParseTuple(args, "sii:transport_new", &model, 
            &ph, &loglevel)) 
        return NULL;
    int n = newTransport(model, ph, loglevel);
    if (n < 0) return reportError(n);
    return Py_BuildValue("i",n);
}


/**
 * Delete the Phase object.
 */
static PyObject*
py_transport_delete(PyObject *self, PyObject *args)
{
    int tr;
    if (!PyArg_ParseTuple(args, "i:transport_delete", &tr)) 
        return NULL;    
    delTransport(tr);
    return Py_BuildValue("i",0);
}


static PyObject*
py_viscosity(PyObject *self, PyObject *args) {
    int n;
    if (!PyArg_ParseTuple(args, "i:py_viscosity", &n)) return NULL;
    double mu = trans_viscosity(n);        
    if (mu < 0.0) return reportError(int(mu));
    return Py_BuildValue("d",mu);        
}

static PyObject*
py_thermalConductivity(PyObject *self, PyObject *args) {
    int n;
    if (!PyArg_ParseTuple(args, "i:py_thermalConductivity", &n)) return NULL;
    double lambda = trans_thermalConductivity(n);
    if (lambda < 0.0) return reportError(int(lambda));
    return Py_BuildValue("d",lambda);
}

static PyObject*
py_thermalDiffCoeffs(PyObject *self, PyObject *args) {
    int n, idt;
    if (!PyArg_ParseTuple(args, "ii:py_thermalDiffCoeffs", &n, &idt)) 
        return NULL;
    PyArrayObject* dt = 
        (PyArrayObject*)PyArray_FromDims(1, &idt, PyArray_DOUBLE);
    int iok = trans_getThermalDiffCoeffs(n, idt, (double*)dt->data);
    if (iok < 0) return reportError(iok);
    return PyArray_Return(dt);
}

static PyObject*
py_binaryDiffCoeffs(PyObject *self, PyObject *args) {
    int n, id;
    if (!PyArg_ParseTuple(args, "ii:py_binaryDiffCoeffs", &n, &id)) 
        return NULL;
    int idim[2];
    idim[0] = id;
    idim[1] = id;
    PyArrayObject* d = 
        (PyArrayObject*)PyArray_FromDims(2, idim, PyArray_DOUBLE);
    int iok = trans_getBinDiffCoeffs(n, id, (double*)d->data);
    if (iok < 0) return reportError(iok);
    return PyArray_Return(d);
}

static PyObject*
py_mixDiffCoeffs(PyObject *self, PyObject *args) {
    int n, id;
    if (!PyArg_ParseTuple(args, "ii:py_mixDiffCoeffs", &n, &id)) 
        return NULL;
    PyArrayObject* d = 
        (PyArrayObject*)PyArray_FromDims(1, &id, PyArray_DOUBLE);
    int iok = trans_getMixDiffCoeffs(n, id, (double*)d->data);
    if (iok < 0) return reportError(iok);
    return PyArray_Return(d);
}

static PyObject*
py_multiDiffCoeffs(PyObject *self, PyObject *args) {
    int n, id;
    if (!PyArg_ParseTuple(args, "ii:py_multiDiffCoeffs", &n, &id)) 
        return NULL;
    //vector_int idim(2,id);
    int idim[2];
    idim[0] = id;
    idim[1] = id;
    PyArrayObject* d = 
        (PyArrayObject*)PyArray_FromDims(2, idim, PyArray_DOUBLE);
    int iok = trans_getMultiDiffCoeffs(n, id, (double*)d->data);
    if (iok < 0) return reportError(iok);
    return PyArray_Return(d);
}


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"Transport", py_transport_new,  METH_VARARGS},
    {"delete", py_transport_delete,  METH_VARARGS},
    {"viscosity", py_viscosity,  METH_VARARGS},
    {"thermalConductivity", py_thermalConductivity,  METH_VARARGS},
    {"thermalDiffCoeffs", py_thermalDiffCoeffs,  METH_VARARGS},
    {"binaryDiffCoeffs", py_binaryDiffCoeffs,  METH_VARARGS},
    {"mixDiffCoeffs", py_mixDiffCoeffs,  METH_VARARGS},
    {"multiDiffCoeffs", py_multiDiffCoeffs,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initcttransport) */

    DL_EXPORT(void) initcttransport(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("cttransport", ct_methods);
        import_array();

        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }
}
