
/**
 * @file ctpyfunc1.cpp
 *
 */

// turn off warnings about long names under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
using namespace std;

#include "Python.h"
#include "Numeric/arrayobject.h"

#include "ct.h"
#include "ctfunc.h" 

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h"


static PyObject*
py_func_new(PyObject *self, PyObject *args)
{
    int type, n;
    PyObject* c;
    if (!PyArg_ParseTuple(args, "iiO:func_new", &type, &n, &c))
        return NULL;
    PyArrayObject* coeffs = (PyArrayObject*)c;
    double* xd = (double*)coeffs->data;
    int lenc = coeffs->dimensions[0];
    int nn = func_new(type, n, lenc, xd);
    if (nn < 0) return reportError(nn);
    return Py_BuildValue("i",nn);
}

static PyObject*
py_func_newcombo(PyObject *self, PyObject *args)
{
    int type, n, m;
    if (!PyArg_ParseTuple(args, "iii:func_newcombo", &type, &n, &m))
        return NULL;
    int nn = func_new(type, n, m, 0);
    if (nn < 0) return reportError(nn);
    return Py_BuildValue("i",nn);
}

static PyObject*
py_func_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:func_del", &n))
        return NULL;
    int iok = func_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_func_value(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:func_value", &n, &t))
        return NULL;
    double r = func_value(n, t);
    return Py_BuildValue("d",r);
}


static PyMethodDef ct_methods[] = {
    {"func_new", py_func_new, METH_VARARGS},
    {"func_newcombo", py_func_newcombo, METH_VARARGS},
    {"func_del", py_func_del, METH_VARARGS},
    {"func_value", py_func_value, METH_VARARGS},
    {NULL, NULL}
};

extern "C" {

    /* Initialization function for the module */

    DL_EXPORT(void) initctfunc(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctfunc", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}

