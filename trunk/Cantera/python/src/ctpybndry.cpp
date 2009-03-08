
/**
 * @file ctpybndry.cpp
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
#include "ctbdry.h" 

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h"
#include <iostream>
using namespace std; 

static PyObject*
py_bndry_new(PyObject *self, PyObject *args)
{
    int itype;
    if (!PyArg_ParseTuple(args, "i:bndry_new", &itype))
        return NULL;
    int iok = bndry_new(itype);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_bndry_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bndry_del", &n))
        return NULL;
    bndry_del(n);
    return Py_BuildValue("i",0);
}

static PyObject*
py_bndry_temperature(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bndry_temperature", &n))
        return NULL;
    double t = bndry_temperature(n);
    return Py_BuildValue("d",t);
}

static PyObject*
py_bndry_settemperature(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:bndry_settemperature", &n, &t))
        return NULL;
    int iok = bndry_settemperature(n, t);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_bndry_spreadrate(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bndry_spreadrate", &n))
        return NULL;
    double v = bndry_spreadrate(n);
    return Py_BuildValue("d",v);
}

static PyObject*
py_bndry_setspreadrate(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:bndry_setspreadrate", &n, &v))
        return NULL;
    int iok = bndry_setspreadrate(n, v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_bndry_mdot(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bndry_mdot", &n))
        return NULL;
    double mdot = bndry_mdot(n);
    return Py_BuildValue("d",mdot);
}

static PyObject*
py_bndry_setmdot(PyObject *self, PyObject *args)
{
    int n;
    double mdot;
    if (!PyArg_ParseTuple(args, "id:bndry_setmdot", &n, &mdot))
        return NULL;
    int iok = bndry_setmdot(n, mdot);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_bndry_setxinbyname(PyObject *self, PyObject *args)
{
    int n;
    char* xin;
    //PyObject* o;
    if (!PyArg_ParseTuple(args, "is:bndry_setxin", &n, &xin))
        return NULL;
    //double* x = (double*)((PyArrayObject*)o)->data;
    int iok = bndry_setxinbyname(n, xin);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}


static PyObject*
py_bndry_setxin(PyObject *self, PyObject *args)
{
    int n;
    PyObject* xin;
    if (!PyArg_ParseTuple(args, "iO:bndry_setxin", &n, &xin))
        return NULL;
    double* x = (double*)((PyArrayObject*)xin)->data;
    int iok = bndry_setxin(n, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}


static PyMethodDef ct_methods[] = {
    {"bndry_temperature", py_bndry_temperature, METH_VARARGS},
    {"bndry_setxin", py_bndry_setxin, METH_VARARGS},
    {"bndry_setxinbyname", py_bndry_setxinbyname, METH_VARARGS},
    {"bndry_settemperature", py_bndry_settemperature, METH_VARARGS},
    {"bndry_setspreadrate", py_bndry_setspreadrate, METH_VARARGS},
    {"bndry_spreadrate", py_bndry_spreadrate, METH_VARARGS},
    {"bndry_new", py_bndry_new, METH_VARARGS},
    {"bndry_del", py_bndry_del, METH_VARARGS},
    {"bndry_mdot", py_bndry_mdot, METH_VARARGS},
    {"bndry_setmdot", py_bndry_setmdot, METH_VARARGS},
    {NULL, NULL}
};

extern "C" {

    /* Initialization function for the module */

    DL_EXPORT(void) initctbndry(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctbndry", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}

