
/**
 * @file ctpyxml.cpp
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
#include <iostream>
using namespace std;

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h"


static PyObject*
py_xml_new(PyObject *self, PyObject *args)
{
    char* nm;
    if (!PyArg_ParseTuple(args, "s:xml_new", &nm))
        return NULL;
    int n = xml_new(nm);
    return Py_BuildValue("i",n);
}

static PyObject*
py_xml_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_del", &n))
        return NULL;
    int iok = xml_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_xml_build(PyObject *self, PyObject *args)
{
    int n;
    char* file;
    if (!PyArg_ParseTuple(args, "is:xml_build", &n, &file))
        return NULL;
    int iok = xml_build(n, file);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_xml_attrib(PyObject *self, PyObject *args)
{
    int n;
    char *key;
    if (!PyArg_ParseTuple(args, "is:xml_attrib", &n, &key))
        return NULL;
    char* val = new char[81];
    int iok = xml_attrib(n, key, val);
    if (iok < 0) return reportError(iok);
    PyObject* r = Py_BuildValue("s",val);
    delete val;
    return r;
}

static PyObject*
py_xml_tag(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_tag", &n))
        return NULL;
    char* val = new char[81];
    int iok = xml_tag(n, val);
    if (iok < 0) return reportError(iok);
    PyObject* r = Py_BuildValue("s",val);
    delete val;
    return r;
}

static PyObject*
py_xml_value(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_value", &n))
        return NULL;
    char* val = new char[81];
    int iok = xml_value(n, val);
    if (iok < 0) return reportError(iok);
    PyObject* r = Py_BuildValue("s",val);
    delete val;
    return r;
}

static PyObject*
py_xml_child(PyObject *self, PyObject *args)
{
    int n;
    char* loc;
    if (!PyArg_ParseTuple(args, "is:xml_child", &n, &loc))
        return NULL;
    int m = xml_child(n, loc);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_childbynumber(PyObject *self, PyObject *args)
{
    int n, k;
    if (!PyArg_ParseTuple(args, "ii:xml_childbynumber", &n, &k))
        return NULL;
    int m = xml_child_bynumber(n, k);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_findID(PyObject *self, PyObject *args)
{
    int n;
    char* id;
    if (!PyArg_ParseTuple(args, "is:xml_findID", &n, &id))
        return NULL;
    int m = xml_findID(n, id);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_findByName(PyObject *self, PyObject *args)
{
    int n;
    char* nm;
    if (!PyArg_ParseTuple(args, "is:xml_findID", &n, &nm))
        return NULL;
    int m = xml_findByName(n, nm);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_nChildren(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_nChildren", &n))
        return NULL;
    int m = xml_nChildren(n);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_addChild(PyObject *self, PyObject *args)
{
    int n;
    char *name, *value;
    if (!PyArg_ParseTuple(args, "iss:xml_addChild", &n, &name, &value))
        return NULL;
    int m = xml_addChild(n, name, value);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}


static PyObject*
py_xml_addChildNode(PyObject *self, PyObject *args)
{
    int n, j;
    if (!PyArg_ParseTuple(args, "ii:xml_addChildNode", &n, &j))
        return NULL;
    int m = xml_addChildNode(n, j);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}


static PyObject*
py_xml_addAttrib(PyObject *self, PyObject *args)
{
    int n;
    char *name, *value;
    if (!PyArg_ParseTuple(args, "iss:xml_addAttrib", &n, &name, &value))
        return NULL;
    int m = xml_addAttrib(n, name, value);
    if (m < 0) return reportError(m);
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_removeChild(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:xml_removeChild", &n, &m))
        return NULL;
    int iok = xml_removeChild(n, m);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}


static PyObject*
py_xml_write(PyObject *self, PyObject *args)
{
    int n;
    char *file;
    if (!PyArg_ParseTuple(args, "is:xml_write", &n, &file))
        return NULL;
    int iok = xml_write(n, file);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_ctml_getFloatArray(PyObject *self, PyObject *args)
{
    int n;
    int iconv, ia;
    if (!PyArg_ParseTuple(args, "iii", &n, &iconv, &ia))
        return NULL;
    
    PyArrayObject* a = 
        (PyArrayObject*)PyArray_FromDims(1, &ia, PyArray_DOUBLE);
    double* x = (double*)a->data;
    int iok = ctml_getFloatArray(n, ia, x, iconv);
    if (iok < 0) return reportError(iok);
    return PyArray_Return(a);
}

static PyMethodDef ct_methods[] = {
    {"xml_attrib", py_xml_attrib, METH_VARARGS},
    {"xml_addAttrib", py_xml_addAttrib, METH_VARARGS},
    {"xml_tag", py_xml_tag, METH_VARARGS},
    {"xml_value", py_xml_value, METH_VARARGS},
    {"xml_new", py_xml_new, METH_VARARGS},
    {"xml_del", py_xml_del, METH_VARARGS},
    {"xml_build", py_xml_build, METH_VARARGS},
    {"xml_child", py_xml_child, METH_VARARGS},
    {"xml_childbynumber", py_xml_childbynumber, METH_VARARGS},
    {"xml_findID", py_xml_findID, METH_VARARGS},
    {"xml_findByName", py_xml_findByName, METH_VARARGS},
    {"xml_nChildren", py_xml_nChildren, METH_VARARGS},
    {"xml_addChild", py_xml_addChild, METH_VARARGS},
    {"xml_addChildNode", py_xml_addChildNode, METH_VARARGS},
    {"xml_removeChild", py_xml_removeChild, METH_VARARGS},
    {"xml_write", py_xml_write, METH_VARARGS},
    {"ctml_getFloatArray", py_ctml_getFloatArray, METH_VARARGS},
    {NULL, NULL}
};

extern "C" {

    /* Initialization function for the module */

    DL_EXPORT(void) initctxml(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctxml", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}

