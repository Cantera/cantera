
/**
 * @file ctpyrpath.cpp
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
#include "ctrpath.h" 

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h"


static PyObject*
py_rdiag_new(PyObject *self, PyObject *args)
{
    int iok = rdiag_new();
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_rdiag_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rdiag_del", &n))
        return NULL;
    int iok = rdiag_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_detailed(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rdiag_detailed", &n))
        return NULL;
    int iok = rdiag_detailed(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_brief(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rdiag_brief", &n))
        return NULL;
    int iok = rdiag_brief(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setThreshold(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setThreshold", &n, &v))
        return NULL;
    int iok = rdiag_setThreshold(n,v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setBoldColor(PyObject *self, PyObject *args)
{
    int n;
    char* color;
    if (!PyArg_ParseTuple(args, "is:rdiag_setBoldColor", &n, &color))
        return NULL;
    int iok = rdiag_setBoldColor(n, color);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setNormalColor(PyObject *self, PyObject *args)
{
    int n;
    char* color;
    if (!PyArg_ParseTuple(args, "is:rdiag_setNormalColor", &n, &color))
        return NULL;
    int iok = rdiag_setNormalColor(n, color);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setDashedColor(PyObject *self, PyObject *args)
{
    int n;
    char* color;
    if (!PyArg_ParseTuple(args, "is:rdiag_setDashedColor", &n, &color))
        return NULL;
    int iok = rdiag_setDashedColor(n,color);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setDotOptions(PyObject *self, PyObject *args)
{
    int n;
    char* opt;
    if (!PyArg_ParseTuple(args, "is:rdiag_setDotOptions", &n, &opt))
        return NULL;
    int iok = rdiag_setDotOptions(n,opt);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setBoldThreshold(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setBoldThreshold", &n, &v))
        return NULL;
    int iok = rdiag_setBoldThreshold(n,v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setNormalThreshold(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setNormalThreshold", &n, &v))
        return NULL;
    int iok = rdiag_setNormalThreshold(n,v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setLabelThreshold(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setLabelThreshold", &n, &v))
        return NULL;
    int iok = rdiag_setLabelThreshold(n,v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setScale(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setScale", &n, &v))
        return NULL;
    int iok = rdiag_setScale(n,v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setFlowType(PyObject *self, PyObject *args)
{
    int n;
    int iflow;
    if (!PyArg_ParseTuple(args, "ii:rdiag_setFlowType", &n, &iflow))
        return NULL;
    int iok = rdiag_setFlowType(n, iflow);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setArrowWidth(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setArrowWidth", &n, &v))
        return NULL;
    int iok = rdiag_setArrowWidth(n,v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}


static PyObject*
py_rdiag_displayOnly(PyObject *self, PyObject *args)
{
    int n, k;
    if (!PyArg_ParseTuple(args, "ii:rdiag_displayOnly", &n, &k))
        return NULL;
    int iok = rdiag_displayOnly(n,k);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setTitle(PyObject *self, PyObject *args)
{
    int n;
    char* t;
    if (!PyArg_ParseTuple(args, "is:rdiag_setTitle", &n, &t))
        return NULL;
    int iok = rdiag_setTitle(n,t);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_add(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:rdiag_add", &n, &m))
        return NULL;
    int iok = rdiag_add(n,m);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_findMajor(PyObject *self, PyObject *args)
{
    int n, m;
    double thresh;
    PyObject* a;
    if (!PyArg_ParseTuple(args, "idO:rdiag_findMajor", &n, &thresh, &a))
        return NULL;
    PyArrayObject* aa = (PyArrayObject*)a;
    int lda = aa->dimensions[0];
    double* x = (double*)aa->data;
    int iok = rdiag_findMajor(n, thresh, lda, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_write(PyObject *self, PyObject *args)
{
    int n, fmt;
    char* nm;
    if (!PyArg_ParseTuple(args, "iis:rdiag_write", &n, &fmt, &nm))
        return NULL;
    int iok = rdiag_write(n, fmt, nm);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rbuild_new(PyObject *self, PyObject *args)
{
    int iok = rbuild_new();
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_rbuild_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rbuild_del", &n))
        return NULL;
    int iok = rbuild_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rbuild_init(PyObject *self, PyObject *args)
{
    int n;
    char* log;
    int k;
    if (!PyArg_ParseTuple(args, "isi:rbuild_init", &n, &log, &k))
        return NULL;
    int iok = rbuild_init(n,log,k);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_rbuild_build(PyObject *self, PyObject *args)
{
    int n;
    int k, idiag, iquiet;
    char *el, *dotfile;
    if (!PyArg_ParseTuple(args, "iissii:rbuild_build", &n, &k, 
            &el, &dotfile, &idiag, &iquiet))
        return NULL;
    int iok = rbuild_build(n,k,el,dotfile,idiag,iquiet);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyMethodDef ct_methods[] = {
    {"rdiag_setDashedColor", py_rdiag_setDashedColor, METH_VARARGS},
    {"rbuild_new", py_rbuild_new, METH_VARARGS},
    {"rdiag_write", py_rdiag_write, METH_VARARGS},
    {"rdiag_setDotOptions", py_rdiag_setDotOptions, METH_VARARGS},
    {"rdiag_setScale", py_rdiag_setScale, METH_VARARGS},
    {"rdiag_setTitle", py_rdiag_setTitle, METH_VARARGS},
    {"rdiag_setArrowWidth", py_rdiag_setArrowWidth, METH_VARARGS},
    {"rdiag_displayOnly", py_rdiag_displayOnly, METH_VARARGS},
    {"rdiag_setThreshold", py_rdiag_setThreshold, METH_VARARGS},
    {"rdiag_setBoldThreshold", py_rdiag_setBoldThreshold, METH_VARARGS},
    {"rdiag_new", py_rdiag_new, METH_VARARGS},
    {"rdiag_del", py_rdiag_del, METH_VARARGS},
    {"rdiag_detailed", py_rdiag_detailed, METH_VARARGS},
    {"rdiag_add", py_rdiag_add, METH_VARARGS},
    {"rdiag_findMajor", py_rdiag_findMajor, METH_VARARGS},
    {"rbuild_build", py_rbuild_build, METH_VARARGS},
    {"rdiag_setNormalThreshold", py_rdiag_setNormalThreshold, METH_VARARGS},
    {"rdiag_brief", py_rdiag_brief, METH_VARARGS},
    {"rbuild_del", py_rbuild_del, METH_VARARGS},
    {"rdiag_setNormalColor", py_rdiag_setNormalColor, METH_VARARGS},
    {"rbuild_init", py_rbuild_init, METH_VARARGS},
    {"rdiag_setBoldColor", py_rdiag_setBoldColor, METH_VARARGS},
    {"rdiag_setFlowType", py_rdiag_setFlowType, METH_VARARGS},
    {"rdiag_setLabelThreshold", py_rdiag_setLabelThreshold, METH_VARARGS},
    {NULL, NULL}
};

extern "C" {

    /* Initialization function for the module */

    DL_EXPORT(void) initctrpath(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctrpath", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}

