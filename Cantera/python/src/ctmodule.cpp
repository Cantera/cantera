/**
 * @file ctmodule.cpp
 * Cantera Python Interface
 *
 */


// turn off warnings about long names under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

//  #define USE_DL_EXPORT

#include "Python.h"
#include "ct.h"

//  constants defined in the module
static PyObject *ErrorObject;

#include <iostream>
#include <string>
using namespace std;

// local includes
#include "pyutils.h" 

static PyObject *
ct_buildSolutionFromXML(PyObject *self, PyObject *args)
{
    int ixml, ith, ikin;
    char *src=0, *id=0;
    if (!PyArg_ParseTuple(args, "sisii:buildSolutionFromXML", &src, &ixml, 
            &id, &ith, &ikin)) 
        return NULL;
    int ok = buildSolutionFromXML(src, ixml, id, ith, ikin);
    if (ok == -1) { return reportCanteraError();}
    return Py_BuildValue("i",ok); 
}

static PyObject *
ct_get_cantera_error(PyObject *self, PyObject *args)
{
    char* buf = new char[400];
    getCanteraError(400, buf);
    PyObject* msg = Py_BuildValue("s",buf);
    delete buf;
    return msg;
}

static PyObject *
ct_print(PyObject *self, PyObject *args)
{
    char* msg;
    if (!PyArg_ParseTuple(args, "s:print", &msg)) 
        return NULL;
    printf(msg);
    return Py_BuildValue("i",0);
}

static PyObject *
ct_readlog(PyObject *self, PyObject *args)
{
    char* msg = 0;
    int n = readlog(-1, msg);
    if (n > 0) {
        msg = new char[n+1];
        int ok = readlog(n, msg);
        PyObject* r = Py_BuildValue("s",msg);
        return r;
    }
    else 
        return Py_BuildValue("s","");
}

static PyObject *
ct_ck2ctml(PyObject *self, PyObject *args)
{
    int iok;
    char *infile, *thermo, *tran, *outfile, *idtag;
    if (!PyArg_ParseTuple(args, "sssss:ck2ctml", &infile, 
            &thermo, &tran, &outfile, &idtag)) 
        return NULL;
    iok = ck_to_ctml(infile, thermo, tran, outfile, idtag);
    if (iok == -1) { return reportCanteraError();}
    return Py_BuildValue("i",iok); 
}

/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"get_Cantera_Error", ct_get_cantera_error,  METH_VARARGS},
    {"ct_print", ct_print,  METH_VARARGS},
    {"readlog", ct_readlog,  METH_VARARGS},
    {"ck2ctml", ct_ck2ctml,  METH_VARARGS},
    {"buildSolutionFromXML", ct_buildSolutionFromXML, METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initctmodule) */

    DL_EXPORT(void) initctmodule(void)
    {
        PyObject *m, *d;
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctmodule", ct_methods);
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }
}
