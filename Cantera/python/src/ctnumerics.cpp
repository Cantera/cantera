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
#include "ct.h"
#include "ctnum.h" 
#include <string>
#include <vector>
#include <iostream>
using namespace std;


//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h" 

static PyObject *
py_new_matrix(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:ct_newMatrix", &n, &m)) 
        return NULL;
    int nn = newMatrix(n,m);
    if (nn < 0) return reportCanteraError();
    return Py_BuildValue("i",nn);
}

static PyObject*
py_matrix_delete(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:matrix_delete", &n)) return NULL;    
    delMatrix(n);
    return Py_BuildValue("i",0);
}

static PyObject*
py_matrix_newcopy(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:matrix_newcopy", &n)) return NULL;
    return Py_BuildValue("i",matrix_copy(n));
}

static PyObject*
py_matrix_assign(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:matrix_assign", &n, &m)) return NULL;
    return Py_BuildValue("i",matrix_assign(n, m));
}

static PyObject*
py_matrix_nrows(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:matrix_nrows", &n)) return NULL;
    return Py_BuildValue("i",matrix_nRows(n));
}

static PyObject*
py_matrix_ncols(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:matrix_ncols", &n)) return NULL;
    return Py_BuildValue("i",matrix_nColumns(n));
}

static PyObject*
py_matrix_value(PyObject *self, PyObject *args)
{
    int nn,m,n;
    if (!PyArg_ParseTuple(args, "iii:matrix_value", &nn, &m, &n)) return NULL;
    return Py_BuildValue("d",matrix_value(nn,m,n));
}

static PyObject*
py_matrix_setvalue(PyObject *self, PyObject *args)
{
    int nn,m,n;
    double v;
    if (!PyArg_ParseTuple(args, "iiid:matrix_setvalue", &nn, &m, &n, &v)) return NULL;
    return Py_BuildValue("d",matrix_setvalue(nn,m,n,v));
}

static PyObject*
py_matrix_solve(PyObject *self, PyObject *args)
{
    int na, nb;
    if (!PyArg_ParseTuple(args, "ii:matrix_solve", &na, &nb)) return NULL;
    int i = matrix_solve(na, nb);
    if (i == -1) return reportCanteraError();
    return Py_BuildValue("i",i);
}

static PyObject*
py_matrix_mult(PyObject *self, PyObject *args)
{
    int na, nb, np;
    if (!PyArg_ParseTuple(args, "iii:matrix_mult", &na, &nb, &np)) return NULL;
    int i = matrix_multiply(na, nb, np);
    if (i == -1) return reportCanteraError();
    return Py_BuildValue("i",i);
}

static PyObject*
py_matrix_invert(PyObject *self, PyObject *args)
{
    int na;
    if (!PyArg_ParseTuple(args, "i:matrix_invert", &na)) return NULL;
    int i = matrix_invert(na);
    if (i == -1) return reportCanteraError();
    return Py_BuildValue("i",i);
}


static PyObject *
py_bandmatrix_new(PyObject *self, PyObject *args)
{
    int n, kl, ku;
    if (!PyArg_ParseTuple(args, "iii:bandmatrix_new", &n, &kl, &ku)) 
        return NULL;
    int nn = bmatrix_new(n,kl,ku);
    if (nn < 0) return reportCanteraError();
    return Py_BuildValue("i",nn);
}

static PyObject*
py_bandmatrix_delete(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bandmatrix_delete", &n)) return NULL;    
    bmatrix_del(n);
    return Py_BuildValue("i",0);
}

static PyObject*
py_bandmatrix_newcopy(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bandmatrix_newcopy", &n)) return NULL;
    return Py_BuildValue("i",bmatrix_copy(n));
}

static PyObject*
py_bandmatrix_assign(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:bandmatrix_assign", &n, &m)) return NULL;
    return Py_BuildValue("i",bmatrix_assign(n, m));
}

static PyObject*
py_bandmatrix_nrows(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bandmatrix_nrows", &n)) return NULL;
    return Py_BuildValue("i",bmatrix_nRows(n));
}

static PyObject*
py_bandmatrix_ncols(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bandmatrix_ncols", &n)) return NULL;
    return Py_BuildValue("i",bmatrix_nColumns(n));
}

static PyObject*
py_bandmatrix_value(PyObject *self, PyObject *args)
{
    int nn,m,n;
    if (!PyArg_ParseTuple(args, "iii:bandmatrix_value", &nn, &m, &n)) return NULL;
    return Py_BuildValue("d",bmatrix_value(nn,m,n));
}

static PyObject*
py_bandmatrix_setvalue(PyObject *self, PyObject *args)
{
    int nn,m,n;
    double v;
    if (!PyArg_ParseTuple(args, "iiid:bandmatrix_setvalue", &nn, &m, &n, &v)) 
        return NULL;
    return Py_BuildValue("d",bmatrix_setvalue(nn,m,n,v));
}

static PyObject*
py_bandmatrix_solve(PyObject *self, PyObject *args)
{
    int na, nb;
    if (!PyArg_ParseTuple(args, "ii:bandmatrix_solve", &na, &nb)) 

        return NULL;
    int i = bmatrix_solve(na, nb);
    if (i == -1) return reportCanteraError();
    return Py_BuildValue("i",i);
}

static PyObject*
py_bandmatrix_mult(PyObject *self, PyObject *args)
{
    int na, nb, np;
    if (!PyArg_ParseTuple(args, "iii:bandmatrix_mult", &na, &nb, &np)) 
        return NULL;
    int i = bmatrix_multiply(na, nb, np);
    if (i == -1) return reportCanteraError();
    return Py_BuildValue("i",i);
}


// static PyObject*
// num_getarray(PyObject *self, PyObject *args)
// {
//     int n;
//     int job;
    
//     if (!PyArg_ParseTuple(args, "ii:num_getarray", &n, &job)) 
//         return NULL;

//     // array attributes
//     int iok = -22;
//     int nrxns = kin_nReactions(kin);
//     int nsp = phase_nSpecies(kin);
//     vector<double> x;
//     switch (job) {
//     case 1:
//         x.resize(nrxns);
//         iok = kin_getFwdRatesOfProgress(kin, nrxns, x.begin());
//         break;
//     case 2:
//         x.resize(nrxns);
//         iok = kin_getRevRatesOfProgress(kin, nrxns, x.begin());
//         break;
//     case 3:
//         x.resize(nrxns);
//         iok = kin_getEquilibriumConstants(kin, nrxns, x.begin());
//         break;
//     default:
//         ;
//     }
//     if (iok >= 0) {
//         return pyNumericTuple_FromVector(x);
//     }
//     else if (iok == -1) {
//         return reportCanteraError();
//     }
//     else {
//         PyErr_SetString(ErrorObject,"Unknown array attribute");
//         return NULL;
//     }
// }


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"Matrix", py_new_matrix,  METH_VARARGS},
    {"matrix_delete", py_matrix_delete,  METH_VARARGS},
    {"matrix_newcopy", py_matrix_newcopy,  METH_VARARGS},
    {"matrix_assign", py_matrix_assign,  METH_VARARGS},
    {"matrix_nrows", py_matrix_nrows,  METH_VARARGS},
    {"matrix_ncols", py_matrix_ncols,  METH_VARARGS},
    {"matrix_value", py_matrix_value,  METH_VARARGS},
    {"matrix_setvalue", py_matrix_setvalue,  METH_VARARGS},
    {"matrix_solve", py_matrix_solve,  METH_VARARGS},
    {"matrix_mult", py_matrix_mult,  METH_VARARGS},
    {"matrix_invert", py_matrix_invert,  METH_VARARGS},
    {"BandMatrix", py_bandmatrix_new,  METH_VARARGS},
    {"bandmatrix_delete", py_bandmatrix_delete,  METH_VARARGS},
    {"bandmatrix_newcopy", py_bandmatrix_newcopy,  METH_VARARGS},
    {"bandmatrix_assign", py_bandmatrix_assign,  METH_VARARGS},
    {"bandmatrix_nrows", py_bandmatrix_nrows,  METH_VARARGS},
    {"bandmatrix_ncols", py_bandmatrix_ncols,  METH_VARARGS},
    {"bandmatrix_value", py_bandmatrix_value,  METH_VARARGS},
    {"bandmatrix_setvalue", py_bandmatrix_setvalue,  METH_VARARGS},
    {"bandmatrix_solve", py_bandmatrix_solve,  METH_VARARGS},
    {"bandmatrix_mult", py_bandmatrix_mult,  METH_VARARGS},
    //{"setfp", thermo_setfp,  METH_VARARGS},
    //{"equil", thermo_equil,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initctnumerics) */

    DL_EXPORT(void) initctnumerics(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctnumerics", ct_methods);
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}
