/**
 * @file ctsurf.cpp
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
#include "ctsurf.h" 
#include <string>
#include <vector>
#include <iostream>

using namespace std;

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h" 

static PyObject*
py_surf_setsitedensity(PyObject *self, PyObject *args)
{
    int n;
    double s0;
    if (!PyArg_ParseTuple(args, "id:surf_setsitedensity", &n, &s0)) 
        return NULL;

    int iok = surf_setsitedensity(n, s0);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_surf_sitedensity(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:surf_sitedensity", &n)) 
        return NULL;

    double s0 = surf_sitedensity(n);
    return Py_BuildValue("d",s0);
}


static PyObject*
py_surf_setcoverages(PyObject *self, PyObject *args)
{
    int n;
    PyObject* cov;
    if (!PyArg_ParseTuple(args, "iO:surf_setcoverages", &n, &cov)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)cov)->data;
    int iok = surf_setcoverages(n, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_surf_setconcentrations(PyObject *self, PyObject *args)
{
    int n;
    PyObject* c;
    if (!PyArg_ParseTuple(args, "iO:surf_setconcentrations", &n, &c)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)c)->data;
    int iok = surf_setconcentrations(n, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_surf_getcoverages(PyObject *self, PyObject *args)
{
    int n;
    PyArrayObject* cov;
    if (!PyArg_ParseTuple(args, "i:surf_getcoverages", &n)) 
        return NULL;
    int nsp = th_nSpecies(n);
    cov = (PyArrayObject*)PyArray_FromDims(1, &nsp, PyArray_DOUBLE);
    double* x = (double*)((PyArrayObject*)cov)->data;
    int iok = surf_getcoverages(n, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("O",cov);
}

static PyObject*
py_surf_getconcentrations(PyObject *self, PyObject *args)
{
    int n;
    PyArrayObject* c;
    if (!PyArg_ParseTuple(args, "i:surf_getconcentrations", &n)) 
        return NULL;
    int nsp = th_nSpecies(n);
    c = (PyArrayObject*)PyArray_FromDims(1, &nsp, PyArray_DOUBLE);
    double* x = (double*)((PyArrayObject*)c)->data;
    int iok = surf_getconcentrations(n, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("O",c);
}

// static PyObject*
// py_surf_doc(PyObject *self, PyObject *args)
// {
//     int n;
//     char *k, *v;
//     if (!PyArg_ParseTuple(args, "iss:surf_doc", &n, &k, &v)) 
//         return NULL;
//     surf_doc(n, k, v);
//     return Py_BuildValue("i",0);
// }

// static PyObject *
// py_surfkin_new(PyObject *self, PyObject *args)
// {
//     int isurf, ib1, ib2;
//     if (!PyArg_ParseTuple(args, "iii:surfkin_new", &isurf, &ib1, &ib2)) 
//         return NULL;
//     int nn = surfkin_new(isurf, ib1, ib2);
//     if (nn < 0) return reportError(nn);
//     return Py_BuildValue("i",nn);
// }

// static PyObject*
// py_surfkin_delete(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:surfkin_delete", &n)) return NULL;    
//     surf_del(n);
//     return Py_BuildValue("i",0);
// }


// static PyObject *
// py_surfkin_addreaction(PyObject *self, PyObject *args)
// {
//     int isurf, nr, np, nrate;
//     PyObject *pyr, *pyrst, *pyro, *pyp, *pypst, *pyrate;
//     if (!PyArg_ParseTuple(args, "iiOOOiOOiO:surfkin_addreaction", 
//             &isurf, &nr, &pyr, &pyrst, &pyro, &np, &pyp, &pypst, 
//             &nrate, &pyrate)) 
//         return NULL;
//     int* r = (int*)((PyArrayObject*)pyr)->data;
//     int* rst = (int*)((PyArrayObject*)pyrst)->data;
//     int* ro = (int*)((PyArrayObject*)pyro)->data;
//     int* p = (int*)((PyArrayObject*)pyp)->data;
//     int* pst = (int*)((PyArrayObject*)pypst)->data;
//     double* rate = (double*)((PyArrayObject*)pyrate)->data;
//     int nn = surfkin_addreaction(isurf, nr, r, rst, ro, np, p, pst, 
//         nrate, rate);
//     if (nn < 0) return reportError(nn);
//     return Py_BuildValue("i",nn);
// }

// static PyObject *
// py_surfkin_nreactions(PyObject *self, PyObject *args)
// {
//     int isurf;
//     if (!PyArg_ParseTuple(args, "i:surfkin_nreactions", &isurf)) 
//         return NULL;
//     int nn = surfkin_nreactions(isurf);
//     return Py_BuildValue("i",nn);
// }

// static PyObject*
// py_surfkin_getratesofprogress(PyObject *self, PyObject *args)
// {
//     int n;
//     PyArrayObject* rop;
//     if (!PyArg_ParseTuple(args, "i:surfkin_getratesofprogress", &n)) 
//         return NULL;
//     int nrxn = surfkin_nreactions(n);
//     rop = (PyArrayObject*)PyArray_FromDims(1, &nrxn, PyArray_DOUBLE);
//     double* x = (double*)((PyArrayObject*)rop)->data;
//     int iok = surfkin_getratesofprogress(n, x);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("O",rop);
// }

// static PyObject*
// py_surfkin_getnetproductionrates(PyObject *self, PyObject *args)
// {
//     int n, nsp;
//     PyArrayObject* sdot;
//     if (!PyArg_ParseTuple(args, "ii:surfkin_getnetproductionrates", &n, &nsp)) 
//         return NULL;
//     sdot = (PyArrayObject*)PyArray_FromDims(1, &nsp, PyArray_DOUBLE);
//     double* x = (double*)((PyArrayObject*)sdot)->data;
//     int iok = surfkin_getnetproductionrates(n, x);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("O",sdot);
// }

// static PyObject *
// py_surfkin_integrate(PyObject *self, PyObject *args)
// {
//     int isurf;
//     double dt;
//     if (!PyArg_ParseTuple(args, "id:surfkin_integrate", &isurf, &dt)) 
//         return NULL;
//     int nn = surfkin_integrate(isurf, dt);
//     if (nn < 0) return reportError(nn);
//     return Py_BuildValue("i",0);
// }

// static PyObject *
// py_surfkin_save(PyObject *self, PyObject *args)
// {
//     int isurf;
//     char *fname, *id, *comment;
//     if (!PyArg_ParseTuple(args, "isss:surfkin_save", &isurf, &fname, &id, &comment)) 
//         return NULL;
//     int nn = surfkin_save(isurf, fname, id, comment);
//     if (nn < 0) return reportError(nn);
//     return Py_BuildValue("i",0);
// }

// static PyObject *
// py_surf1d_new(PyObject *self, PyObject *args)
// {
//     int ikin;
//     if (!PyArg_ParseTuple(args, "i:surf1d_new", &ikin)) 
//         return NULL;
//     int nn = surface_new(ikin);
//     if (nn < 0) return reportError(nn);
//     return Py_BuildValue("i",nn);
// }

// static PyObject*
// py_surf1d_delete(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:surf1d_delete", &n)) return NULL;    
//     surface_del(n);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_surf1d_settolerances(PyObject *self, PyObject *args)
// {
//     int n, nr, na;
//     PyObject *prtol, *patol;
//     if (!PyArg_ParseTuple(args, "iiOiO:surf1d_settolerances", &n, &nr, &prtol, 
//             &na, &patol)) 
//         return NULL;
//     double* rtol = (double*)((PyArrayObject*)prtol)->data;
//     double* atol = (double*)((PyArrayObject*)patol)->data;
//     int iok = surface_settolerances(n, nr, rtol, na, atol);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",iok);
// }

// static PyObject*
// py_surf1d_settemperature(PyObject *self, PyObject *args)
// {
//     int n;
//     double t;
//     if (!PyArg_ParseTuple(args, "id:surf1d_settemperature", &n, &t)) 
//         return NULL;
//     int iok = surface_settemperature(n, t);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",iok);
// }

// static PyObject*
// py_surf1d_temperature(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:surf1d_temperature", &n)) 
//         return NULL;
//     return Py_BuildValue("d",surface_temperature(n));
// }

// static PyObject*
// py_surf1d_setcoverages(PyObject *self, PyObject *args)
// {
//     int n;
//     PyObject* cov;
//     if (!PyArg_ParseTuple(args, "iO:surf1d_setcoverages", &n, &cov)) 
//         return NULL;
//     double* x = (double*)((PyArrayObject*)cov)->data;
//     int iok = surface_setcoverages(n, x);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_surf1d_setmultiplier(PyObject *self, PyObject *args)
// {
//     int n, k;
//     double f;
//     if (!PyArg_ParseTuple(args, "iid:surf1d_setmultiplier", &n, &k, &f)) 
//         return NULL;
//     int iok = surface_setmultiplier(n, k, f);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",iok);
// }

// static PyObject*
// py_surf1d_multiplier(PyObject *self, PyObject *args)
// {
//     int n, k;
//     if (!PyArg_ParseTuple(args, "ii:surf1d_multiplier", &n, &k)) 
//         return NULL;
//     return Py_BuildValue("d",surface_multiplier(n, k));
// }

// static PyObject*
// py_surf1d_fixspecies(PyObject *self, PyObject *args)
// {
//     int n, k;
//     double c;
//     if (!PyArg_ParseTuple(args, "iid:surf1d_fixspecies", &n, &k, &c)) 
//         return NULL;
//     int iok = surface_fixspecies(n, k, c);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",iok);
// }

// static PyObject*
// py_surf1d_solvespecies(PyObject *self, PyObject *args)
// {
//     int n, slen;
//     PyObject* s;    
//     if (!PyArg_ParseTuple(args, "iiO:surf1d_solvespecies", &n, &slen, &s)) 
//         return NULL;
//     double* x = (double*)((PyArrayObject*)s)->data;
//     for (int i = 0; i < slen; i++) {
//         if (x[i] <= 0.0) 
//             surface_fixspecies(n, i);
//         else
//             surface_solvespecies(n, i);
//     }
//     return Py_BuildValue("i",0);
// }


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    //    {"surf_new", py_surf_new,  METH_VARARGS},
    //{"surf_delete", py_surf_delete,  METH_VARARGS},
    //{"surf_addspecies", py_surf_addspecies,  METH_VARARGS},
    {"surf_setsitedensity", py_surf_setsitedensity,  METH_VARARGS},
    {"surf_sitedensity", py_surf_sitedensity,  METH_VARARGS},
    //{"surf_nspecies", py_surf_nspecies,  METH_VARARGS},
    {"surf_setcoverages", py_surf_setcoverages,  METH_VARARGS},
    {"surf_getcoverages", py_surf_getcoverages,  METH_VARARGS},
    {"surf_setconcentrations", py_surf_setconcentrations,  METH_VARARGS},
    {"surf_getconcentrations", py_surf_getconcentrations,  METH_VARARGS},
    //{"surf_doc", py_surf_doc,  METH_VARARGS},
//     {"surfkin_new", py_surfkin_new,  METH_VARARGS},
//     {"surfkin_delete", py_surfkin_delete,  METH_VARARGS},
//     {"surfkin_addreaction", py_surfkin_addreaction,  METH_VARARGS},
//     {"surfkin_nreactions", py_surfkin_nreactions,  METH_VARARGS},
//     {"surfkin_getratesofprogress", py_surfkin_getratesofprogress,  METH_VARARGS},
//     {"surfkin_getsdot", py_surfkin_getnetproductionrates,  METH_VARARGS},
//     {"surfkin_integrate", py_surfkin_integrate,  METH_VARARGS},
//     {"surfkin_save", py_surfkin_save,  METH_VARARGS},
//     {"surf1d_new", py_surf1d_new,  METH_VARARGS},
//     {"surf1d_delete", py_surf1d_delete,  METH_VARARGS},
//     {"surf1d_settolerances", py_surf1d_settolerances,  METH_VARARGS},
//     {"surf1d_settemperature", py_surf1d_settemperature,  METH_VARARGS},
//     {"surf1d_temperature", py_surf1d_temperature,  METH_VARARGS},
//     {"surf1d_setcoverages", py_surf1d_setcoverages,  METH_VARARGS},
//     {"surf1d_setmultiplier", py_surf1d_setmultiplier,  METH_VARARGS},
//     {"surf1d_multiplier", py_surf1d_multiplier,  METH_VARARGS},
//     {"surf1d_fixspecies", py_surf1d_fixspecies,  METH_VARARGS},
//     {"surf1d_solvespecies", py_surf1d_solvespecies,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initctsurf) */

    DL_EXPORT(void) initctsurf(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctsurf", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}
