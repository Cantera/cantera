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
#include "ctstagn.h" 
#include <string>
#include <vector>
#include <iostream>

using namespace std;

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h" 

static PyObject *
py_flow_new(PyObject *self, PyObject *args)
{
    int itype, iph, np;
    if (!PyArg_ParseTuple(args, "iii:flow_new", 
            &itype, &iph, &np)) 
        return NULL;
    int nn = flow_new(itype,iph,np);
    if (nn < 0) return reportError(nn);
    return Py_BuildValue("i",nn);
}

static PyObject*
py_flow_delete(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flow_delete", &n)) return NULL;    
    flow_del(n);
    return Py_BuildValue("i",0);
}

static PyObject*
py_flow_setupgrid(PyObject *self, PyObject *args)
{
    int n;
    PyObject* grid;
    if (!PyArg_ParseTuple(args, "iO:flow_setupgrid", &n, &grid)) 
        return NULL;

    //vector<double> z;
    //int iok = pyNumericSequence_ToVector(grid, z);
    PyArrayObject* g = (PyArrayObject*)grid;
    double* xd = (double*)g->data;
    int glen = g->dimensions[0];
    int 
    //    if (iok == -1) {
    //    PyErr_SetString(ErrorObject, "Third argument must be a sequence");
    //    return NULL;
    //}
    //if (iok == -2) {
    //    PyErr_SetString(ErrorObject, 
    //        "Sequence must contain only numeric values");
    //    return NULL;
    //}

    iok = flow_setupgrid(n, glen, xd);
    if (iok == -1) return reportCanteraError();
    else if (iok < 0) return NULL;
    return Py_BuildValue("i",0);
}

static PyObject*
py_flow_setkinetics(PyObject *self, PyObject *args)
{
    int nn,k;
    if (!PyArg_ParseTuple(args, "ii:flow_setkinetics", &nn, &k)) return NULL;
    return Py_BuildValue("i",flow_setkinetics(nn,k));
}

static PyObject*
py_flow_settransport(PyObject *self, PyObject *args)
{
    int nn,k,soret;
    if (!PyArg_ParseTuple(args, "iii:flow_settransport", 
            &nn, &k, &soret)) return NULL;
    return Py_BuildValue("i",flow_settransport(nn,k,soret));
}

static PyObject*
py_flow_setthermo(PyObject *self, PyObject *args)
{
    int nn,k;
    if (!PyArg_ParseTuple(args, "ii:flow_setthermo", &nn, &k)) return NULL;
    return Py_BuildValue("i",flow_setthermo(nn,k));
}

static PyObject*
py_flow_setpressure(PyObject *self, PyObject *args)
{
    int nn;
    double p;
    if (!PyArg_ParseTuple(args, "id:flow_setpressure", &nn, &p)) 
        return NULL;
    return Py_BuildValue("i",flow_setpressure(nn,p));
}

static PyObject*
py_flow_settemperature(PyObject *self, PyObject *args)
{
    int n, j;
    double t;
    if (!PyArg_ParseTuple(args, "iid:flow_settemperature", &n, &j, &t)) 
        return NULL;
    return Py_BuildValue("i",flow_settemperature(n,j,t));
}

static PyObject*
py_flow_setenergyfactor(PyObject *self, PyObject *args)
{
    int n;
    double e;
    if (!PyArg_ParseTuple(args, "id:flow_setenergyfactor", &n, &e)) 
        return NULL;
    return Py_BuildValue("i",flow_setenergyfactor(n,e));
}

static PyObject*
py_flow_setmassfraction(PyObject *self, PyObject *args)
{
    int n, j, k;
    double y;
    if (!PyArg_ParseTuple(args, "iiid:flow_setinlet_v", &n, &j, &k, &y)) 
        return NULL;
    return Py_BuildValue("i",flow_setmassfraction(n,j,k,y));
}

static PyObject*
py_flow_showsolution(PyObject *self, PyObject *args)
{
    int n;
    char* fname;
    PyObject* soln;
    if (!PyArg_ParseTuple(args, "isO:flow_showsolution", &n, &fname, &soln)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)soln)->data;
    return Py_BuildValue("i",flow_showsolution(n,fname,x));
}

static PyObject*
py_flow_solvespecies(PyObject *self, PyObject *args)
{
    int n, slen;
    PyObject* s;
    if (!PyArg_ParseTuple(args, "iiO:flow_solvespecies", &n, &slen, &s)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)s)->data;
    for (int i = 0; i < slen; i++) {
        if (x[i] <= 0.0) 
            flow_fixspecies(n, i);
        else
            flow_solvespecies(n, i);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_copy(PyObject *self, PyObject *args)
{
    int n;
    PyObject *s, *snew;
    if (!PyArg_ParseTuple(args, "iOO:copy", &n, &s, &snew)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)s)->data;
    double* xnew = (double*)((PyArrayObject*)snew)->data;
    for (int i = 0; i < n; i++) xnew[i] = x[i];
    return Py_BuildValue("i",0);
}


static PyObject*
py_flow_settolerances(PyObject *self, PyObject *args)
{
    int n, nr, na;
    PyObject *prtol, *patol;
    if (!PyArg_ParseTuple(args, "iiOiO:flow_solve", &n, &nr, &prtol, 
            &na, &patol)) 
        return NULL;
    double* rtol = (double*)((PyArrayObject*)prtol)->data;
    double* atol = (double*)((PyArrayObject*)patol)->data;
    int iok = flow_settolerances(n, nr, rtol, na, atol);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

// static PyObject*
// py_flow_integratechem(PyObject *self, PyObject *args)
// {
//     int n;
//     PyObject *px;
//     double dt;
//     if (!PyArg_ParseTuple(args, "iOd:flow_solve", &n, &px, &dt)) 
//         return NULL;
//     double* x = (double*)((PyArrayObject*)px)->data;
//     int iok = flow_integratechem(n, x, dt);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",iok);
// }

static PyObject*
py_flow_outputtec(PyObject *self, PyObject *args)
{
    int n;
    PyObject *px;
    char *fname, *title;
    int zone;
    if (!PyArg_ParseTuple(args, "iOssi:flow_outputtec", &n, &px, 
            &fname, &title, &zone)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    int iok = flow_outputtec(n, x, fname, title, zone);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_flow_resize(PyObject *self, PyObject *args)
{
    int n, points;
    if (!PyArg_ParseTuple(args, "ii:flow_resize", &n, &points)) 
        return NULL;
    int iok = flow_resize(n, points);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_flow_energy(PyObject *self, PyObject *args)
{
    int n, j, flag, iok;
    if (!PyArg_ParseTuple(args, "iii:flow_energy", &n, &j, &flag)) 
        return NULL;
    if (flag == 1)
        iok = flow_solveenergyeqn(n, j);
    else
        iok = flow_fixtemperature(n, j);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_flow_setfixedpoint(PyObject *self, PyObject *args)
{
    int n, j0;
    double t0;
    if (!PyArg_ParseTuple(args, "iid:flow_setfixedpoint", &n, &j0, &t0)) 
        return NULL;
    int iok = flow_setfixedpoint(n, j0, t0);
    return Py_BuildValue("i",iok);
}

// static PyObject*
// py_flow_eval(PyObject *self, PyObject *args)
// {
//     int n, j;
//     PyObject *px, *pr;
//     if (!PyArg_ParseTuple(args, "iiOO:flow_eval", &n, &j, &px, &pr)) 
//         return NULL;
//     double* x = (double*)((PyArrayObject*)px)->data;
//     double* r = (double*)((PyArrayObject*)pr)->data;
//     int iok = flow_eval(n, j, x, r);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",iok);
// }

static PyObject*
py_flow_restore(PyObject *self, PyObject *args)
{
    int n, job, iz, isoln;
    char *fname, *id;
    PyArrayObject *pz, *psoln;
    if (!PyArg_ParseTuple(args, "iiss:flow_restore", 
            &n, &job, &fname, &id)) 
        return NULL;
    int iok;
    double *z=0, *soln=0;
    iok = flow_restore(n, -1, fname, id, iz, z, isoln, soln);
    if (iok < 0) return reportError(iok);
    if (job < 0) {
        return Py_BuildValue("(ii)",iz,isoln);
    }
    pz = (PyArrayObject*)PyArray_FromDims(1, &iz, PyArray_DOUBLE);
    int sdim[2];
    sdim[0] = iz;
    sdim[1] = isoln/iz;
    psoln = (PyArrayObject*)PyArray_FromDims(2, sdim, PyArray_DOUBLE);
    z = (double*)((PyArrayObject*)pz)->data;
    soln = (double*)((PyArrayObject*)psoln)->data;
    iok = flow_restore(n, 0, fname, id, iz, z, isoln, soln);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("(OO)",pz,psoln);
}

static PyObject*
py_flow_setboundaries(PyObject *self, PyObject *args)
{
    int n, nleft, nright;
    if (!PyArg_ParseTuple(args, "iii:flow_setboundaries", &n, &nleft, 
            &nright)) 
        return NULL;
    int iok = flow_setboundaries(n, nleft, nright);
    return Py_BuildValue("i",iok);
}


/* flow boundary objects */


static PyObject *
py_bdry_new(PyObject *self, PyObject *args)
{
    int itype, ip, kin;
    if (!PyArg_ParseTuple(args, "iii:bdry_new", &itype, &ip, &kin)) 
        return NULL;
    int nn = bdry_new(itype,ip,kin);
    if (nn < 0) return reportError(nn);
    return Py_BuildValue("i",nn);
}

static PyObject*
py_bdry_delete(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bdry_delete", &n)) return NULL;    
    bdry_del(n);
    return Py_BuildValue("i",0);
}

static PyObject*
py_bdry_set(PyObject *self, PyObject *args)
{
    int n, i;
    double v;
    PyObject* px;
    double* x;
    if (!PyArg_ParseTuple(args, "iidO:bdry_set", &n, &i, &v, &px)) 
        return NULL;
    if (i < 4)
        bdry_set(n, i, &v);
    else {
        x = (double*)((PyArrayObject*)px)->data;
        bdry_set(n, i, x);
    }
    return Py_BuildValue("i",0);
}

static PyObject *
py_onedim_new(PyObject *self, PyObject *args)
{
    int n;
    PyObject *pydom, *pytype;
    if (!PyArg_ParseTuple(args, "iOO:onedim_new", &n, &pydom, &pytype)) 
        return NULL;
    int* dom = (int*)((PyArrayObject*)pydom)->data;
    int* typ = (int*)((PyArrayObject*)pytype)->data;
    int nn = onedim_new(n, dom, typ);
    if (nn == -1) return reportCanteraError();
    return Py_BuildValue("i",nn);
}

static PyObject*
py_onedim_delete(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:onedim_delete", &n)) return NULL;    
    onedim_del(n);
    return Py_BuildValue("i",0);
}

static PyObject*
py_onedim_solve(PyObject *self, PyObject *args)
{
    int n, loglevel;
    PyObject *s, *snew;
    if (!PyArg_ParseTuple(args, "iOOi:onedim_solve", &n, &s, &snew, &loglevel)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)s)->data;
    double* xnew = (double*)((PyArrayObject*)snew)->data;
    int iok = onedim_solve(n,x,xnew,loglevel);
    if (iok == -1) return reportCanteraError();
    return Py_BuildValue("i",iok);
}


static PyObject*
py_onedim_ssnorm(PyObject *self, PyObject *args)
{
    int n;
    PyObject *ps, *pr;
    if (!PyArg_ParseTuple(args, "iOO:flow_solve", &n, &ps, &pr)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)ps)->data;
    double* r = (double*)((PyArrayObject*)pr)->data;
    double ss = onedim_ssnorm(n,x,r);
    return Py_BuildValue("d",ss);
}

static PyObject*
py_onedim_setnewtonoptions(PyObject *self, PyObject *args)
{
    int n, age;
    if (!PyArg_ParseTuple(args, "ii:onedim_setnewtonoptions", &n, &age)) 
        return NULL;
    int iok = onedim_setnewtonoptions(n, age);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_onedim_setsteadymode(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:onedim_setsteadymode", &n)) 
        return NULL;
    int iok = onedim_setsteadymode(n);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_onedim_settransientmode(PyObject *self, PyObject *args)
{
    int n;
    double dt;
    PyArrayObject* px;
    if (!PyArg_ParseTuple(args, "idO:onedim_settransientmode", &n, &dt, &px)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    int iok = onedim_settransientmode(n, dt, x);
    return Py_BuildValue("i",iok);
}


static PyObject*
py_onedim_eval(PyObject *self, PyObject *args)
{
    int n;
    PyObject *px, *pr;
    if (!PyArg_ParseTuple(args, "iOO:onedim_eval", &n, &px, &pr)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    double* r = (double*)((PyArrayObject*)pr)->data;
    int iok = onedim_eval(n, x, r);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_onedim_addflow(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:onedim_addflow", &n, &m)) 
        return NULL;
    int iok = onedim_addFlow(n, m);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

// static PyObject*
// py_onedim_addsurf(PyObject *self, PyObject *args)
// {
//     int n, m;
//     if (!PyArg_ParseTuple(args, "ii:onedim_addflow", &n, &m)) 
//         return NULL;
//     int iok = onedim_addSurf(n, m);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",iok);
// }


static PyObject*
py_onedim_resize(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:onedim_resize", &n)) 
        return NULL;
    int iok = onedim_resize(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_onedim_writestats(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:onedim_writeStats", &n)) 
        return NULL;
    int iok = onedim_writeStats(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_onedim_timestep(PyObject *self, PyObject *args)
{
    int n, nsteps, loglevel;
    double dt;
    PyObject *px, *pr;
    if (!PyArg_ParseTuple(args, "iidOOi:onedim_timestep", &n, &nsteps,
            &dt, &px, &pr, &loglevel)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    double* r = (double*)((PyArrayObject*)pr)->data;
    double newdt = onedim_timestep(n, nsteps, dt, x, r, loglevel);
    if (newdt < 0.0) return reportError(-1);
    return Py_BuildValue("d",newdt);
}

static PyObject*
py_onedim_save(PyObject *self, PyObject *args)
{
    int n;
    char *fname, *id, *desc;
    PyArrayObject* px;
    if (!PyArg_ParseTuple(args, "isssO:onedim_save", &n, &fname, &id, &desc, &px)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    int iok = onedim_save(n, fname, id, desc, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"Flow", py_flow_new,  METH_VARARGS},
    {"flow_delete", py_flow_delete,  METH_VARARGS},
    {"flow_setupgrid", py_flow_setupgrid,  METH_VARARGS},
    {"flow_setthermo", py_flow_setthermo,  METH_VARARGS},
    {"flow_setkinetics", py_flow_setkinetics,  METH_VARARGS},
    {"flow_settransport", py_flow_settransport,  METH_VARARGS},
    {"flow_setpressure", py_flow_setpressure,  METH_VARARGS},
    //{"flow_setinletstate", py_flow_setinletstate,  METH_VARARGS},
    //    {"flow_setinlet_u", py_flow_setinlet_u,  METH_VARARGS},
    //{"flow_setinlet_v", py_flow_setinlet_v,  METH_VARARGS},
    //{"flow_setsurface_t", py_flow_setsurface_t,  METH_VARARGS},
    {"flow_solvespecies", py_flow_solvespecies,  METH_VARARGS},
    {"flow_settemperature", py_flow_settemperature,  METH_VARARGS},
    {"flow_setenergyfactor", py_flow_setenergyfactor,  METH_VARARGS},
    {"flow_setmassfraction", py_flow_setmassfraction,  METH_VARARGS},
    {"flow_settolerances", py_flow_settolerances,  METH_VARARGS},
    {"flow_energy", py_flow_energy,  METH_VARARGS},
    {"flow_showsolution", py_flow_showsolution,  METH_VARARGS},
    //    {"flow_integratechem", py_flow_integratechem,  METH_VARARGS},
    {"flow_resize", py_flow_resize,  METH_VARARGS},
    {"flow_outputtec", py_flow_outputtec,  METH_VARARGS},
    //    {"flow_eval", py_flow_eval,  METH_VARARGS},
    {"flow_restore", py_flow_restore,  METH_VARARGS},
    {"flow_setfixedpoint", py_flow_setfixedpoint,  METH_VARARGS},
    {"flow_setboundaries", py_flow_setboundaries,  METH_VARARGS},
    {"copy", py_copy,  METH_VARARGS},
    {"bdry_new", py_bdry_new,  METH_VARARGS},
    {"bdry_del", py_bdry_delete,  METH_VARARGS},
    {"bdry_set", py_bdry_set,  METH_VARARGS},
    {"onedim_solve", py_onedim_solve,  METH_VARARGS},
    {"onedim_new", py_onedim_new,  METH_VARARGS},
    {"onedim_del", py_onedim_delete,  METH_VARARGS},
    {"onedim_setnewtonoptions", py_onedim_setnewtonoptions,  METH_VARARGS},
    {"onedim_ssnorm", py_onedim_ssnorm,  METH_VARARGS},
    {"onedim_setsteadymode", py_onedim_setsteadymode,  METH_VARARGS},
    {"onedim_settransientmode", py_onedim_settransientmode,  METH_VARARGS},
    {"onedim_eval", py_onedim_eval,  METH_VARARGS},
    {"onedim_addflow", py_onedim_addflow,  METH_VARARGS},
    //{"onedim_addsurf", py_onedim_addsurf,  METH_VARARGS},
    {"onedim_resize", py_onedim_resize,  METH_VARARGS},
    {"onedim_writestats", py_onedim_writestats,  METH_VARARGS},
    {"onedim_timestep", py_onedim_timestep,  METH_VARARGS},
    {"onedim_save", py_onedim_save,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initctflow) */

    DL_EXPORT(void) initctflow(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctflow", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}
