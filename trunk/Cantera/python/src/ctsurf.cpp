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

// static PyObject*
// py_flow_newcopy(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:flow_newcopy", &n)) return NULL;
//     return Py_BuildValue("i",flow_copy(n));
// }

// static PyObject*
// py_flow_assign(PyObject *self, PyObject *args)
// {
//     int n, m;
//     if (!PyArg_ParseTuple(args, "ii:flow_assign", &n, &m)) return NULL;
//     return Py_BuildValue("i",flow_assign(n, m));
// }

// static PyObject*
// py_flow_readinputs(PyObject *self, PyObject *args)
// {
//     int n;
//     char* infile;
//     if (!PyArg_ParseTuple(args, "is:flow_readinputs", &n, &infile)) 
//         return NULL;
//     return Py_BuildValue("i",flow_readinputs(n,infile));
// }

static PyObject*
py_flow_setupgrid(PyObject *self, PyObject *args)
{
    int n;
    PyObject* grid;
    if (!PyArg_ParseTuple(args, "iO:flow_setupgrid", &n, &grid)) 
        return NULL;

    vector<double> z;
    int iok = pyNumericSequence_ToVector(grid, z);
    if (iok == -1) {
        PyErr_SetString(ErrorObject, "Third argument must be a sequence");
        return NULL;
    }
    if (iok == -2) {
        PyErr_SetString(ErrorObject, 
            "Sequence must contain only numeric values");
        return NULL;
    }
    iok = flow_setupgrid(n, z.size(), z.begin());
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
py_flow_setinletstate(PyObject *self, PyObject *args)
{
    int nn,gas;
    if (!PyArg_ParseTuple(args, "ii:flow_setinletstate", &nn, &gas)) 
        return NULL;
    return Py_BuildValue("i",flow_setinletstate(nn,gas));
}

static PyObject*
py_flow_setinlet_u(PyObject *self, PyObject *args)
{
    int nn;
    double u;
    if (!PyArg_ParseTuple(args, "id:flow_setinlet_u", &nn, &u)) 
        return NULL;
    return Py_BuildValue("i",flow_setinlet_u(nn,u));
}

static PyObject*
py_flow_setinlet_v(PyObject *self, PyObject *args)
{
    int nn;
    double v;
    if (!PyArg_ParseTuple(args, "id:flow_setinlet_v", &nn, &v)) 
        return NULL;
    return Py_BuildValue("i",flow_setinlet_v(nn,v));
}

static PyObject*
py_flow_setsurface_t(PyObject *self, PyObject *args)
{
    int nn;
    double t;
    if (!PyArg_ParseTuple(args, "id:flow_setsurface_t", &nn, &t)) 
        return NULL;
    return Py_BuildValue("i",flow_setsurface_t(nn,t));
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
    int n, j, slen;
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
py_flow_solve(PyObject *self, PyObject *args)
{
    int n, loglevel;
    PyObject *s, *snew;
    if (!PyArg_ParseTuple(args, "iOOi:flow_solve", &n, &s, &snew, &loglevel)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)s)->data;
    double* xnew = (double*)((PyArrayObject*)snew)->data;
    int iok = flow_solve(n,x,xnew,loglevel);
    if (iok == -1) return reportCanteraError();
    return Py_BuildValue("i",iok);
}

static PyObject*
py_flow_ssnorm(PyObject *self, PyObject *args)
{
    int n, loglevel;
    PyObject *ps, *pr;
    if (!PyArg_ParseTuple(args, "iOO:flow_solve", &n, &ps, &pr)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)ps)->data;
    double* r = (double*)((PyArrayObject*)pr)->data;
    double ss = flow_ssnorm(n,x,r);
    return Py_BuildValue("d",ss);
}


static PyObject*
py_flow_timeinteg(PyObject *self, PyObject *args)
{
    int n, nsteps, loglevel;
    double dt;
    PyObject *s, *snew;
    if (!PyArg_ParseTuple(args, "iidOOi:flow_solve", &n, &nsteps, &dt, 
            &s, &snew, &loglevel)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)s)->data;
    double* xnew = (double*)((PyArrayObject*)snew)->data;
    int iok = flow_timeinteg(n,nsteps,dt,x,xnew,loglevel);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
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


static PyObject*
py_flow_eval(PyObject *self, PyObject *args)
{
    int n, j;
    PyObject *px, *pr;
    if (!PyArg_ParseTuple(args, "iiOO:flow_solve", &n, &j, &px, &pr)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    double* r = (double*)((PyArrayObject*)pr)->data;
    int iok = flow_eval(n, j, x, r);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}


static PyObject*
py_flow_integratechem(PyObject *self, PyObject *args)
{
    int n, j;
    PyObject *px;
    double dt;
    if (!PyArg_ParseTuple(args, "iOd:flow_solve", &n, &px, &dt)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    int iok = flow_integratechem(n, x, dt);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}

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
    PyObject *px;
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
py_flow_setnewtonoptions(PyObject *self, PyObject *args)
{
    int n, age;
    double ratio;
    if (!PyArg_ParseTuple(args, "iid:flow_setnewtonoptions", &n, &age, &ratio)) 
        return NULL;
    int iok = flow_setnewtonoptions(n, age, ratio);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_flow_setsteadymode(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flow_setsteadymode", &n)) 
        return NULL;
    int iok = flow_setsteadymode(n);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_flow_settransientmode(PyObject *self, PyObject *args)
{
    int n;
    double dt;
    PyArrayObject* px;
    if (!PyArg_ParseTuple(args, "idO:flow_settransientmode", &n, &dt, &px)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    int iok = flow_settransientmode(n, dt, x);
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

static PyObject*
py_flow_save(PyObject *self, PyObject *args)
{
    int n;
    char *fname, *id;
    PyArrayObject* px;
    if (!PyArg_ParseTuple(args, "issO:flow_solve", &n, &fname, &id, &px)) 
        return NULL;
    double* x = (double*)((PyArrayObject*)px)->data;
    int iok = flow_save(n, fname, id, x);
    return Py_BuildValue("i",iok);
}

static PyObject*
py_flow_restore(PyObject *self, PyObject *args)
{
    int n, job, iz, isoln;
    char *fname, *id;
    PyArrayObject *pz, *psoln;
    if (!PyArg_ParseTuple(args, "iiss:flow_solve", 
            &n, &job, &fname, &id)) 
        return NULL;
    int iok;
    double *z=0, *soln=0;

    iok = flow_restore(n, -1, fname, id, iz, z, isoln, soln);
    if (job < 0) {
        return Py_BuildValue("(ii)",iz,isoln);
    }
    pz = (PyArrayObject*)PyArray_FromDims(1, &iz, PyArray_DOUBLE);
    vector<int> sdim(2);
    sdim[0] = iz;
    sdim[1] = isoln/iz;
    psoln = (PyArrayObject*)PyArray_FromDims(2, sdim.begin(), PyArray_DOUBLE);
    z = (double*)((PyArrayObject*)pz)->data;
    soln = (double*)((PyArrayObject*)psoln)->data;
    iok = flow_restore(n, 0, fname, id, iz, z, isoln, soln);

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
    int itype, ip;
    if (!PyArg_ParseTuple(args, "ii:bdry_new", &itype, &ip)) 
        return NULL;
    int nn = bdry_new(itype,ip);
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


/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"Flow", py_flow_new,  METH_VARARGS},
    {"flow_delete", py_flow_delete,  METH_VARARGS},
    {"flow_setupgrid", py_flow_setupgrid,  METH_VARARGS},
    {"flow_setthermo", py_flow_setthermo,  METH_VARARGS},
    {"flow_setkinetics", py_flow_setkinetics,  METH_VARARGS},
    {"flow_settransport", py_flow_settransport,  METH_VARARGS},
    {"flow_setpressure", py_flow_setpressure,  METH_VARARGS},
    {"flow_setinletstate", py_flow_setinletstate,  METH_VARARGS},
    {"flow_setinlet_u", py_flow_setinlet_u,  METH_VARARGS},
    {"flow_setinlet_v", py_flow_setinlet_v,  METH_VARARGS},
    {"flow_setsurface_t", py_flow_setsurface_t,  METH_VARARGS},
    {"flow_solvespecies", py_flow_solvespecies,  METH_VARARGS},
    {"flow_settemperature", py_flow_settemperature,  METH_VARARGS},
    {"flow_setmassfraction", py_flow_setmassfraction,  METH_VARARGS},
    {"flow_settolerances", py_flow_settolerances,  METH_VARARGS},
    {"flow_energy", py_flow_energy,  METH_VARARGS},
    {"flow_showsolution", py_flow_showsolution,  METH_VARARGS},
    {"flow_eval", py_flow_eval,  METH_VARARGS},
    {"flow_solve", py_flow_solve,  METH_VARARGS},
    {"flow_timeinteg", py_flow_timeinteg,  METH_VARARGS},
    {"flow_integratechem", py_flow_integratechem,  METH_VARARGS},
    {"flow_setnewtonoptions", py_flow_setnewtonoptions,  METH_VARARGS},
    {"flow_resize", py_flow_resize,  METH_VARARGS},
    {"flow_outputtec", py_flow_outputtec,  METH_VARARGS},
    {"flow_setsteadymode", py_flow_setsteadymode,  METH_VARARGS},
    {"flow_settransientmode", py_flow_settransientmode,  METH_VARARGS},
    {"flow_save", py_flow_save,  METH_VARARGS},
    {"flow_restore", py_flow_restore,  METH_VARARGS},
    {"flow_setfixedpoint", py_flow_setfixedpoint,  METH_VARARGS},
    {"flow_setboundaries", py_flow_setboundaries,  METH_VARARGS},
    {"flow_ssnorm", py_flow_ssnorm,  METH_VARARGS},
    {"copy", py_copy,  METH_VARARGS},
    {"bdry_new", py_bdry_new,  METH_VARARGS},
    {"bdry_del", py_bdry_delete,  METH_VARARGS},
    {"bdry_set", py_bdry_set,  METH_VARARGS},
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
