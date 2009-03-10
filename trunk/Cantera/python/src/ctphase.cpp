/**
 * @file ctphase.cpp
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
static PyObject *OneAtmos;
static PyObject *GasCon;

// local includes
#include "pyutils.h" 



// /**
//  * Create a new Phase object.
//  */
// static PyObject *
// ct_newPhase(PyObject *self, PyObject *args) {
//     int n = newPhase();
//     return Py_BuildValue("i",n);
// }


// /**
//  * Delete the Phase object.
//  */
// static PyObject*
// phase_delete(PyObject *self, PyObject *args)
// {
//     int ph;
//     if (!PyArg_ParseTuple(args, "i:phase_delete", &ph)) 
//         return NULL;    
//     delPhase(ph);
//     return Py_BuildValue("i",0);
// }


static PyObject*
py_temperature(PyObject *self, PyObject *args) {
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_temperature", &ph)) return NULL;
    return Py_BuildValue("d",phase_temperature(ph));        
}

static PyObject*
py_density(PyObject *self, PyObject *args) {
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_density", &ph)) return NULL;
    return Py_BuildValue("d",phase_density(ph));        
}

static PyObject*
py_molardensity(PyObject *self, PyObject *args) {
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_molardensity", &ph)) return NULL;
    return Py_BuildValue("d",phase_molarDensity(ph));        
}

static PyObject*
py_meanmolwt(PyObject *self, PyObject *args) {
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_meanmolwt", &ph)) return NULL;
    return Py_BuildValue("d",phase_meanMolecularWeight(ph));        
}

static PyObject*
py_molefraction(PyObject *self, PyObject *args) {
    int ph, k;
    if (!PyArg_ParseTuple(args, "ii:py_molefraction", &ph, &k)) return NULL;
    return Py_BuildValue("d",phase_moleFraction(ph, k));        
}

static PyObject*
py_massfraction(PyObject *self, PyObject *args) {
    int ph, k;
    if (!PyArg_ParseTuple(args, "ii:py_massfraction", &ph, &k)) return NULL;
    return Py_BuildValue("d",phase_massFraction(ph, k));        
}

static PyObject*
py_nelements(PyObject *self, PyObject *args) {
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_nelements", &ph)) return NULL;
    return Py_BuildValue("i",phase_nElements(ph));        
}

static PyObject*
py_nspecies(PyObject *self, PyObject *args) {
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_nspecies", &ph)) return NULL;
    return Py_BuildValue("i",phase_nSpecies(ph));        
}

static PyObject*
py_natoms(PyObject *self, PyObject *args) {
    int ph, k, m;
    if (!PyArg_ParseTuple(args, "iii:py_natoms", &ph, &k, &m)) return NULL;
    return Py_BuildValue("d",phase_nAtoms(ph, k, m));        
}

static PyObject*
py_addelement(PyObject *self, PyObject *args) {
    int ph;
    char* name;
    double wt;
    if (!PyArg_ParseTuple(args, "isd:py_addelement", &ph, &name, &wt)) 
        return NULL;
    int ok = phase_addElement(ph, name, wt);
    if (ok < 0) return reportError(ok);
    else return Py_BuildValue("i",0);
}

static PyObject*
py_elementindex(PyObject *self, PyObject *args) {
    int ph;
    char* nm;
    if (!PyArg_ParseTuple(args, "is:py_elementindex", &ph, &nm)) return NULL;
    int k = phase_elementIndex(ph,nm);
    if (k >= 0) 
        return Py_BuildValue("i",k);
    else {
        PyErr_SetString(ErrorObject,(
                            "Unknown element ("+string(nm)+")").c_str());
        return NULL;
    }
}

static PyObject*
py_speciesindex(PyObject *self, PyObject *args) {
    int ph;
    char* nm;
    if (!PyArg_ParseTuple(args, "is:py_speciesindex", &ph, &nm)) return NULL;
    int k = phase_speciesIndex(ph,nm);
    if (k >= 0) 
        return Py_BuildValue("i",k);
    else {
        PyErr_SetString(ErrorObject,(
                            "Unknown species ("+string(nm)+")").c_str());
        return NULL;
    }
}

static PyObject*
py_report(PyObject *self, PyObject *args) {
    int th, show_thermo;
    int buflen = 400;
    char* output_buf = new char[buflen];    
    if (!PyArg_ParseTuple(args, "ii:py_report", &th, &show_thermo)) 
        return NULL;
    int iok = phase_report(th, buflen, output_buf, show_thermo);
    if (iok < -1 && iok != -999) {
        delete output_buf;
        output_buf = new char[-iok];
        iok = phase_report(th, -iok, output_buf, show_thermo);
    }
    if (iok < 0) return reportError(iok);
    PyObject* s = Py_BuildValue("s",output_buf);
    delete output_buf;
    return s;
}


static PyObject*
phase_getarray(PyObject *self, PyObject *args)
{
    int ph;
    int job;
    
    if (!PyArg_ParseTuple(args, "ii:phase_getarray", &ph, &job)) 
        return NULL;

    // array attributes
    int iok = -22;
    PyArrayObject* x = 0;
    double* xd = 0;
    if (job > 10) {

        int nsp = phase_nSpecies(ph);
        x = (PyArrayObject*)PyArray_FromDims(1, &nsp, PyArray_DOUBLE);
        xd = (double*)x->data;
        switch (job) {
        case 20:
            iok = phase_getMoleFractions(ph,nsp,xd);
            break;
        case 21:
            iok = phase_getMassFractions(ph,nsp,xd);
            break;
        case 22:
            iok = phase_getMolecularWeights(ph,nsp,xd);
            break;
        default:
            ;
        }
    }
    else {

        int nel = phase_nElements(ph);
        x = (PyArrayObject*)PyArray_FromDims(1, &nel, PyArray_DOUBLE);
        xd = (double*)x->data;
        switch (job) {
        case 1:
            iok = phase_getAtomicWeights(ph,nel,xd);
            break;
        default:
            ;
        }
    }

    if (iok >= 0) {
        return PyArray_Return(x);
    }
    else {
        PyErr_SetString(ErrorObject,"Unknown array attribute");
        return NULL;
    }
}

// string attributes 
static PyObject*
phase_getstring(PyObject *self, PyObject *args)
{
    int ph, job, iok = -1;
    int k;
    int buflen;
    char* output_buf = 0;
    if (!PyArg_ParseTuple(args, "iii:phase_getstring", &ph, &job, &k)) 
        return NULL;
    switch (job) {
    case 1:
        buflen = 20;
        output_buf = new char[buflen];
        iok = phase_getElementName(ph, k, buflen, output_buf);
        break;
    case 2:
        buflen = 40;
        output_buf = new char[buflen];
        iok = phase_getSpeciesName(ph, k, buflen, output_buf);
        break;
    default:
        iok = -10;
    }
    if (iok >= 0) {
        PyObject* str = Py_BuildValue("s",output_buf);
        delete output_buf;
        return str;
    }
    delete output_buf;
    if (iok == -1) 
        return reportCanteraError();
    else {
        PyErr_SetString(ErrorObject,"Unknown string attribute");
        return NULL;
    }        
}

static PyObject*
phase_setfp(PyObject *self, PyObject *args)
{
    double vv;
    int iok = -2;
    int ph;
    int job;
    
    if (!PyArg_ParseTuple(args, "iid:phase_getfp", &ph, &job, &vv)) 
        return NULL;

    // set floating-point attributes
    switch (job) {
    case 1:
        iok = phase_setTemperature(ph, vv); break; 
    case 2:
        iok = phase_setDensity(ph, vv); break;
    default:
        iok = -10; 
    }
    if (iok >= 0) 
        return Py_BuildValue("i",iok);        
    else {
        PyErr_SetString(ErrorObject,"Unknown floating-point attribute");
        return NULL;
    }        
}


static PyObject*
phase_setarray(PyObject *self, PyObject *args)
{
    int ph;
    int job;
    int norm;
    int iok;
    PyObject* seq;
    if (!PyArg_ParseTuple(args, "iiiO:phase_setarray", &ph, &job, &norm, &seq)) 
        return NULL;
    //vector_fp v;
    PyArrayObject* a = (PyArrayObject*)seq;
    //iok = pyNumericSequence_ToVector(seq, v);
    double* xd = (double*)a->data;
    int len = a->dimensions[0];
    //if (iok == -1) {
    //    PyErr_SetString(ErrorObject, "Fourth argument must be a sequence");
    //    return NULL;
    //}
    switch (job) {
    case 1:
        iok = phase_setMoleFractions(ph, len, xd, norm); 
        break;
    case 2:
        iok = phase_setMassFractions(ph, len, xd, norm); 
        break;
    default:
        iok = -10;
    }
    if (iok >= 0)
        return Py_BuildValue("i",iok);
    if (iok == -1) 
        return reportCanteraError();
    else {
        PyErr_SetString(ErrorObject, "Error in phase_setarray");
        return NULL;        
    }
}

static PyObject*
phase_setstring(PyObject *self, PyObject *args)
{
    int ph;
    int job;
    int iok;
    char* str;
    if (!PyArg_ParseTuple(args, "iis:phase_setstring", &ph, &job, &str)) 
        return NULL;
    switch (job) {
    case 1:
        iok = phase_setMoleFractionsByName(ph, str); 
        break;
    case 2:
        iok = phase_setMassFractionsByName(ph, str); 
        break;
    default:
        iok = -10;
    }
    if (iok >= 0)
        return Py_BuildValue("i",iok);
    if (iok == -1) 
        return reportCanteraError();
    else {
        PyErr_SetString(ErrorObject, "Error in phase_setstring");
        return NULL;        
    }
}

/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"temperature", py_temperature,  METH_VARARGS},
    {"density", py_density,  METH_VARARGS},
    {"molardensity", py_molardensity,  METH_VARARGS},
    {"meanmolwt", py_meanmolwt,  METH_VARARGS},
    {"molefraction", py_molefraction,  METH_VARARGS},
    {"massfraction", py_massfraction,  METH_VARARGS},
    {"nelements", py_nelements,  METH_VARARGS},
    {"nspecies", py_nspecies,  METH_VARARGS},
    {"natoms", py_natoms,  METH_VARARGS},
    {"addelement", py_addelement,  METH_VARARGS},
    {"elementindex", py_elementindex,  METH_VARARGS},
    {"speciesindex", py_speciesindex,  METH_VARARGS},
    {"getarray", phase_getarray,  METH_VARARGS},
    {"getstring", phase_getstring,  METH_VARARGS},
    {"setfp", phase_setfp,  METH_VARARGS},
    {"setarray", phase_setarray,  METH_VARARGS},
    {"setstring", phase_setstring,  METH_VARARGS},
    {"report", py_report,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initctphase) */

    DL_EXPORT(void) initctphase(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctphase", ct_methods);
        import_array();

        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }
}
