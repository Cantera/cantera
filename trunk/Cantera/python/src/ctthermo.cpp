/**
 * @file ctthermo.cpp
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


//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h" 

static PyObject *
ct_newThermoFromXML(PyObject *self, PyObject *args)
{
    int mxml;
    //char* id;
    if (!PyArg_ParseTuple(args, "i:ct_newThermoFromXML", &mxml)) 
        return NULL;
    int n = newThermoFromXML(mxml);
    if (n < 0) return reportCanteraError();
    //int p = th_phase(n);
    return Py_BuildValue("i",n);
}

static PyObject*
thermo_delete(PyObject *self, PyObject *args)
{
    int th;
    if (!PyArg_ParseTuple(args, "i:thermo_delete", &th)) 
        return NULL;    
    delThermo(th);
    return Py_BuildValue("i",0);
}

static PyObject*
thermo_index(PyObject *self, PyObject *args) {
    char* id;
    if (!PyArg_ParseTuple(args, "s:index", &id)) return NULL;
    return Py_BuildValue("i",th_thermoIndex(id));        
}

static PyObject*
thermo_refpressure(PyObject *self, PyObject *args) {
    int th;
    if (!PyArg_ParseTuple(args, "i:refpressure", &th)) return NULL;
    return Py_BuildValue("d",th_refPressure(th));        
}

static PyObject*
thermo_mintemp(PyObject *self, PyObject *args) {
    int th, k;
    if (!PyArg_ParseTuple(args, "ii:mintemp", &th, &k)) return NULL;
    return Py_BuildValue("d",th_minTemp(th,k));        
}

static PyObject*
thermo_maxtemp(PyObject *self, PyObject *args) {
    int th, k;
    if (!PyArg_ParseTuple(args, "ii:maxtemp", &th, &k)) return NULL;
    return Py_BuildValue("d",th_maxTemp(th,k));        
}

// static PyObject*
// thermo_geteos(PyObject *self, PyObject *args) {
//     char *fname, *id;
//     if (!PyArg_ParseTuple(args, "ss:geteos", &fname, &id)) return NULL;
//     return Py_BuildValue("i",get_eos(fname, id));        
// }

static PyObject*
thermo_import(PyObject *self, PyObject *args) {
    int n, mxml;
    char* id;
    if (!PyArg_ParseTuple(args, "iis:import", &n, &mxml, &id)) return NULL;
    int iok = import_phase(n, mxml, id);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}


static PyObject*
thermo_getfp(PyObject *self, PyObject *args)
{
    double vv;
    bool ok = true;
    int th;
    int job;
    
    if (!PyArg_ParseTuple(args, "ii:thermo_getfp", &th, &job)) 
        return NULL;

    // floating-point attributes
    switch (job) {
    case 1:
        vv = th_enthalpy_mole(th); break; 
    case 2:
        vv = th_intEnergy_mole(th); break;
    case 3:
        vv = th_entropy_mole(th); break;
    case 4:
        vv = th_gibbs_mole(th); break;
    case 5:
        vv = th_cp_mole(th); break;
    case 6:
        vv = th_cv_mole(th); break;
    case 7:
        vv = th_pressure(th); break;
    case 8:
        vv = th_enthalpy_mass(th); break; 
    case 9:
        vv = th_intEnergy_mass(th); break;
    case 10:
        vv = th_entropy_mass(th); break;
    case 11:
        vv = th_gibbs_mass(th); break;
    case 12:
        vv = th_cp_mass(th); break;
    case 13:
        vv = th_cv_mass(th); break;
    default:
        ok = false;
    }
    if (ok) {
        if (vv == -999.999) {
            return reportCanteraError();
        }
        return Py_BuildValue("d",vv);        
    }
    else {
        PyErr_SetString(ErrorObject,"Unknown floating-point attribute");
        return NULL;
    }        
}

static PyObject*
thermo_setfp(PyObject *self, PyObject *args)
{
    double v1, v2;
    int iok = -2;
    int th;
    int job;
    
    if (!PyArg_ParseTuple(args, "iidd:thermo_setfp", &th, &job, &v1, &v2)) 
        return NULL;

    //vector_fp v(2);
    double* v = new double[2];
    v[0] = v1; v[1] = v2;

    // set floating-point attributes
    switch (job) {
    case 1:
        iok = th_setPressure(th, v1); break; 
    case 2:
        iok = th_set_HP(th, v); break;
    case 3:
        iok = th_set_UV(th, v); break;
    case 4:
        iok = th_set_SV(th, v); break;
    case 5:
        iok = th_set_SP(th, v); break;
    default:
        iok = -10; 
    }
    delete v;
    if (iok >= 0) 
        return Py_BuildValue("i",iok);
    if (iok == -1) return reportCanteraError();
    else {
        PyErr_SetString(ErrorObject,"Error in thermo_setfp");
        return NULL;
    }        
}


static PyObject*
thermo_getarray(PyObject *self, PyObject *args)
{
    int th;
    int job;

    if (!PyArg_ParseTuple(args, "ii:thermo_getarray", &th, &job)) 
        return NULL;

    int nsp = th_nSpecies(th);

    // array attributes
    int iok = -22;
    PyArrayObject* x = 
        (PyArrayObject*)PyArray_FromDims(1, &nsp, PyArray_DOUBLE);
    double* xd = (double*)x->data;
    switch (job) {
    case 20:
        iok = th_chemPotentials(th,nsp,xd);
        break;
    case 23:
        iok = th_getEnthalpies_RT(th,nsp,xd);
        break;
    case 24:
        iok = th_getEntropies_R(th,nsp,xd);
        break;
    case 25:
        iok = th_getCp_R(th,nsp,xd);
        break;

    default:
        ;
    }
    if (iok >= 0) {
        return PyArray_Return(x);
    }
    else {
        PyErr_SetString(ErrorObject,"Unknown array attribute");
        return NULL;
    }
}


static PyObject*
thermo_equil(PyObject *self, PyObject *args)
{
    int iok = -2;
    int th;
    int XY;
    
    if (!PyArg_ParseTuple(args, "ii:thermo_equil", &th, &XY)) 
        return NULL;

    iok = th_equil(th, XY);
    if (iok >= 0) 
        return Py_BuildValue("i",iok);
    if (iok == -1) return reportCanteraError();
    else {
        PyErr_SetString(ErrorObject,"Error in thermo_equil");
        return NULL;
    }
}

/* List of functions defined in the module */

static PyMethodDef ct_methods[] = {
    {"ThermoFromXML", ct_newThermoFromXML,  METH_VARARGS},
    {"delete", thermo_delete,  METH_VARARGS},
    {"mintemp", thermo_mintemp,  METH_VARARGS},
    {"maxtemp", thermo_maxtemp,  METH_VARARGS},
    {"thermoIndex", thermo_index,  METH_VARARGS},
    {"refpressure", thermo_refpressure,  METH_VARARGS},
    {"getfp", thermo_getfp,  METH_VARARGS},
    {"setfp", thermo_setfp,  METH_VARARGS},
    {"getarray", thermo_getarray,  METH_VARARGS},
    {"equil", thermo_equil,  METH_VARARGS},
    {"import_xml", thermo_import,  METH_VARARGS},
    {NULL,  NULL}		   /* sentinel */
};


extern "C" {

    /* Initialization function for the module (*must* be called initctthermo) */

    DL_EXPORT(void) initctthermo(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctthermo", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}
