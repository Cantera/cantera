
/**
 * @file ctpyreactor.cpp
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
#include "ctreactor.h" 

//  constants defined in the module
static PyObject *ErrorObject;

// local includes
#include "pyutils.h"


static PyObject*
py_reactor_new(PyObject *self, PyObject *args)
{
    int type;
    if (!PyArg_ParseTuple(args, "i:reactor_new", &type))
        return NULL;
    int n = reactor_new(type);
    return Py_BuildValue("i",n);
}

static PyObject*
py_reactor_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_del", &n))
        return NULL;
    int iok = reactor_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setInitialVolume(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:reactor_setInitialVolume", &n, &v))
        return NULL;
    int iok = reactor_setInitialVolume(n,v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setInitialTime(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:reactor_setInitialTime", &n, &t))
        return NULL;
    int iok = reactor_setInitialTime(n, t);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setEnergy(PyObject *self, PyObject *args)
{
    int n, eflag;
    if (!PyArg_ParseTuple(args, "ii:reactor_setEnergy", &n, &eflag))
        return NULL;
    int iok = reactor_setEnergy(n, eflag);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setThermoMgr(PyObject *self, PyObject *args)
{
    int n;
    int th;
    if (!PyArg_ParseTuple(args, "ii:reactor_setThermoMgr", &n, &th))
        return NULL;
    int iok = reactor_setThermoMgr(n, th);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setKineticsMgr(PyObject *self, PyObject *args)
{
    int n;
    int kin;
    if (!PyArg_ParseTuple(args, "ii:reactor_setKineticsMgr", &n, &kin))
        return NULL;
    int iok = reactor_setKineticsMgr(n, kin);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_advance(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:reactor_advance", &n, &t))
        return NULL;
    int iok = reactor_advance(n, t);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_step(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:reactor_step", &n, &t))
        return NULL;
    return Py_BuildValue("d",reactor_step(n, t));
}

static PyObject*
py_reactor_time(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_time", &n))
        return NULL;
    double t = reactor_time(n);
    return Py_BuildValue("d",t);
}

static PyObject*
py_reactor_mass(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_mass", &n))
        return NULL;
    double m = reactor_mass(n);
    return Py_BuildValue("d",m);
}

static PyObject*
py_reactor_volume(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_volume", &n))
        return NULL;
    double v = reactor_volume(n);
    return Py_BuildValue("d",v);
}

static PyObject*
py_reactor_density(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_density", &n))
        return NULL;
    double rho = reactor_density(n);
    return Py_BuildValue("d",rho);
}

static PyObject*
py_reactor_temperature(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_temperature", &n))
        return NULL;
    double t = reactor_temperature(n);
    return Py_BuildValue("d",t);
}

static PyObject*
py_reactor_enthalpy_mass(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_enthalpy_mass", &n))
        return NULL;
    double h = reactor_enthalpy_mass(n);
    return Py_BuildValue("d",h);
}

static PyObject*
py_reactor_intEnergy_mass(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_intEnergy_mass", &n))
        return NULL;
    double u = reactor_intEnergy_mass(n);
    return Py_BuildValue("d",u);
}

static PyObject*
py_reactor_pressure(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_pressure", &n))
        return NULL;
    double p = reactor_pressure(n);
    return Py_BuildValue("d",p);
}

static PyObject*
py_reactor_massFraction(PyObject *self, PyObject *args)
{
    int n;
    int k;
    if (!PyArg_ParseTuple(args, "ii:reactor_massFraction", &n, &k))
        return NULL;
    double y = reactor_massFraction(n, k);
    return Py_BuildValue("d",y);
}

// static PyObject*
// py_reactor_setArea(PyObject *self, PyObject *args)
// {
//     int n;
//     double a;
//     if (!PyArg_ParseTuple(args, "id:reactor_setArea", &n, &a))
//         return NULL;
//     int iok = reactor_setArea(n,a);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_reactor_setExtTemp(PyObject *self, PyObject *args)
// {
//     int n;
//     double t;
//     if (!PyArg_ParseTuple(args, "id:reactor_setExtTemp", &n, &t))
//         return NULL;
//     int iok = reactor_setExtTemp(n, t);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_reactor_setExtRadTemp(PyObject *self, PyObject *args)
// {
//     int n;
//     double t;
//     if (!PyArg_ParseTuple(args, "id:reactor_setExtRadTemp", &n, &t))
//         return NULL;
//     int iok = reactor_setExtRadTemp(n, t);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_reactor_setVDotCoeff(PyObject *self, PyObject *args)
// {
//     int n;
//     double vdt;
//     if (!PyArg_ParseTuple(args, "id:reactor_setVDotCoeff", &n, &vdt))
//         return NULL;
//     int iok = reactor_setVDotCoeff(n, vdt);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_reactor_setHeatTransferCoeff(PyObject *self, PyObject *args)
// {
//     int n;
//     double h;
//     if (!PyArg_ParseTuple(args, "id:reactor_setHeatTransferCoeff", &n, &h))
//         return NULL;
//     int iok = reactor_setHeatTransferCoeff(n, h);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_reactor_setEmissivity(PyObject *self, PyObject *args)
// {
//     int n;
//     double e;
//     if (!PyArg_ParseTuple(args, "id:reactor_setEmissivity", &n, &e))
//         return NULL;
//     int iok = reactor_setEmissivity(n, e);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_reactor_setExtPressure(PyObject *self, PyObject *args)
// {
//     int n;
//     double p;
//     if (!PyArg_ParseTuple(args, "id:reactor_setExtPressure", &n, &p))
//         return NULL;
//     int iok = reactor_setExtPressure(n,p);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

static PyObject*
py_flowdev_new(PyObject *self, PyObject *args)
{
    int type;
    if (!PyArg_ParseTuple(args, "i:flowdev_new", &type))
        return NULL;
    int n = flowdev_new(type);
    return Py_BuildValue("i",n);
}

static PyObject*
py_flowdev_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flowdev_del", &n))
        return NULL;
    int iok = flowdev_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_install(PyObject *self, PyObject *args)
{
    int n, r1, r2;
    if (!PyArg_ParseTuple(args, "iii:flowdev_install", &n, &r1, &r2))
        return NULL;
    int iok = flowdev_install(n, r1, r2);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_massFlowRate(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flowdev_massFlowRate", &n))
        return NULL;
    double mdot = flowdev_massFlowRate(n);
    return Py_BuildValue("d",mdot);
}

static PyObject*
py_flowdev_setpoint(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flowdev_setpoint", &n))
        return NULL;
    double v = flowdev_setpoint(n);
    return Py_BuildValue("d",v);
}

static PyObject*
py_flowdev_setSetpoint(PyObject *self, PyObject *args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:flowdev_setSetpoint", &n, &v))
        return NULL;
    int iok = flowdev_setSetpoint(n, v);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

// static PyObject*
// py_flowdev_setGains(PyObject *self, PyObject *args)
// {
//     int n, sz;
//     PyObject* gains;
//     if (!PyArg_ParseTuple(args, "iiO:flowdev_setGains", &n, &sz, &gains))
//         return NULL;
//     double* x = (double*)((PyArrayObject*)gains)->data;
//     int iok = flowdev_setGains(n, sz, x);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

static PyObject*
py_flowdev_setParameters(PyObject *self, PyObject *args)
{
    int n, sz;
    PyObject* c;
    if (!PyArg_ParseTuple(args, "iiO:flowdev_setParameters", &n, &sz, &c))
        return NULL;
    double* x = (double*)((PyArrayObject*)c)->data;
    int iok = flowdev_setParameters(n, sz, x);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_setFunction(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:flowdev_setFunction", &n, &m))
        return NULL;
    int iok = flowdev_setFunction(n, m);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

// static PyObject*
// py_flowdev_getGains(PyObject *self, PyObject *args)
// {
//     int n, sz;
//     if (!PyArg_ParseTuple(args, "ii:flowdev_getGains", &n, &sz))
//         return NULL;
//     PyArrayObject* x = 
//         (PyArrayObject*)PyArray_FromDims(1, &sz, PyArray_DOUBLE);
//     double* xd = (double*)x->data;
//     int iok = flowdev_getGains(n, sz, xd);
//     if (iok < 0) return reportError(iok);
//     return PyArray_Return(x);
// }

// static PyObject*
// py_flowdev_reset(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:flowdev_reset", &n))
//         return NULL;
//     int iok = flowdev_reset(n);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_flowdev_update(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:flowdev_update", &n))
//         return NULL;
//     int iok = flowdev_update(n);
//     if (iok < 0) return reportError(iok);
//     return Py_BuildValue("i",0);
// }

// static PyObject*
// py_flowdev_maxError(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:flowdev_maxError", &n))
//         return NULL;
//     double e = flowdev_maxError(n);
//     return Py_BuildValue("d",e);
// }

static PyObject*
py_flowdev_ready(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flowdev_ready", &n))
        return NULL;
    int iok = flowdev_ready(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",iok);
}


static PyObject*
py_wall_new(PyObject *self, PyObject *args)
{
    int type;
    if (!PyArg_ParseTuple(args, "i:wall_new", &type))
        return NULL;
    int n = wall_new(type);
    return Py_BuildValue("i",n);
}

static PyObject*
py_wall_del(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:wall_del", &n))
        return NULL;
    int iok = wall_del(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_install(PyObject *self, PyObject *args)
{
    int n, r1, r2;
    if (!PyArg_ParseTuple(args, "iii:wall_install", &n, &r1, &r2))
        return NULL;
    int iok = wall_install(n, r1, r2);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_vdot(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:wall_vdot", &n, &t))
        return NULL;
    double vdt = wall_vdot(n,t);
    return Py_BuildValue("d",vdt);
}

static PyObject*
py_wall_Q(PyObject *self, PyObject *args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:wall_Q", &n, &t))
        return NULL;
    return Py_BuildValue("d",wall_Q(n, t));
}

static PyObject*
py_wall_area(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:wall_area", &n))
        return NULL;
    return Py_BuildValue("d",wall_area(n));
}

static PyObject*
py_wall_setArea(PyObject *self, PyObject *args)
{
    int n;
    double area;
    if (!PyArg_ParseTuple(args, "id:wall_setArea", &n, &area))
        return NULL;
    int iok = wall_setArea(n, area);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setThermalResistance(PyObject *self, PyObject *args)
{
    int n;
    double rth;
    if (!PyArg_ParseTuple(args, "id:wall_setThermalResistance", &n, &rth))
        return NULL;
    int iok = wall_setThermalResistance(n,rth);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setHeatTransferCoeff(PyObject *self, PyObject *args)
{
    int n;
    double u;
    if (!PyArg_ParseTuple(args, "id:wall_setHeatTransferCoeff", &n, &u))
        return NULL;
    int iok = wall_setHeatTransferCoeff(n,u);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setExpansionRateCoeff(PyObject *self, PyObject *args)
{
    int n;
    double k;
    if (!PyArg_ParseTuple(args, "id:wall_setExpansionRateCoeff", &n, &k))
        return NULL;
    int iok = wall_setExpansionRateCoeff(n,k);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setExpansionRate(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:wall_setExpansionRate", &n, &m))
        return NULL;
    int iok = wall_setExpansionRate(n,m);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setHeatFlux(PyObject *self, PyObject *args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:wall_setHeatFlux", &n, &m))
        return NULL;
    int iok = wall_setHeatFlux(n,m);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_ready(PyObject *self, PyObject *args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:wall_ready", &n))
        return NULL;
    int iok = wall_ready(n);
    if (iok < 0) return reportError(iok);
    return Py_BuildValue("i",0);
}


static PyMethodDef ct_methods[] = {
    //    {"reactor_setExtTemp", py_reactor_setExtTemp, METH_VARARGS},
    {"flowdev_ready", py_flowdev_ready, METH_VARARGS},
    //{"reactor_setEmissivity", py_reactor_setEmissivity, METH_VARARGS},
    {"reactor_setInitialTime", py_reactor_setInitialTime, METH_VARARGS},
    {"flowdev_new", py_flowdev_new, METH_VARARGS},
    {"flowdev_massFlowRate", py_flowdev_massFlowRate, METH_VARARGS},
    {"flowdev_del", py_flowdev_del, METH_VARARGS},
    {"flowdev_setpoint", py_flowdev_setpoint, METH_VARARGS},
    {"reactor_temperature", py_reactor_temperature, METH_VARARGS},
    {"flowdev_setSetpoint", py_flowdev_setSetpoint, METH_VARARGS},
    //{"flowdev_reset", py_flowdev_reset, METH_VARARGS},
    //{"reactor_setExtRadTemp", py_reactor_setExtRadTemp, METH_VARARGS},
    {"flowdev_install", py_flowdev_install, METH_VARARGS},
    {"reactor_setThermoMgr", py_reactor_setThermoMgr, METH_VARARGS},
    {"reactor_setEnergy", py_reactor_setEnergy, METH_VARARGS},
    {"reactor_volume", py_reactor_volume, METH_VARARGS},
    {"reactor_time", py_reactor_time, METH_VARARGS},
    {"reactor_advance", py_reactor_advance, METH_VARARGS},
    {"reactor_step", py_reactor_step, METH_VARARGS},
    //{"reactor_setExtPressure", py_reactor_setExtPressure, METH_VARARGS},
    //{"flowdev_setGains", py_flowdev_setGains, METH_VARARGS},
    {"flowdev_setParameters", py_flowdev_setParameters, METH_VARARGS},
    {"flowdev_setFunction", py_flowdev_setFunction, METH_VARARGS},
    {"reactor_mass", py_reactor_mass, METH_VARARGS},
    {"reactor_new", py_reactor_new, METH_VARARGS},
    //{"reactor_setVDotCoeff", py_reactor_setVDotCoeff, METH_VARARGS},
    //{"reactor_setHeatTransferCoeff", py_reactor_setHeatTransferCoeff, METH_VARARGS},
    {"reactor_enthalpy_mass", py_reactor_enthalpy_mass, METH_VARARGS},
    //{"flowdev_maxError", py_flowdev_maxError, METH_VARARGS},
    //{"flowdev_getGains", py_flowdev_getGains, METH_VARARGS},
    //{"flowdev_update", py_flowdev_update, METH_VARARGS},
    //{"reactor_setArea", py_reactor_setArea, METH_VARARGS},
    {"reactor_pressure", py_reactor_pressure, METH_VARARGS},
    {"reactor_setInitialVolume", py_reactor_setInitialVolume, METH_VARARGS},
    {"reactor_density", py_reactor_density, METH_VARARGS},
    {"reactor_setKineticsMgr", py_reactor_setKineticsMgr, METH_VARARGS},
    {"reactor_del", py_reactor_del, METH_VARARGS},
    {"reactor_intEnergy_mass", py_reactor_intEnergy_mass, METH_VARARGS},
    {"reactor_massFraction", py_reactor_massFraction, METH_VARARGS},
    {"wall_install", py_wall_install, METH_VARARGS},
    {"wall_area", py_wall_area, METH_VARARGS},
    {"wall_setArea", py_wall_setArea, METH_VARARGS},
    {"wall_setThermalResistance", py_wall_setThermalResistance, METH_VARARGS},
    {"wall_setHeatTransferCoeff", py_wall_setHeatTransferCoeff, METH_VARARGS},
    {"wall_setHeatFlux", py_wall_setHeatFlux, METH_VARARGS},
    {"wall_Q", py_wall_Q, METH_VARARGS},
    {"wall_new", py_wall_new, METH_VARARGS},
    {"wall_vdot", py_wall_vdot, METH_VARARGS},
    {"wall_del", py_wall_del, METH_VARARGS},
    {"wall_setExpansionRate", py_wall_setExpansionRate, METH_VARARGS},
    {"wall_setExpansionRateCoeff", py_wall_setExpansionRateCoeff, METH_VARARGS},
    {"wall_ready", py_wall_ready, METH_VARARGS},
    {NULL, NULL}
};

extern "C" {

    /* Initialization function for the module */

    DL_EXPORT(void) initctreactor(void)
    {
        PyObject *m, *d;

        /* Initialize the type of the new type object here; doing it here
         * is required for portability to Windows without requiring C++. */
 
        /* Create the module and add the functions */
        m = Py_InitModule("ctreactor", ct_methods);
        import_array();
 
        /* Add some symbolic constants to the module */
        d = PyModule_GetDict(m);
        ErrorObject = PyErr_NewException("cantera.error", NULL, NULL);
        PyDict_SetItemString(d, "error", ErrorObject);
    }

}

