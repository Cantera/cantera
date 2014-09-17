
static PyObject*
py_reactor_new(PyObject* self, PyObject* args)
{
    int type;
    if (!PyArg_ParseTuple(args, "i:reactor_new", &type)) {
        return NULL;
    }
    int n = reactor_new(type);
    return Py_BuildValue("i",n);
}

static PyObject*
py_reactor_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_del", &n)) {
        return NULL;
    }
    int iok = reactor_del(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setInitialVolume(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:reactor_setInitialVolume", &n, &v)) {
        return NULL;
    }
    int iok = reactor_setInitialVolume(n,v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setEnergy(PyObject* self, PyObject* args)
{
    int n, eflag;
    if (!PyArg_ParseTuple(args, "ii:reactor_setEnergy", &n, &eflag)) {
        return NULL;
    }
    int iok = reactor_setEnergy(n, eflag);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setThermoMgr(PyObject* self, PyObject* args)
{
    int n;
    int th;
    if (!PyArg_ParseTuple(args, "ii:reactor_setThermoMgr", &n, &th)) {
        return NULL;
    }
    int iok = reactor_setThermoMgr(n, th);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_setKineticsMgr(PyObject* self, PyObject* args)
{
    int n;
    int kin;
    if (!PyArg_ParseTuple(args, "ii:reactor_setKineticsMgr", &n, &kin)) {
        return NULL;
    }
    int iok = reactor_setKineticsMgr(n, kin);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactor_nSensParams(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:reactor_nSensParams", &i)) {
        return NULL;
    }

    _val = int(reactor_nSensParams(i));
    return Py_BuildValue("i",_val);
}


static PyObject*
py_reactor_addSensitivityReaction(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int rxn;
    if (!PyArg_ParseTuple(args, "ii:reactor_addSensitivityReaction", &i, &rxn)) {
        return NULL;
    }

    _val = reactor_addSensitivityReaction(i,rxn);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_flowReactor_setMassFlowRate(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double mdot;
    if (!PyArg_ParseTuple(args, "id:flowReactor_setMassFlowRate", &i, &mdot)) {
        return NULL;
    }

    _val = flowReactor_setMassFlowRate(i,mdot);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_reactor_mass(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_mass", &n)) {
        return NULL;
    }
    double m = reactor_mass(n);
    return Py_BuildValue("d",m);
}

static PyObject*
py_reactor_volume(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_volume", &n)) {
        return NULL;
    }
    double v = reactor_volume(n);
    return Py_BuildValue("d",v);
}

static PyObject*
py_reactor_density(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_density", &n)) {
        return NULL;
    }
    double rho = reactor_density(n);
    return Py_BuildValue("d",rho);
}

static PyObject*
py_reactor_temperature(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_temperature", &n)) {
        return NULL;
    }
    double t = reactor_temperature(n);
    return Py_BuildValue("d",t);
}

static PyObject*
py_reactor_enthalpy_mass(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_enthalpy_mass", &n)) {
        return NULL;
    }
    double h = reactor_enthalpy_mass(n);
    return Py_BuildValue("d",h);
}

static PyObject*
py_reactor_intEnergy_mass(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_intEnergy_mass", &n)) {
        return NULL;
    }
    double u = reactor_intEnergy_mass(n);
    return Py_BuildValue("d",u);
}

static PyObject*
py_reactor_pressure(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactor_pressure", &n)) {
        return NULL;
    }
    double p = reactor_pressure(n);
    return Py_BuildValue("d",p);
}

static PyObject*
py_reactor_massFraction(PyObject* self, PyObject* args)
{
    int n;
    int k;
    if (!PyArg_ParseTuple(args, "ii:reactor_massFraction", &n, &k)) {
        return NULL;
    }
    double y = reactor_massFraction(n, k);
    return Py_BuildValue("d",y);
}


static PyObject*
py_flowdev_new(PyObject* self, PyObject* args)
{
    int type;
    if (!PyArg_ParseTuple(args, "i:flowdev_new", &type)) {
        return NULL;
    }
    int n = flowdev_new(type);
    return Py_BuildValue("i",n);
}

static PyObject*
py_flowdev_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flowdev_del", &n)) {
        return NULL;
    }
    int iok = flowdev_del(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_install(PyObject* self, PyObject* args)
{
    int n, r1, r2;
    if (!PyArg_ParseTuple(args, "iii:flowdev_install", &n, &r1, &r2)) {
        return NULL;
    }
    int iok = flowdev_install(n, r1, r2);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_setMaster(PyObject* self, PyObject* args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:flowdev_setMaster", &n, &m)) {
        return NULL;
    }
    int iok = flowdev_setMaster(n, m);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_massFlowRate(PyObject* self, PyObject* args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:flowdev_massFlowRate", &n, &t)) {
        return NULL;
    }
    double mdot = flowdev_massFlowRate(n, t);
    return Py_BuildValue("d",mdot);
}

// static PyObject*
// py_flowdev_setpoint(PyObject *self, PyObject *args)
// {
//     int n;
//     if (!PyArg_ParseTuple(args, "i:flowdev_setpoint", &n))
//         return NULL;
//     double v = flowdev_setpoint(n);
//     return Py_BuildValue("d",v);
// }

static PyObject*
py_flowdev_setMassFlowRate(PyObject* self, PyObject* args)
{
    int n;
    double mdot;
    if (!PyArg_ParseTuple(args, "id:flowdev_setMassFlowRate", &n, &mdot)) {
        return NULL;
    }
    int iok = flowdev_setMassFlowRate(n, mdot);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_setParameters(PyObject* self, PyObject* args)
{
    int n, sz;
    PyObject* c;
    if (!PyArg_ParseTuple(args, "iiO:flowdev_setParameters", &n, &sz, &c)) {
        return NULL;
    }
    PyArrayObject* ca = (PyArrayObject*)
                        PyArray_ContiguousFromObject(c, PyArray_DOUBLE, 1, 1);
    double* x = (double*)ca->data;
    int iok = flowdev_setParameters(n, sz, x);
    Py_DECREF(ca);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_setFunction(PyObject* self, PyObject* args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:flowdev_setFunction", &n, &m)) {
        return NULL;
    }
    int iok = flowdev_setFunction(n, m);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_flowdev_ready(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:flowdev_ready", &n)) {
        return NULL;
    }
    int iok = flowdev_ready(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",iok);
}


static PyObject*
py_wall_new(PyObject* self, PyObject* args)
{
    int type;
    if (!PyArg_ParseTuple(args, "i:wall_new", &type)) {
        return NULL;
    }
    int n = wall_new(type);
    return Py_BuildValue("i",n);
}

static PyObject*
py_wall_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:wall_del", &n)) {
        return NULL;
    }
    int iok = wall_del(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_install(PyObject* self, PyObject* args)
{
    int n, r1, r2;
    if (!PyArg_ParseTuple(args, "iii:wall_install", &n, &r1, &r2)) {
        return NULL;
    }
    int iok = wall_install(n, r1, r2);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setkinetics(PyObject* self, PyObject* args)
{
    int n, k1, k2;
    if (!PyArg_ParseTuple(args, "iii:wall_setkinetics", &n, &k1, &k2)) {
        return NULL;
    }
    int iok = wall_setkinetics(n, k1, k2);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_vdot(PyObject* self, PyObject* args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:wall_vdot", &n, &t)) {
        return NULL;
    }
    double vdt = wall_vdot(n,t);
    return Py_BuildValue("d",vdt);
}

static PyObject*
py_wall_Q(PyObject* self, PyObject* args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:wall_Q", &n, &t)) {
        return NULL;
    }
    return Py_BuildValue("d",wall_Q(n, t));
}

static PyObject*
py_wall_area(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:wall_area", &n)) {
        return NULL;
    }
    return Py_BuildValue("d",wall_area(n));
}

static PyObject*
py_wall_setArea(PyObject* self, PyObject* args)
{
    int n;
    double area;
    if (!PyArg_ParseTuple(args, "id:wall_setArea", &n, &area)) {
        return NULL;
    }
    int iok = wall_setArea(n, area);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setThermalResistance(PyObject* self, PyObject* args)
{
    int n;
    double rth;
    if (!PyArg_ParseTuple(args, "id:wall_setThermalResistance", &n, &rth)) {
        return NULL;
    }
    int iok = wall_setThermalResistance(n,rth);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setHeatTransferCoeff(PyObject* self, PyObject* args)
{
    int n;
    double u;
    if (!PyArg_ParseTuple(args, "id:wall_setHeatTransferCoeff", &n, &u)) {
        return NULL;
    }
    int iok = wall_setHeatTransferCoeff(n,u);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setEmissivity(PyObject* self, PyObject* args)
{
    int n;
    double epsilon;
    if (!PyArg_ParseTuple(args, "id:wall_setEmissivity", &n, &epsilon)) {
        return NULL;
    }
    int iok = wall_setEmissivity(n,epsilon);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setExpansionRateCoeff(PyObject* self, PyObject* args)
{
    int n;
    double k;
    if (!PyArg_ParseTuple(args, "id:wall_setExpansionRateCoeff", &n, &k)) {
        return NULL;
    }
    int iok = wall_setExpansionRateCoeff(n,k);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setVelocity(PyObject* self, PyObject* args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:wall_setVelocity", &n, &m)) {
        return NULL;
    }
    int iok = wall_setVelocity(n,m);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_setHeatFlux(PyObject* self, PyObject* args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:wall_setHeatFlux", &n, &m)) {
        return NULL;
    }
    int iok = wall_setHeatFlux(n,m);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_ready(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:wall_ready", &n)) {
        return NULL;
    }
    int iok = wall_ready(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_wall_addSensitivityReaction(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int lr;
    int rxn;
    if (!PyArg_ParseTuple(args, "iii:wall_addSensitivityReaction", &i, &lr, &rxn)) {
        return NULL;
    }

    _val = wall_addSensitivityReaction(i,lr,rxn);
    if (int(_val) == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_reactornet_new(PyObject* self, PyObject* args)
{
    int n = reactornet_new();
    return Py_BuildValue("i",n);
}

static PyObject*
py_reactornet_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactornet_del", &n)) {
        return NULL;
    }
    int iok = reactornet_del(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactornet_setTolerances(PyObject* self, PyObject* args)
{
    int n;
    double rtol, atol;
    if (!PyArg_ParseTuple(args, "idd:reactornet_setTolerances", &n, &rtol, &atol)) {
        return NULL;
    }
    int iok = reactornet_setTolerances(n, rtol, atol);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactornet_setSensitivityTolerances(PyObject* self, PyObject* args)
{
    int n;
    double rtol, atol;
    if (!PyArg_ParseTuple(args, "idd:reactornet_setSensitivityTolerances", &n, &rtol, &atol)) {
        return NULL;
    }
    int iok = reactornet_setSensitivityTolerances(n, rtol, atol);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactornet_setInitialTime(PyObject* self, PyObject* args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:reactornet_setInitialTime", &n, &t)) {
        return NULL;
    }
    int iok = reactornet_setInitialTime(n, t);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactornet_time(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:reactornet_time", &n)) {
        return NULL;
    }
    double t = reactornet_time(n);
    return Py_BuildValue("d",t);
}

static PyObject*
py_reactornet_addreactor(PyObject* self, PyObject* args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:reactornet_addreactor", &n, &m)) {
        return NULL;
    }
    int iok = reactornet_addreactor(n, m);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactornet_advance(PyObject* self, PyObject* args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:reactornet_advance", &n, &t)) {
        return NULL;
    }
    int iok = reactornet_advance(n, t);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_reactornet_step(PyObject* self, PyObject* args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:reactornet_step", &n, &t)) {
        return NULL;
    }
    double ret = reactornet_step(n, t);
    if (ret == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d", ret);
}


static PyObject*
py_reactornet_rtol(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:reactornet_rtol", &i)) {
        return NULL;
    }

    _val = reactornet_rtol(i);
    return Py_BuildValue("d",_val);
}


static PyObject*
py_reactornet_atol(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:reactornet_atol", &i)) {
        return NULL;
    }

    _val = reactornet_atol(i);
    return Py_BuildValue("d",_val);
}


static PyObject*
py_reactornet_sensitivity(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    char* v;
    int p, r;
    if (!PyArg_ParseTuple(args, "isii:reactornet_sensitivity", &i, &v, &p, &r)) {
        return NULL;
    }

    _val = reactornet_sensitivity(i,v,p,r);
    return Py_BuildValue("d",_val);
}

