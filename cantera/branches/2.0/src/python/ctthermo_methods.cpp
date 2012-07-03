
static PyObject*
ct_newThermoFromXML(PyObject* self, PyObject* args)
{
    int mxml;
    //char* id;
    if (!PyArg_ParseTuple(args, "i:ct_newThermoFromXML", &mxml)) {
        return NULL;
    }
    int n = int(newThermoFromXML(mxml));
    if (n < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",n);
}

static PyObject*
thermo_delete(PyObject* self, PyObject* args)
{
    int th;
    if (!PyArg_ParseTuple(args, "i:thermo_delete", &th)) {
        return NULL;
    }
    delThermo(th);
    return Py_BuildValue("i",0);
}

static PyObject*
thermo_refpressure(PyObject* self, PyObject* args)
{
    int th;
    if (!PyArg_ParseTuple(args, "i:refpressure", &th)) {
        return NULL;
    }
    return Py_BuildValue("d",th_refPressure(th));
}

static PyObject*
thermo_mintemp(PyObject* self, PyObject* args)
{
    int th, k;
    if (!PyArg_ParseTuple(args, "ii:mintemp", &th, &k)) {
        return NULL;
    }
    return Py_BuildValue("d",th_minTemp(th,k));
}

static PyObject*
thermo_maxtemp(PyObject* self, PyObject* args)
{
    int th, k;
    if (!PyArg_ParseTuple(args, "ii:maxtemp", &th, &k)) {
        return NULL;
    }
    return Py_BuildValue("d",th_maxTemp(th,k));
}

static PyObject*
thermo_import(PyObject* self, PyObject* args)
{
    int n, mxml;
    char* id;
    if (!PyArg_ParseTuple(args, "iis:import", &n, &mxml, &id)) {
        return NULL;
    }
    int iok = import_phase(n, mxml, id);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}


static PyObject*
thermo_getfp(PyObject* self, PyObject* args)
{
    double vv = -999.999;
    bool ok = true;
    int th;
    int job;

    if (!PyArg_ParseTuple(args, "ii:thermo_getfp", &th, &job)) {
        return NULL;
    }

    // try {

    // floating-point attributes
    switch (job) {
    case 1:
        vv = th_enthalpy_mole(th);
        break;
    case 2:
        vv = th_intEnergy_mole(th);
        break;
    case 3:
        vv = th_entropy_mole(th);
        break;
    case 4:
        vv = th_gibbs_mole(th);
        break;
    case 5:
        vv = th_cp_mole(th);
        break;
    case 6:
        vv = th_cv_mole(th);
        break;
    case 7:
        vv = th_pressure(th);
        break;
    case 8:
        vv = th_enthalpy_mass(th);
        break;
    case 9:
        vv = th_intEnergy_mass(th);
        break;
    case 10:
        vv = th_entropy_mass(th);
        break;
    case 11:
        vv = th_gibbs_mass(th);
        break;
    case 12:
        vv = th_cp_mass(th);
        break;
    case 13:
        vv = th_cv_mass(th);
        break;
    case 25:
        vv = th_electricPotential(th);
        break;
    case 50:
        vv = th_critTemperature(th);
        break;
    case 51:
        vv = th_critPressure(th);
        break;
    case 52:
        vv = th_critDensity(th);
        break;
    case 53:
        vv = th_vaporFraction(th);
        break;

    default:
        ok = false;
    }
    if (ok) {
        if (vv == -999.999) {
            return reportCanteraError();
        }
        return Py_BuildValue("d",vv);
    } else {
        PyErr_SetString(ErrorObject,"Unknown floating-point attribute");
        return NULL;
    }
    //}
    //catch (CanteraError) {
    // return reportCanteraError();
    //}
}

static PyObject*
thermo_setfp(PyObject* self, PyObject* args)
{
    double v1 = -1.0, v2 = -1.0;
    int iok = -2;
    int th;
    int job;

    if (!PyArg_ParseTuple(args, "iidd:thermo_setfp", &th, &job, &v1, &v2)) {
        return NULL;
    }

    //vector_fp v(2);
    double v[2];
    v[0] = v1;
    v[1] = v2;
    // set floating-point attributes
    switch (job) {
    case 1:
        iok = th_setPressure(th, v1);
        break;
    case 2:
        iok = th_set_HP(th, v);
        break;
    case 3:
        iok = th_set_UV(th, v);
        break;
    case 4:
        iok = th_set_SV(th, v);
        break;
    case 5:
        iok = th_set_SP(th, v);
        break;
    case 6:
        iok = th_setElectricPotential(th, v[0]);
        break;
    case 7:
        iok = th_setState_Tsat(th, v1, v2);
        break;
    case 8:
        iok = th_setState_Psat(th, v1, v2);
        break;
    default:
        iok = -10;
    }
    //delete v;
    if (iok >= 0) {
        return Py_BuildValue("i",iok);
    }
    if (iok == -1) {
        return reportCanteraError();
    } else {
        PyErr_SetString(ErrorObject,"Error in thermo_setfp");
        return NULL;
    }
}


static PyObject*
thermo_getarray(PyObject* self, PyObject* args)
{
    int th;
    int job;

    if (!PyArg_ParseTuple(args, "ii:thermo_getarray", &th, &job)) {
        return NULL;
    }

    size_t nsp = th_nSpecies(th);
    size_t nel = phase_nElements(th);
    size_t xlen = (job == 21 ? nel : nsp);

    // array attributes
    int iok = -22;

#ifdef HAS_NUMPY
    npy_intp nnn = xlen;
    PyArrayObject* x =
        (PyArrayObject*)PyArray_SimpleNew(1, &nnn, PyArray_DOUBLE);
#else
    int nnn = int(xlen);
    PyArrayObject* x =
        (PyArrayObject*)PyArray_FromDims(1, &nnn, PyArray_DOUBLE);
#endif
    double* xd = (double*)x->data;
    switch (job) {
    case 20:
        iok = th_chemPotentials(th,nsp,xd);
        break;
    case 21:
        iok = th_elementPotentials(th,nel,xd);
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
    } else if (iok == -1) {
        return reportCanteraError();
    } else {
        PyErr_SetString(ErrorObject,"Unknown array attribute");
        return NULL;
    }
}


static PyObject*
thermo_equil(PyObject* self, PyObject* args)
{
    int iok = -2;
    int th;
    char* XY;
    int solver;
    double rtol;
    int maxsteps;
    int maxiter;
    int loglevel;

    if (!PyArg_ParseTuple(args, "isidiii:thermo_equil", &th, &XY,
                          &solver, &rtol, &maxsteps, &maxiter, &loglevel)) {
        return NULL;
    }

    iok = th_equil(th, XY, solver, rtol, maxsteps, maxiter, loglevel);
    if (iok >= 0) {
        return Py_BuildValue("i",iok);
    }
    if (iok == -1) {
        return reportCanteraError();
    } else {
        PyErr_SetString(ErrorObject,"Error in thermo_equil");
        return NULL;
    }
}





