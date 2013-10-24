
static PyObject*
py_mix_new(PyObject* self, PyObject* args)
{
    int _val;
    _val = mix_new();
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_del(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_del", &i)) {
        return NULL;
    }

    _val = mix_del(i);
    if (int(_val) < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_addPhase(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int j;
    double moles;
    if (!PyArg_ParseTuple(args, "iid:mix_addPhase", &i, &j, &moles)) {
        return NULL;
    }

    _val = mix_addPhase(i,j,moles);
    if (int(_val) < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_init(PyObject* self, PyObject* args)
{
    int i;
    int _val;
    if (!PyArg_ParseTuple(args, "i:mix_init", &i)) {
        return NULL;
    }

    _val = mix_init(i);
    if (_val < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_nElements(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_nElements", &i)) {
        return NULL;
    }

    _val = int(mix_nElements(i));
    if (_val <= -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_mix_elementIndex(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* name;
    if (!PyArg_ParseTuple(args, "is:mix_elementIndex", &i, &name)) {
        return NULL;
    }

    _val = int(mix_elementIndex(i,name));
    if (_val <= -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_mix_nSpecies(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_nSpecies", &i)) {
        return NULL;
    }

    _val = int(mix_nSpecies(i));
    if (_val <= -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_mix_speciesIndex(PyObject* self, PyObject* args)
{
    size_t _val;
    int i, k, p;
    if (!PyArg_ParseTuple(args, "iii:mix_speciesIndex", &i, &k, &p)) {
        return NULL;
    }

    _val = mix_speciesIndex(i,k,p);
    if (int(_val) <= -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_mix_nAtoms(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int k;
    int m;
    if (!PyArg_ParseTuple(args, "iii:mix_nAtoms", &i, &k, &m)) {
        return NULL;
    }

    _val = mix_nAtoms(i,k,m);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}

static PyObject*
py_mix_setTemperature(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double t;
    if (!PyArg_ParseTuple(args, "id:mix_setTemperature", &i, &t)) {
        return NULL;
    }

    _val = mix_setTemperature(i,t);
    if (int(_val) <= -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_minTemp(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_minTemp", &i)) {
        return NULL;
    }

    _val = mix_minTemp(i);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_maxTemp(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_maxTemp", &i)) {
        return NULL;
    }

    _val = mix_maxTemp(i);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_charge(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_charge", &i)) {
        return NULL;
    }

    _val = mix_charge(i);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_phaseCharge(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int p;
    if (!PyArg_ParseTuple(args, "ii:mix_phaseCharge", &i, &p)) {
        return NULL;
    }

    _val = mix_phaseCharge(i,p);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_temperature(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_temperature", &i)) {
        return NULL;
    }

    _val = mix_temperature(i);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_setPressure(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    double p;
    if (!PyArg_ParseTuple(args, "id:mix_setPressure", &i, &p)) {
        return NULL;
    }

    _val = mix_setPressure(i,p);
    if (int(_val) <= -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_pressure(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    if (!PyArg_ParseTuple(args, "i:mix_pressure", &i)) {
        return NULL;
    }

    _val = mix_pressure(i);
    if (int(_val) <= -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}

static PyObject*
py_mix_phaseMoles(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int n;
    if (!PyArg_ParseTuple(args, "ii:mix_phaseMoles", &i, &n)) {
        return NULL;
    }

    _val = mix_phaseMoles(i,n);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_setPhaseMoles(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "iid:mix_setPhaseMoles", &i, &n, &v)) {
        return NULL;
    }

    _val = mix_setPhaseMoles(i,n,v);
    if (int(_val) < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_speciesMoles(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int k;
    if (!PyArg_ParseTuple(args, "ii:mix_speciesMoles", &i, &k)) {
        return NULL;
    }

    _val = mix_speciesMoles(i,k);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_elementMoles(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    int m;
    if (!PyArg_ParseTuple(args, "ii:mix_elementMoles", &i, &m)) {
        return NULL;
    }

    _val = mix_elementMoles(i,m);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_setMoles(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    PyObject* n;
    if (!PyArg_ParseTuple(args, "iO:mix_setMoles", &i, &n)) {
        return NULL;
    }


    PyArrayObject* n_array = (PyArrayObject*)n;
    double* n_data = (double*)n_array->data;
    size_t n_len = n_array->dimensions[0];

    _val = mix_setMoles(i,n_len,n_data);
    if (int(_val) < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}


static PyObject*
py_mix_setMolesByName(PyObject* self, PyObject* args)
{
    int _val;
    int i;
    char* n;
    if (!PyArg_ParseTuple(args, "is:mix_setMolesByName", &i, &n)) {
        return NULL;
    }

    _val = mix_setMolesByName(i,n);
    if (int(_val) < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

static PyObject*
py_mix_equilibrate(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    char* XY;
    double err;
    int maxsteps, maxiter, loglevel;
    if (!PyArg_ParseTuple(args, "isdiii:mix_equilibrate", &i, &XY, &err,
                          &maxsteps, &maxiter, &loglevel)) {
        return NULL;
    }

    _val = mix_equilibrate(i,XY,err,maxsteps,maxiter,loglevel);
    if (int(_val) == DERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d",_val);
}


static PyObject*
py_mix_vcs_equilibrate(PyObject* self, PyObject* args)
{
    double _val;
    int i;
    char* XY;
    int estimateEquil;
    int printLvl;
    int solver;
    double rtol;
    int maxsteps;
    int maxiter;
    int loglevel;
    if (!PyArg_ParseTuple(args, "isiiidiii:mix_vcs_equilibrate", &i, &XY,
                          &estimateEquil, &printLvl, &solver,
                          &rtol, &maxsteps, &maxiter, &loglevel)) {
        return NULL;
    }

    _val = mix_vcs_equilibrate(i, XY, estimateEquil, printLvl, solver,
                               rtol, maxsteps, maxiter, loglevel);

    if (int(_val) == ERR) {
        return reportCanteraError();
    }
    return Py_BuildValue("d", _val);
}


static PyObject*
py_mix_getChemPotentials(PyObject* self, PyObject* args)
{
    int i;
    int _val;
    PyObject* mu;
    if (!PyArg_ParseTuple(args, "iO:mix_getChemPotentials", &i, &mu)) {
        return NULL;
    }

    PyArrayObject* mu_array = (PyArrayObject*)mu;
    double* mu_data = (double*)mu_array->data;
    size_t mu_len = mu_array->dimensions[0];

    _val = mix_getChemPotentials(i, mu_len, mu_data);
    if (int(_val) < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",_val);
}

