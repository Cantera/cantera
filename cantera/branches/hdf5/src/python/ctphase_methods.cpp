
static PyObject*
py_temperature(PyObject* self, PyObject* args)
{
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_temperature", &ph)) {
        return NULL;
    }
    return Py_BuildValue("d",phase_temperature(ph));
}

static PyObject*
py_density(PyObject* self, PyObject* args)
{
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_density", &ph)) {
        return NULL;
    }
    return Py_BuildValue("d",phase_density(ph));
}

static PyObject*
py_molardensity(PyObject* self, PyObject* args)
{
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_molardensity", &ph)) {
        return NULL;
    }
    return Py_BuildValue("d",phase_molarDensity(ph));
}

static PyObject*
py_meanmolwt(PyObject* self, PyObject* args)
{
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_meanmolwt", &ph)) {
        return NULL;
    }
    return Py_BuildValue("d",phase_meanMolecularWeight(ph));
}

static PyObject*
py_molefraction(PyObject* self, PyObject* args)
{
    int ph, k;
    if (!PyArg_ParseTuple(args, "ii:py_molefraction", &ph, &k)) {
        return NULL;
    }
    return Py_BuildValue("d",phase_moleFraction(ph, k));
}

static PyObject*
py_massfraction(PyObject* self, PyObject* args)
{
    int ph, k;
    if (!PyArg_ParseTuple(args, "ii:py_massfraction", &ph, &k)) {
        return NULL;
    }
    return Py_BuildValue("d",phase_massFraction(ph, k));
}

static PyObject*
py_nelements(PyObject* self, PyObject* args)
{
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_nelements", &ph)) {
        return NULL;
    }
    return Py_BuildValue("i",phase_nElements(ph));
}

static PyObject*
py_nspecies(PyObject* self, PyObject* args)
{
    int ph;
    if (!PyArg_ParseTuple(args, "i:py_nspecies", &ph)) {
        return NULL;
    }
    return Py_BuildValue("i",phase_nSpecies(ph));
}

static PyObject*
py_natoms(PyObject* self, PyObject* args)
{
    int ph, k, m;
    if (!PyArg_ParseTuple(args, "iii:py_natoms", &ph, &k, &m)) {
        return NULL;
    }
    return Py_BuildValue("d",phase_nAtoms(ph, k, m));
}

// static PyObject*
// py_addelement(PyObject *self, PyObject *args) {
//     int ph;
//     char* name;
//     double wt;
//     if (!PyArg_ParseTuple(args, "isd:py_addelement", &ph, &name, &wt))
//         return NULL;
//     int ok = phase_addElement(ph, name, wt);
//     if (ok < 0) return reportError(ok);
//     else return Py_BuildValue("i",0);
// }

static PyObject*
py_elementindex(PyObject* self, PyObject* args)
{
    int ph;
    char* nm;
    if (!PyArg_ParseTuple(args, "is:py_elementindex", &ph, &nm)) {
        return NULL;
    }
    size_t k = phase_elementIndex(ph,nm);
    return Py_BuildValue("i",k);
}

static PyObject*
py_speciesindex(PyObject* self, PyObject* args)
{
    int ph;
    char* nm;
    if (!PyArg_ParseTuple(args, "is:py_speciesindex", &ph, &nm)) {
        return NULL;
    }
    size_t k = phase_speciesIndex(ph,nm);
    return Py_BuildValue("i",k);
}

static PyObject*
py_report(PyObject* self, PyObject* args)
{
    int th, show_thermo;
    int buflen = 400;
    char* output_buf = new char[buflen];
    if (!PyArg_ParseTuple(args, "ii:py_report", &th, &show_thermo)) {
        return NULL;
    }
    int iok = phase_report(th, buflen, output_buf, show_thermo);
    if (iok < -1 && iok != -999) {
        delete output_buf;
        output_buf = new char[-iok];
        iok = phase_report(th, -iok, output_buf, show_thermo);
    }
    if (iok < 0) {
        return reportError(iok);
    }
    PyObject* s = Py_BuildValue("s",output_buf);
    delete output_buf;
    return s;
}


static PyObject*
phase_getarray(PyObject* self, PyObject* args)
{
    int ph;
    int job;

    if (!PyArg_ParseTuple(args, "ii:phase_getarray", &ph, &job)) {
        return NULL;
    }

    // array attributes
    int iok = -22;
    PyArrayObject* x = 0;
    double* xd = 0;
    if (job > 10) {

        size_t nsp = phase_nSpecies(ph);
#ifdef HAS_NUMPY
        npy_intp nnn = nsp;
        x = (PyArrayObject*)PyArray_SimpleNew(1,  &nnn, PyArray_DOUBLE);
#else
        int nnn = int(nsp);
        x = (PyArrayObject*)PyArray_FromDims(1, &nnn, PyArray_DOUBLE);
#endif
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
    } else {

        size_t nel = phase_nElements(ph);
#ifdef HAS_NUMPY
        npy_intp nnn = nel;
        x = (PyArrayObject*)PyArray_SimpleNew(1, &nnn, PyArray_DOUBLE);
#else
        int nnn = int(nel);
        x = (PyArrayObject*)PyArray_FromDims(1, &nnn, PyArray_DOUBLE);
#endif
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
    } else {
        PyErr_SetString(ErrorObject,"Unknown array attribute");
        return NULL;
    }
}

// string attributes
static PyObject*
phase_getstring(PyObject* self, PyObject* args)
{
    int ph, job, iok = -1;
    int k;
    int buflen;
    char* output_buf = 0;
    if (!PyArg_ParseTuple(args, "iii:phase_getstring", &ph, &job, &k)) {
        return NULL;
    }
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
        delete[] output_buf;
        return str;
    }
    delete output_buf;
    if (iok == -1) {
        return reportCanteraError();
    } else {
        PyErr_SetString(ErrorObject,"Unknown string attribute");
        return NULL;
    }
}

static PyObject*
phase_setfp(PyObject* self, PyObject* args)
{
    double vv;
    int iok = -2;
    int ph;
    int job;

    if (!PyArg_ParseTuple(args, "iid:phase_getfp", &ph, &job, &vv)) {
        return NULL;
    }

    // set floating-point attributes
    switch (job) {
    case 1:
        iok = phase_setTemperature(ph, vv);
        break;
    case 2:
        iok = phase_setDensity(ph, vv);
        break;
    case 3:
        iok = phase_setMolarDensity(ph, vv);
        break;
    default:
        iok = -10;
    }
    if (iok >= 0) {
        return Py_BuildValue("i",iok);
    } else {
        PyErr_SetString(ErrorObject,"Unknown floating-point attribute");
        return NULL;
    }
}


static PyObject*
phase_setarray(PyObject* self, PyObject* args)
{
    int ph;
    int job;
    int norm;
    int iok;
    PyObject* seq;
    if (!PyArg_ParseTuple(args, "iiiO:phase_setarray", &ph, &job, &norm, &seq)) {
        return NULL;
    }
    PyArrayObject* a = (PyArrayObject*)
                       PyArray_ContiguousFromObject(seq, PyArray_DOUBLE, 1, 1);
    double* xd = (double*)a->data;
    size_t len = a->dimensions[0];
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
    Py_DECREF(a);
    if (iok >= 0) {
        return Py_BuildValue("i",iok);
    }
    if (iok == -1) {
        return reportCanteraError();
    } else {
        PyErr_SetString(ErrorObject, "Error in phase_setarray");
        return NULL;
    }
}

static PyObject*
phase_setstring(PyObject* self, PyObject* args)
{
    int ph;
    int job;
    int iok;
    char* str;
    if (!PyArg_ParseTuple(args, "iis:phase_setstring", &ph, &job, &str)) {
        return NULL;
    }
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
    if (iok >= 0) {
        return Py_BuildValue("i",iok);
    }
    if (iok == -1) {
        return reportCanteraError();
    } else {
        PyErr_SetString(ErrorObject, "Error in phase_setstring");
        return NULL;
    }
}

