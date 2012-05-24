/**
 * @file ctkinetics_methods.cpp
 *
 */
static PyObject*
kin_newFromXML(PyObject* self, PyObject* args)
{
    int mxml, iphase, neighbor1, neighbor2, neighbor3, neighbor4;
    if (!PyArg_ParseTuple(args, "iiiiii:newFromXML", &mxml,
                          &iphase, &neighbor1, &neighbor2, &neighbor3, &neighbor4)) {
        return NULL;
    }
    int n = int(newKineticsFromXML(mxml, iphase, neighbor1, neighbor2,
                                   neighbor3, neighbor4));
    if (n < 0) {
        return reportError(n);
    }
    return Py_BuildValue("i",n);
}

static PyObject*
kin_delete(PyObject* self, PyObject* args)
{
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_delete", &kin)) {
        return NULL;
    }
    delKinetics(kin);
    return Py_BuildValue("i",0);
}

static PyObject*
kin_phase(PyObject* self, PyObject* args)
{
    int kin, n;
    if (!PyArg_ParseTuple(args, "ii:kin_phase", &kin, &n)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_phase(kin, n));
}

static PyObject*
kin_nspecies(PyObject* self, PyObject* args)
{
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_nspecies", &kin)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_nSpecies(kin));
}

static PyObject*
kin_rstoichcoeff(PyObject* self, PyObject* args)
{
    int kin, i, k;
    if (!PyArg_ParseTuple(args, "iii:kin_rstoichcoeff", &kin, &k, &i)) {
        return NULL;
    }
    return Py_BuildValue("d",kin_reactantStoichCoeff(kin, k, i));
}

static PyObject*
kin_pstoichcoeff(PyObject* self, PyObject* args)
{
    int kin, i, k;
    if (!PyArg_ParseTuple(args, "iii:kin_pstoichcoeff", &kin, &k, &i)) {
        return NULL;
    }
    return Py_BuildValue("d",kin_productStoichCoeff(kin, k, i));
}

static PyObject*
kin_nrxns(PyObject* self, PyObject* args)
{
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_nreactions", &kin)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_nReactions(kin));
}

static PyObject*
kin_nPhases(PyObject* self, PyObject* args)
{
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_nPhases", &kin)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_nPhases(kin));
}

static PyObject*
kin_phaseIndex(PyObject* self, PyObject* args)
{
    int kin;
    char* ph;
    if (!PyArg_ParseTuple(args, "is:kin_phaseIndex", &kin, &ph)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_phaseIndex(kin, ph));
}

static PyObject*
kin_reactionPhaseIndex(PyObject* self, PyObject* args)
{
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_reactionPhaseIndex", &kin)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_reactionPhaseIndex(kin));
}

static PyObject*
kin_isrev(PyObject* self, PyObject* args)
{
    int kin, i;
    if (!PyArg_ParseTuple(args, "ii:kin_isrev", &kin, &i)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_isReversible(kin,i));
}

static PyObject*
kin_rxntype(PyObject* self, PyObject* args)
{
    int kin, i;
    if (!PyArg_ParseTuple(args, "ii:kin_rxntype", &kin, &i)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_reactionType(kin,i));
}

static PyObject*
kin_multiplier(PyObject* self, PyObject* args)
{
    int kin, i;
    if (!PyArg_ParseTuple(args, "ii:kin_multiplier", &kin, &i)) {
        return NULL;
    }
    return Py_BuildValue("d",kin_multiplier(kin,i));
}

static PyObject*
kin_setMultiplier(PyObject* self, PyObject* args)
{
    int kin, i;
    double v;
    if (!PyArg_ParseTuple(args, "iid:kin_setMultiplier", &kin, &i, &v)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_setMultiplier(kin,i,v));
}


static PyObject*
kin_type(PyObject* self, PyObject* args)
{
    int kin;
    if (!PyArg_ParseTuple(args, "i:kin_type", &kin)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_type(kin));
}

static PyObject*
kin_start(PyObject* self, PyObject* args)
{
    int kin, p;
    if (!PyArg_ParseTuple(args, "ii:kin_start", &kin, &p)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_start(kin,p));
}

static PyObject*
kin_speciesIndex(PyObject* self, PyObject* args)
{
    int kin;
    char* nm, *ph;
    if (!PyArg_ParseTuple(args, "iss:kin_speciesIndex", &kin, &nm, &ph)) {
        return NULL;
    }
    return Py_BuildValue("i",kin_speciesIndex(kin,nm,ph));
}

static PyObject*
kin_advanceCoverages(PyObject* self, PyObject* args)
{
    int kin;
    double dt;
    if (!PyArg_ParseTuple(args, "id:kin_advanceCoverages", &kin, &dt)) {
        return NULL;
    }
    int iok = kin_advanceCoverages(kin, dt);
    if (iok < 0) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",0);
}

static PyObject*
kin_getarray(PyObject* self, PyObject* args)
{
    int kin;
    int job;

    if (!PyArg_ParseTuple(args, "ii:kin_getarray", &kin, &job)) {
        return NULL;
    }

    // array attributes
    int iok = -22;
    size_t nrxns = kin_nReactions(kin);
    size_t nsp = kin_nSpecies(kin);
    size_t ix;
    if (job < 45 || job >= 90) {
        ix = nrxns;
    } else {
        ix = nsp;
    }

#ifdef HAS_NUMPY
    npy_intp nix = ix;
    PyArrayObject* x = (PyArrayObject*)PyArray_SimpleNew(1, &nix, PyArray_DOUBLE);
#else
    int nix = int(ix);
    PyArrayObject* x =
        (PyArrayObject*)PyArray_FromDims(1, &nix, PyArray_DOUBLE);
#endif
    double* xd = (double*)x->data;

    switch (job) {
    case 10:
        iok = kin_getFwdRatesOfProgress(kin, nrxns, xd);
        break;
    case 20:
        iok = kin_getRevRatesOfProgress(kin, nrxns, xd);
        break;
    case 30:
        iok = kin_getNetRatesOfProgress(kin, nrxns, xd);
        break;
    case 32:
        iok = kin_getActivationEnergies(kin, nrxns, xd);
        break;
    case 34:
        iok = kin_getFwdRateConstants(kin, nrxns, xd);
        break;
    case 35:
        iok = kin_getRevRateConstants(kin, 1, nrxns, xd);
        break;
    case 36:
        iok = kin_getRevRateConstants(kin, 0, nrxns, xd);
        break;
    case 40:
        iok = kin_getEquilibriumConstants(kin, nrxns, xd);
        break;
    case 50:
        iok = kin_getCreationRates(kin, nsp, xd);
        break;
    case 60:
        iok = kin_getDestructionRates(kin, nsp, xd);
        break;
    case 70:
        iok = kin_getNetProductionRates(kin, nsp, xd);
        break;
    case 80:
        iok = kin_getSourceTerms(kin, nsp, xd);
        break;
    case 90:
        iok = kin_getDelta(kin, 0, nrxns, xd);
        break;
    case 91:
        iok = kin_getDelta(kin, 1, nrxns, xd);
        break;
    case 92:
        iok = kin_getDelta(kin, 2, nrxns, xd);
        break;
    case 93:
        iok = kin_getDelta(kin, 3, nrxns, xd);
        break;
    case 94:
        iok = kin_getDelta(kin, 4, nrxns, xd);
        break;
    case 95:
        iok = kin_getDelta(kin, 5, nrxns, xd);
        break;
    default:
        ;
    }
    if (iok >= 0) {
        return PyArray_Return(x);
    } else {
        return reportError(iok);
    }
}



// string attributes
static PyObject*
kin_getstring(PyObject* self, PyObject* args)
{
    int kin, job, i, iok = -3;
    int buflen;
    char* output_buf = 0;
    if (!PyArg_ParseTuple(args, "iii:kin_getstring", &kin, &job, &i)) {
        return NULL;
    }
    switch (job) {
    case 1:
        buflen = 80;
        output_buf = new char[buflen];
        iok = kin_getReactionString(kin, i, buflen, output_buf);
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
    if (iok == -1) {
        return reportCanteraError();
    } else {
        PyErr_SetString(ErrorObject,"Unknown string attribute");
        return NULL;
    }
}
