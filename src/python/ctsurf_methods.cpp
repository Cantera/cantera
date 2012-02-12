
static PyObject*
py_surf_setsitedensity(PyObject* self, PyObject* args)
{
    int n;
    double s0;
    if (!PyArg_ParseTuple(args, "id:surf_setsitedensity", &n, &s0)) {
        return NULL;
    }

    int iok = surf_setsitedensity(n, s0);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_surf_sitedensity(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:surf_sitedensity", &n)) {
        return NULL;
    }

    double s0 = surf_sitedensity(n);
    return Py_BuildValue("d",s0);
}


static PyObject*
py_surf_setcoverages(PyObject* self, PyObject* args)
{
    int n;
    PyObject* cov;
    if (!PyArg_ParseTuple(args, "iO:surf_setcoverages", &n, &cov)) {
        return NULL;
    }
    double* x = (double*)((PyArrayObject*)cov)->data;
    int iok = surf_setcoverages(n, x);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_surf_setconcentrations(PyObject* self, PyObject* args)
{
    int n;
    PyObject* c;
    if (!PyArg_ParseTuple(args, "iO:surf_setconcentrations", &n, &c)) {
        return NULL;
    }
    double* x = (double*)((PyArrayObject*)c)->data;
    int iok = surf_setconcentrations(n, x);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_surf_getcoverages(PyObject* self, PyObject* args)
{
    int n;
    PyArrayObject* cov;
    if (!PyArg_ParseTuple(args, "i:surf_getcoverages", &n)) {
        return NULL;
    }
    size_t nsp = th_nSpecies(n);
#ifdef HAS_NUMPY
    npy_intp nnsp = nsp;
    cov = (PyArrayObject*)PyArray_SimpleNew(1, &nnsp, PyArray_DOUBLE);
#else
    int nnsp = int(nsp);
    cov = (PyArrayObject*)PyArray_FromDims(1, &nnsp, PyArray_DOUBLE);
#endif
    double* x = (double*)((PyArrayObject*)cov)->data;
    int iok = surf_getcoverages(n, x);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("O",cov);
}

static PyObject*
py_surf_getconcentrations(PyObject* self, PyObject* args)
{
    int n;
    PyArrayObject* c;
    if (!PyArg_ParseTuple(args, "i:surf_getconcentrations", &n)) {
        return NULL;
    }
    size_t nsp = th_nSpecies(n);
#ifdef HAS_NUMPY
    npy_intp nnsp = nsp;
    c = (PyArrayObject*)PyArray_SimpleNew(1, &nnsp, PyArray_DOUBLE);
#else
    int nnsp = int(nsp);
    c = (PyArrayObject*)PyArray_FromDims(1, &nnsp, PyArray_DOUBLE);
#endif
    double* x = (double*)((PyArrayObject*)c)->data;
    int iok = surf_getconcentrations(n, x);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("O",c);
}



