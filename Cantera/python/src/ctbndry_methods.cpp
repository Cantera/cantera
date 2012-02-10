
static PyObject*
py_bndry_new(PyObject* self, PyObject* args)
{
    int itype;
    if (!PyArg_ParseTuple(args, "i:bndry_new", &itype)) {
        return NULL;
    }
    int iok = bndry_new(itype);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",iok);
}

static PyObject*
py_bndry_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bndry_del", &n)) {
        return NULL;
    }
    bndry_del(n);
    return Py_BuildValue("i",0);
}

static PyObject*
py_bndry_temperature(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bndry_temperature", &n)) {
        return NULL;
    }
    double t = bndry_temperature(n);
    return Py_BuildValue("d",t);
}

static PyObject*
py_bndry_settemperature(PyObject* self, PyObject* args)
{
    int n;
    double t;
    if (!PyArg_ParseTuple(args, "id:bndry_settemperature", &n, &t)) {
        return NULL;
    }
    int iok = bndry_settemperature(n, t);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}


static PyObject*
py_bndry_setspreadrate(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:bndry_setspreadrate", &n, &v)) {
        return NULL;
    }
    int iok = bndry_setSpreadRate(n, v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_bndry_mdot(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:bndry_mdot", &n)) {
        return NULL;
    }
    double mdot = bndry_mdot(n);
    return Py_BuildValue("d",mdot);
}

static PyObject*
py_bndry_setmdot(PyObject* self, PyObject* args)
{
    int n;
    double mdot;
    if (!PyArg_ParseTuple(args, "id:bndry_setmdot", &n, &mdot)) {
        return NULL;
    }
    int iok = bndry_setmdot(n, mdot);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_bndry_setxinbyname(PyObject* self, PyObject* args)
{
    int n;
    char* xin;
    //PyObject* o;
    if (!PyArg_ParseTuple(args, "is:bndry_setxin", &n, &xin)) {
        return NULL;
    }
    //double* x = (double*)((PyArrayObject*)o)->data;
    int iok = bndry_setxinbyname(n, xin);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}


static PyObject*
py_bndry_setxin(PyObject* self, PyObject* args)
{
    int n;
    PyObject* xin;
    if (!PyArg_ParseTuple(args, "iO:bndry_setxin", &n, &xin)) {
        return NULL;
    }
    double* x = (double*)((PyArrayObject*)xin)->data;
    int iok = bndry_setxin(n, x);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

