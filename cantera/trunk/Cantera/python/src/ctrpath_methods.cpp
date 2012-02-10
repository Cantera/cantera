
static PyObject*
py_rdiag_new(PyObject* self, PyObject* args)
{
    int iok = rdiag_new();
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",iok);
}

static PyObject*
py_rdiag_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rdiag_del", &n)) {
        return NULL;
    }
    int iok = rdiag_del(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_detailed(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rdiag_detailed", &n)) {
        return NULL;
    }
    int iok = rdiag_detailed(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_brief(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rdiag_brief", &n)) {
        return NULL;
    }
    int iok = rdiag_brief(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setThreshold(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setThreshold", &n, &v)) {
        return NULL;
    }
    int iok = rdiag_setThreshold(n,v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setFont(PyObject* self, PyObject* args)
{
    int n;
    char* font;
    if (!PyArg_ParseTuple(args, "is:rdiag_setFont", &n, &font)) {
        return NULL;
    }
    int iok = rdiag_setFont(n, font);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setBoldColor(PyObject* self, PyObject* args)
{
    int n;
    char* color;
    if (!PyArg_ParseTuple(args, "is:rdiag_setBoldColor", &n, &color)) {
        return NULL;
    }
    int iok = rdiag_setBoldColor(n, color);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setNormalColor(PyObject* self, PyObject* args)
{
    int n;
    char* color;
    if (!PyArg_ParseTuple(args, "is:rdiag_setNormalColor", &n, &color)) {
        return NULL;
    }
    int iok = rdiag_setNormalColor(n, color);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setDashedColor(PyObject* self, PyObject* args)
{
    int n;
    char* color;
    if (!PyArg_ParseTuple(args, "is:rdiag_setDashedColor", &n, &color)) {
        return NULL;
    }
    int iok = rdiag_setDashedColor(n,color);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setDotOptions(PyObject* self, PyObject* args)
{
    int n;
    char* opt;
    if (!PyArg_ParseTuple(args, "is:rdiag_setDotOptions", &n, &opt)) {
        return NULL;
    }
    int iok = rdiag_setDotOptions(n,opt);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setBoldThreshold(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setBoldThreshold", &n, &v)) {
        return NULL;
    }
    int iok = rdiag_setBoldThreshold(n,v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setNormalThreshold(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setNormalThreshold", &n, &v)) {
        return NULL;
    }
    int iok = rdiag_setNormalThreshold(n,v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setLabelThreshold(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setLabelThreshold", &n, &v)) {
        return NULL;
    }
    int iok = rdiag_setLabelThreshold(n,v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setScale(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setScale", &n, &v)) {
        return NULL;
    }
    int iok = rdiag_setScale(n,v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setFlowType(PyObject* self, PyObject* args)
{
    int n;
    int iflow;
    if (!PyArg_ParseTuple(args, "ii:rdiag_setFlowType", &n, &iflow)) {
        return NULL;
    }
    int iok = rdiag_setFlowType(n, iflow);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setArrowWidth(PyObject* self, PyObject* args)
{
    int n;
    double v;
    if (!PyArg_ParseTuple(args, "id:rdiag_setArrowWidth", &n, &v)) {
        return NULL;
    }
    int iok = rdiag_setArrowWidth(n,v);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}


static PyObject*
py_rdiag_displayOnly(PyObject* self, PyObject* args)
{
    int n, k;
    if (!PyArg_ParseTuple(args, "ii:rdiag_displayOnly", &n, &k)) {
        return NULL;
    }
    int iok = rdiag_displayOnly(n,k);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_setTitle(PyObject* self, PyObject* args)
{
    int n;
    char* t;
    if (!PyArg_ParseTuple(args, "is:rdiag_setTitle", &n, &t)) {
        return NULL;
    }
    int iok = rdiag_setTitle(n,t);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_add(PyObject* self, PyObject* args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:rdiag_add", &n, &m)) {
        return NULL;
    }
    int iok = rdiag_add(n,m);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_findMajor(PyObject* self, PyObject* args)
{
    int n;
    double thresh;
    PyObject* a;
    if (!PyArg_ParseTuple(args, "idO:rdiag_findMajor", &n, &thresh, &a)) {
        return NULL;
    }
    PyArrayObject* aa = (PyArrayObject*)a;
    size_t lda = aa->dimensions[0];
    double* x = (double*)aa->data;
    int iok = rdiag_findMajor(n, thresh, lda, x);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rdiag_write(PyObject* self, PyObject* args)
{
    int n, fmt;
    char* nm;
    if (!PyArg_ParseTuple(args, "iis:rdiag_write", &n, &fmt, &nm)) {
        return NULL;
    }
    int iok = rdiag_write(n, fmt, nm);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rbuild_new(PyObject* self, PyObject* args)
{
    int iok = rbuild_new();
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",iok);
}

static PyObject*
py_rbuild_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:rbuild_del", &n)) {
        return NULL;
    }
    int iok = rbuild_del(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rbuild_init(PyObject* self, PyObject* args)
{
    int n;
    char* log;
    int k;
    if (!PyArg_ParseTuple(args, "isi:rbuild_init", &n, &log, &k)) {
        return NULL;
    }
    int iok = rbuild_init(n,log,k);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_rbuild_build(PyObject* self, PyObject* args)
{
    int n;
    int k, idiag, iquiet;
    char* el, *dotfile;
    if (!PyArg_ParseTuple(args, "iissii:rbuild_build", &n, &k,
                          &el, &dotfile, &idiag, &iquiet)) {
        return NULL;
    }
    int iok = rbuild_build(n,k,el,dotfile,idiag,iquiet);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}


