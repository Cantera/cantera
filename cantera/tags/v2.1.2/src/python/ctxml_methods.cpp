
static PyObject*
py_xml_new(PyObject* self, PyObject* args)
{
    char* nm;
    if (!PyArg_ParseTuple(args, "s:xml_new", &nm)) {
        return NULL;
    }
    int n = xml_new(nm);
    return Py_BuildValue("i",n);
}

static PyObject*
py_xml_get_XML_File(PyObject* self, PyObject* args)
{
    char* file;
    int debug;
    if (!PyArg_ParseTuple(args, "si:xml_get_XML_File", &file, &debug)) {
        return NULL;
    }
    int n = xml_get_XML_File(file, debug);
    if (n < 0) {
        return reportError(n);
    }
    return Py_BuildValue("i",n);
}


static PyObject*
py_xml_del(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_del", &n)) {
        return NULL;
    }
    int iok = xml_del(n);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_xml_clear(PyObject* self, PyObject* args)
{
    int iok = xml_clear();
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",0);
}

static PyObject*
py_xml_attrib(PyObject* self, PyObject* args)
{
    int n;
    char* key;
    if (!PyArg_ParseTuple(args, "is:xml_attrib", &n, &key)) {
        return NULL;
    }
    char* val = new char[81];
    int iok = xml_attrib(n, key, val);
    if (iok < 0) {
        return reportError(iok);
    }
    PyObject* r = Py_BuildValue("s",val);
    delete[] val;
    return r;
}

static PyObject*
py_xml_tag(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_tag", &n)) {
        return NULL;
    }
    char* val = new char[81];
    int iok = xml_tag(n, val);
    if (iok < 0) {
        return reportError(iok);
    }
    PyObject* r = Py_BuildValue("s",val);
    delete val;
    return r;
}

static PyObject*
py_xml_value(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_value", &n)) {
        return NULL;
    }
    char* val = new char[81];
    int iok = xml_value(n, val);
    if (iok < 0) {
        return reportError(iok);
    }
    PyObject* r = Py_BuildValue("s",val);
    delete val;
    return r;
}

static PyObject*
py_xml_child(PyObject* self, PyObject* args)
{
    int n;
    char* loc;
    if (!PyArg_ParseTuple(args, "is:xml_child", &n, &loc)) {
        return NULL;
    }
    int m = xml_child(n, loc);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_childbynumber(PyObject* self, PyObject* args)
{
    int n, k;
    if (!PyArg_ParseTuple(args, "ii:xml_childbynumber", &n, &k)) {
        return NULL;
    }
    int m = xml_child_bynumber(n, k);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_findID(PyObject* self, PyObject* args)
{
    int n;
    char* id;
    if (!PyArg_ParseTuple(args, "is:xml_findID", &n, &id)) {
        return NULL;
    }
    int m = xml_findID(n, id);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_findByName(PyObject* self, PyObject* args)
{
    int n;
    char* nm;
    if (!PyArg_ParseTuple(args, "is:xml_findID", &n, &nm)) {
        return NULL;
    }
    int m = xml_findByName(n, nm);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_nChildren(PyObject* self, PyObject* args)
{
    int n;
    if (!PyArg_ParseTuple(args, "i:xml_nChildren", &n)) {
        return NULL;
    }
    int m = xml_nChildren(n);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_addChild(PyObject* self, PyObject* args)
{
    int n;
    char* name, *value;
    if (!PyArg_ParseTuple(args, "iss:xml_addChild", &n, &name, &value)) {
        return NULL;
    }
    int m = xml_addChild(n, name, value);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}


static PyObject*
py_xml_addChildNode(PyObject* self, PyObject* args)
{
    int n, j;
    if (!PyArg_ParseTuple(args, "ii:xml_addChildNode", &n, &j)) {
        return NULL;
    }
    int m = xml_addChildNode(n, j);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}


static PyObject*
py_xml_addAttrib(PyObject* self, PyObject* args)
{
    int n;
    char* name, *value;
    if (!PyArg_ParseTuple(args, "iss:xml_addAttrib", &n, &name, &value)) {
        return NULL;
    }
    int m = xml_addAttrib(n, name, value);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_addComment(PyObject* self, PyObject* args)
{
    int n;
    char* comment;
    if (!PyArg_ParseTuple(args, "is:xml_addComment", &n, &comment)) {
        return NULL;
    }
    int m = xml_addComment(n, comment);
    if (m < 0) {
        return reportError(m);
    }
    return Py_BuildValue("i",m);
}

static PyObject*
py_xml_removeChild(PyObject* self, PyObject* args)
{
    int n, m;
    if (!PyArg_ParseTuple(args, "ii:xml_removeChild", &n, &m)) {
        return NULL;
    }
    int iok = xml_removeChild(n, m);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",iok);
}


static PyObject*
py_xml_write(PyObject* self, PyObject* args)
{
    int n;
    char* file;
    if (!PyArg_ParseTuple(args, "is:xml_write", &n, &file)) {
        return NULL;
    }
    int iok = xml_write(n, file);
    if (iok < 0) {
        return reportError(iok);
    }
    return Py_BuildValue("i",iok);
}

static PyObject*
py_ctml_getFloatArray(PyObject* self, PyObject* args)
{
    int n;
    int iconv, ia;
    if (!PyArg_ParseTuple(args, "iii", &n, &iconv, &ia)) {
        return NULL;
    }

#ifdef HAS_NUMPY
    npy_intp nia = ia;
    PyArrayObject* a =
        (PyArrayObject*)PyArray_SimpleNew(1, &nia, PyArray_DOUBLE);
#else
    PyArrayObject* a =
        (PyArrayObject*)PyArray_FromDims(1, &ia, PyArray_DOUBLE);
#endif
    double* x = (double*)a->data;
    int iok = ctml_getFloatArray(n, ia, x, iconv);
    if (iok < 0) {
        return reportError(iok);
    }
    return PyArray_Return(a);
}

