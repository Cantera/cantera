
// static PyObject *
// ct_buildSolutionFromXML(PyObject *self, PyObject *args)
// {
//     int ixml, ith, ikin;
//     char *src=0, *id=0;
//     if (!PyArg_ParseTuple(args, "sisii:buildSolutionFromXML", &src, &ixml,
//             &id, &ith, &ikin))
//         return NULL;
//     int ok = buildSolutionFromXML(src, ixml, id, ith, ikin);
//     if (ok == -1) { return reportCanteraError();}
//     return Py_BuildValue("i",ok);
// }

static PyObject*
ct_get_cantera_error(PyObject* self, PyObject* args)
{
    char* buf = new char[400];
    getCanteraError(400, buf);
    PyObject* msg = Py_BuildValue("s",buf);
    delete buf;
    return msg;
}

static PyObject*
ct_refcnt(PyObject* self, PyObject* args)
{
    PyObject* o;
    if (!PyArg_ParseTuple(args, "O", &o)) {
        return NULL;
    }
    return Py_BuildValue("i",o->ob_refcnt);
}

// static PyObject *
// ct_print(PyObject *self, PyObject *args)
// {
//     char* msg;
//     if (!PyArg_ParseTuple(args, "s:print", &msg))
//         return NULL;
//     printf(msg);
//     return Py_BuildValue("i",0);
// }

static PyObject*
ct_addDirectory(PyObject* self, PyObject* args)
{
    char* dir;
    if (!PyArg_ParseTuple(args, "s:addDirectory", &dir)) {
        return NULL;
    }
    size_t n = strlen(dir);
    addCanteraDirectory(n, dir);
    return Py_BuildValue("i",0);
}

// static PyObject *
// ct_readlog(PyObject *self, PyObject *args)
// {
//     char* msg = 0;
//     int n = readlog(-1, msg);
//     if (n > 0) {
//         msg = new char[n+1];
//         readlog(n, msg);
//         PyObject* r = Py_BuildValue("s",msg);
//         return r;
//     }
//     else
//        return Py_BuildValue("s","");
//}

static PyObject*
ct_ck2cti(PyObject* self, PyObject* args)
{
    int iok;
    char* infile, *thermo, *tran, *idtag;
    int debug, validate;
    if (!PyArg_ParseTuple(args, "ssssii:ck2cti", &infile,
                          &thermo, &tran, &idtag, &debug, &validate)) {
        return NULL;
    }
    iok = ck_to_cti(infile, thermo, tran, idtag, debug, validate);
    if (iok == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",iok);
}

static PyObject*
ct_writelogfile(PyObject* self, PyObject* args)
{
    int iok;
    char* logfile;
    if (!PyArg_ParseTuple(args, "s:writelogfile", &logfile)) {
        return NULL;
    }
    iok = writelogfile(logfile);
    if (iok == -1) {
        return reportCanteraError();
    }
    return Py_BuildValue("i",iok);
}

static PyObject* ct_get_version(PyObject* self, PyObject* args)
{
    return Py_BuildValue("s",std::string(CANTERA_VERSION).c_str());
}
