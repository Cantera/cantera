#ifndef CT_PYUTILS_H
#define CT_PYUTILS_H

#include "Python.h"

static PyObject* reportCanteraError()
{
    char* buf = 0;
    int buflen = getCanteraError(0, buf) + 1;
    buf = new char[buflen+1];
    getCanteraError(buflen, buf);
    PyErr_SetString(ErrorObject,buf);
    delete buf;
    //writelogfile("log.html");
    return NULL;
}

static PyObject* reportError(int n)
{
    if (n == -1) {
        return reportCanteraError();
    } else if (n < 0) {
        PyErr_SetString(ErrorObject,"Exception occurred.");
    }
    return NULL;
}


#endif
