#ifndef CTPY_UTILS
#define CTPY_UTILS

#include "Python.h"

static PyObject* reportCanteraError() {
    char* buf = 0;
    int buflen = getCanteraError(0, buf);
    buf = new char[buflen+1];
    getCanteraError(buflen, buf);
    PyErr_SetString(ErrorObject,buf);
    delete buf;
    return NULL;
}

static PyObject* reportError(int n) {
    if (n == -1) return reportCanteraError();
    else if (n < 0) {
        PyErr_SetString(ErrorObject,"Exception occurred.");
    }
    return NULL;
}

// template<class V>
// int pyNumericSequence_ToVector(PyObject* seq, V& vec) 
// {
//     if (!PySequence_Check(seq)) return -1;
//     int n = PySequence_Size(seq);
//     vec.resize(n);
//     PyObject* item;
//     double x;
//     for (int i = 0; i < n; i++) {
//         item = PySequence_GetItem(seq, i);
//         if (PyFloat_Check(item))
//             x = PyFloat_AsDouble(item);
//         else if (PyInt_Check(item)) 
//             x = 1.0*PyInt_AsLong(item);
//         else
//             return -2;
//         vec[i] = x;
//     }
//     return 1;
// }


// template<class V>
// PyObject* pyNumericTuple_FromVector(const V& vec) 
// {
//     int n = vec.size();
//     PyObject* seq = PyTuple_New(n);
//     PyObject* item;
//     for (int i = 0; i < n; i++) {
//         item = PyFloat_FromDouble(1.0*vec[i]);
//         if (PyTuple_SetItem(seq, i, item) < 0) return NULL;
//     }
//     return seq;
// }


// template<class V>
// PyObject* pyNumericList_FromVector(const V& vec) 
// {
//     int n = vec.size();
//     PyObject* seq = PyList_New(n);
//     PyObject* item;
//     for (int i = 0; i < n; i++) {
//         item = PyFloat_FromDouble(1.0*vec[i]);
//         if (PyList_SetItem(seq, i, item) < 0) return NULL;
//     }
//     return seq;
// }

// template<class V>
// PyObject* pyStringTuple_FromVector(const V& vec) 
// {
//     int n = vec.size();
//     PyObject* seq = PyTuple_New(n);
//     PyObject* item;
//     for (int i = 0; i < n; i++) {
//         item = PyString_FromString(vec[i].c_str());
//         if (PyTuple_SetItem(seq, i, item) < 0) return NULL;
//     }
//     return seq;
// }


// template<class M>
// int pyStrNumMapping_ToMap(PyObject* d, M& m) 
// {
//     if (!PyMapping_Check(d)) return -1;
//     int n = PyMapping_Length(d);
//     PyObject *keys, *values, *ky, *v;
//     double x;
//     string key;
//     keys = PyMapping_Keys(d);
//     values = PyMapping_Values(d);
//     for (int i = 0; i < n; i++) {
//         ky = PySequence_GetItem(keys, i);
//         v = PySequence_GetItem(values, i);
//         key = string(PyString_AsString(ky));
//         if (PyFloat_Check(v))
//             x = PyFloat_AsDouble(v);
//         else if (PyInt_Check(v)) 
//             x = 1.0*PyInt_AsLong(v);
//         else
//             return -2;
//         m[key] = x;
//     }
//     return 1;
// }


//  template<class M>
//  PyObject* pyStrNumMapping_FromMap(M& m) 
//  {
//      if (!PyMapping_Check(d)) return -1;
//      int n = m.size();
//      double value;
//      char* key;
//      PyObject* x;
//      TYPENAME_KEYWORD M::const_iter ptr;
//      PyObject* map = PyDict_New();
//      for (int i = 0; i < n; i++) {
//          key = ptr->first.c_str();
//          value = ptr->second;
//          x = PyFloat_FromDouble(value);
//          PyMapping_SetItemString(map, key, x);
//      }
//      return map;
//  }


// static int lengthError() {
//     PyErr_SetString(PyExc_AttributeError,
//         "wrong length for sequence");
//     return -1;
// }

#endif
