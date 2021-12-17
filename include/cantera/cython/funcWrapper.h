// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CYTHON_FUNC_WRAPPER
#define CT_CYTHON_FUNC_WRAPPER

#include "cantera/numerics/Func1.h"
#include "cantera/base/ctexceptions.h"
#include <stdexcept>

#define CANTERA_USE_INTERNAL
#include "cantera/clib/clib_defs.h"
#include "Python.h"

typedef double(*callback_wrapper)(double, void*, void**);

//! A class to hold information needed to call Python functions from delgated
//! methods (see class Delegator).
class PyFuncInfo {
public:
    PyFuncInfo()
        : m_func(nullptr)
        , m_exception_type(nullptr)
        , m_exception_value(nullptr)
    {
    }

    PyFuncInfo(const PyFuncInfo& other)
        : m_func(other.m_func)
        , m_exception_type(other.m_exception_type)
        , m_exception_value(other.m_exception_value)
    {
        Py_XINCREF(m_exception_type);
        Py_XINCREF(m_exception_value);
    }

    ~PyFuncInfo() {
        Py_XDECREF(m_exception_type);
        Py_XDECREF(m_exception_value);
    }

    PyObject* func() {
        return m_func;
    }
    void setFunc(PyObject* f) {
        m_func = f;
    }

    PyObject* exceptionType() {
        return m_exception_type;
    }
    void setExceptionType(PyObject* obj) {
        Py_XDECREF(m_exception_type);
        Py_XINCREF(obj);
        m_exception_type = obj;
    }

    PyObject* exceptionValue() {
        return m_exception_value;
    }
    void setExceptionValue(PyObject* obj) {
        Py_XDECREF(m_exception_value);
        Py_XINCREF(obj);
        m_exception_value = obj;
    }

private:
    PyObject* m_func;
    PyObject* m_exception_type;
    PyObject* m_exception_value;
};


// A C++ exception that holds a Python exception so that it can be re-raised
// by translate_exception()
class CallbackError : public Cantera::CanteraError
{
public:
    //! Constructor used by Func1Py
    CallbackError(void* type, void* value) :
        CanteraError("Python callback function"),
        m_type((PyObject*) type),
        m_value((PyObject*) value)
    {
        Py_XINCREF(m_type);
        Py_XINCREF(m_value);
    }

    //! Constructor used by pyOverride()
    explicit CallbackError(PyFuncInfo& info) :
        CanteraError("Python callback function"),
        m_type(info.exceptionType()),
        m_value(info.exceptionValue())
    {
        Py_XINCREF(m_type);
        Py_XINCREF(m_value);
        info.setExceptionType(0);
        info.setExceptionValue(0);
    }

    ~CallbackError() {
        Py_XDECREF(m_type);
        Py_XDECREF(m_value);
    }

    std::string getMessage() const {
        std::string msg;

        PyObject* name = PyObject_GetAttrString(m_type, "__name__");
        PyObject* value_str = PyObject_Str(m_value);

        PyObject* name_bytes = PyUnicode_AsASCIIString(name);
        PyObject* value_bytes = PyUnicode_AsASCIIString(value_str);

        if (name_bytes) {
            msg += PyBytes_AsString(name_bytes);
            Py_DECREF(name_bytes);
        } else {
            msg += "<error determining exception type>";
        }

        msg += ": ";

        if (value_bytes) {
            msg += PyBytes_AsString(value_bytes);
            Py_DECREF(value_bytes);
        } else {
            msg += "<error determining exception message>";
        }

        Py_XDECREF(name);
        Py_XDECREF(value_str);
        return msg;
    }

    virtual std::string getClass() const {
        return "Exception";
    }

    PyObject* m_type;
    PyObject* m_value;
};


// A function of one variable implemented as a callable Python object
class Func1Py : public Cantera::Func1
{
public:
    Func1Py(callback_wrapper callback, void* pyobj) :
        m_callback(callback),
        m_pyobj(pyobj) {
    }

    double eval(double t) const {
        void* err[2] = {0, 0};
        double y = m_callback(t, m_pyobj, err);
        if (err[0]) {
            throw CallbackError(err[0], err[1]);
        }
        return y;
    }

private:
    callback_wrapper m_callback;
    void* m_pyobj;
};

extern "C" {
    CANTERA_CAPI PyObject* pyCanteraError;
}



//! Take a function which requires Python function information (as a PyFuncInfo
//! object) and capture that object to generate a function that does not reqire
//! any Python-specific arguments.
//!
//! The inner function is responsible for catching Python exceptions and
//! stashing their details in the PyFuncInfo object. The wrapper function
//! generated here examines the stashed exception and throws a C++ exception
//!
//! The caller of pyOverride must continue holding a reference to pyFunc for the
//! lifetime of the pyOverride object to prevent it from being prematurely garbage
//! collected.
template <class ... Args>
std::function<void(Args ...)> pyOverride(PyObject* pyFunc, void func(PyFuncInfo&, Args ... args)) {
    PyFuncInfo func_info;
    func_info.setFunc(pyFunc);
    return [func_info, func](Args ... args) mutable {
        func(func_info, args ...);
        if (func_info.exceptionType()) {
            throw CallbackError(func_info);
        }
    };
}

//! Same as above, but for functions that return an int.
template <class ... Args>
std::function<int(Args ...)> pyOverride(PyObject* pyFunc, int func(PyFuncInfo&, Args ... args)) {
    PyFuncInfo func_info;
    func_info.setFunc(pyFunc);
    return [func_info, func](Args ... args) mutable {
        int ret = func(func_info, args ...);
        if (func_info.exceptionType()) {
            throw CallbackError(func_info);
        }
        return ret;
    };
}

// Translate C++ Exceptions generated by Cantera to appropriate Python
// exceptions. Used with Cython function declarations, e.g:
//     cdef double eval(double) except +translate_exception
inline int translate_exception()
{
    try {
        if (!PyErr_Occurred()) {
            // Let the latest Python exception pass through and ignore the
            // current one.
            throw;
        }
    } catch (const CallbackError& exn) {
        // Re-raise a Python exception generated in a callback
        PyErr_SetObject(exn.m_type, exn.m_value);
    } catch (const std::out_of_range& exn) {
        PyErr_SetString(PyExc_IndexError, exn.what());
    } catch (const Cantera::CanteraError& exn) {
        PyErr_SetString(pyCanteraError, exn.what());
    } catch (const std::exception& exn) {
        PyErr_SetString(PyExc_RuntimeError, exn.what());
    } catch (...) {
        PyErr_SetString(PyExc_Exception, "Unknown exception");
    }
    return 0;
}

#endif
