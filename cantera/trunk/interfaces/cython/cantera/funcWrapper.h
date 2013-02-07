#ifndef CT_CYTHON_FUNC_WRAPPER
#define CT_CYTHON_FUNC_WRAPPER

#include "cantera/numerics/Func1.h"

typedef double(*callback_wrapper)(double, void*, void**);

// A C++ exception that holds a Python exception so that it can be re-raised
// by translate_exception()
class CallbackError : public Cantera::CanteraError
{
public:
    CallbackError(void* type, void* value) :
        m_type((PyObject*) type),
        m_value((PyObject*) value)
    {}
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
    } catch (CallbackError& exn) {
        // Re-raise a Python exception generated in a callback
        PyErr_SetObject(exn.m_type, exn.m_value);
    } catch (const std::out_of_range& exn) {
        PyErr_SetString(PyExc_IndexError, exn.what());
    } catch (const std::exception& exn) {
        PyErr_SetString(PyExc_Exception, exn.what());
    } catch (...) {
        PyErr_SetString(PyExc_Exception, "Unknown exception");
    }
    return 0;
}

#endif
