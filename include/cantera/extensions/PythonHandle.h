//! @file PythonHandle.h

#ifndef CT_PYTHONHANDLE_H
#define CT_PYTHONHANDLE_H

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Python.h"
#include "cantera/base/ExtensionManager.h"

namespace Cantera
{

//! Class that holds an owned or weak (borrowed) reference to a Python object
class PythonHandle : public ExternalHandle
{
public:
    //! Create a handle to hold a Python object
    //! @param obj  The Python object to be held
    //! @param weak  `true` if this is a weak reference to the Python object and this
    //!    handle is not responsible for deleting the Python object, or `false` if this
    //!    handle should own a reference to the Python object
    PythonHandle(PyObject* obj, bool weak) : m_obj(obj), m_weak(weak) {
        if (!weak) {
            Py_XINCREF(obj);
        }
    }
    PythonHandle(const PythonHandle&) = delete;
    PythonHandle& operator=(const PythonHandle&) = delete;

    ~PythonHandle() {
        if (!m_weak) {
            Py_XDECREF(m_obj);
        }
    }

    void* get() override {
        return m_obj;
    }

private:
    PyObject* m_obj;
    bool m_weak;
};

}

#endif
