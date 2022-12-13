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
    PythonHandle(PyObject* obj, bool weak) : m_obj(obj), m_weak(weak) {}

    ~PythonHandle() {
        if (!m_weak) {
            Py_XDECREF(m_obj);
        }
    }

    void* get() {
        return m_obj;
    }

private:
    PyObject* m_obj;
    bool m_weak;
};

}

#endif
