# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *
from .func1 cimport *
from .units cimport CxxUnitStack

cdef extern from "<array>" namespace "std" nogil:
    cdef cppclass size_array1 "std::array<size_t, 1>":
        size_array1() except+
        size_t& operator[](size_t)

    cdef cppclass size_array2 "std::array<size_t, 2>":
        size_array2() except+
        size_t& operator[](size_t)

    cdef cppclass size_array3 "std::array<size_t, 3>":
        size_array3() except+
        size_t& operator[](size_t)


cdef extern from "cantera/extensions/PythonHandle.h" namespace "Cantera":
    cdef cppclass CxxExternalHandle "Cantera::ExternalHandle":
        pass

    cdef cppclass CxxPythonHandle "Cantera::PythonHandle" (CxxExternalHandle):
        CxxPythonHandle(PyObject*, cbool)
        PyObject* get()


cdef extern from "cantera/base/Delegator.h" namespace "Cantera":
    cdef cppclass CxxDelegator "Cantera::Delegator":
        Delegator()

        void setDelegatorName(string&)
        void holdExternalHandle(string&, shared_ptr[CxxExternalHandle]&)
        shared_ptr[CxxExternalHandle] getExternalHandle(string&)

        void setDelegate(string&, function[void()], string&) except +translate_exception
        void setDelegate(string&, function[void(cbool)], string&) except +translate_exception
        void setDelegate(string&, function[void(double)], string&) except +translate_exception
        void setDelegate(string&, function[void(CxxAnyMap&)], string&) except +translate_exception
        void setDelegate(string&, function[void(const CxxAnyMap&, const CxxUnitStack&)], string&) except +translate_exception
        void setDelegate(string&, function[void(const string&, void*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array1, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array1, double, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array2, double, double*, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array3, double*, double*, double*)], string&) except +translate_exception
        void setDelegate(string&, function[int(double&, void*)], string&) except +translate_exception
        void setDelegate(string&, function[int(string&, size_t)], string&) except +translate_exception
        void setDelegate(string&, function[int(size_t&, string&)], string&) except +translate_exception


cdef extern from "cantera/cython/funcWrapper.h":
    # pyOverride is actually a templated function, but we have to specify the individual
    # instantiations because Cython doesn't understand variadic templates
    cdef function[void(double)] pyOverride(PyObject*, void(PyFuncInfo&, double))
    cdef function[void(cbool)] pyOverride(PyObject*, void(PyFuncInfo&, cbool))
    cdef function[void()] pyOverride(PyObject*, void(PyFuncInfo&))
    cdef function[void(CxxAnyMap&)] pyOverride(PyObject*, void(PyFuncInfo&, CxxAnyMap&))
    cdef function[void(const CxxAnyMap&, const CxxUnitStack&)] pyOverride(
        PyObject*, void(PyFuncInfo&, const CxxAnyMap&, const CxxUnitStack&))
    cdef function[void(const string&, void*)] pyOverride(
        PyObject*, void(PyFuncInfo&, const string&, void*))
    cdef function[void(size_array1, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array1, double*))
    cdef function[void(size_array1, double, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array1, double, double*))
    cdef function[void(size_array2, double, double*, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array2, double, double*, double*))
    cdef function[void(size_array3, double*, double*, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array3, double*, double*, double*))
    cdef function[int(double&, void*)] pyOverride(PyObject*, int(PyFuncInfo&, double&, void*))
    cdef function[int(string&, size_t)] pyOverride(PyObject*, int(PyFuncInfo&, string&, size_t))
    cdef function[int(size_t&, const string&)] pyOverride(
        PyObject*, int(PyFuncInfo&, size_t&, const string&))

cdef extern from "cantera/base/ExtensionManager.h" namespace "Cantera":
    cdef cppclass CxxExtensionManager "Cantera::ExtensionManager":
        void registerRateBuilder(string&, string&, string&) except +translate_exception
        void registerRateDataBuilder(string&, string&, string&) except +translate_exception

        shared_ptr[CxxExtensionManager] build(string&)


cdef extern from "cantera/base/ExtensionManagerFactory.h" namespace "Cantera":
    cdef cppclass CxxExtensionManagerFactory "Cantera::ExtensionManagerFactory":
        @staticmethod
        shared_ptr[CxxExtensionManager] build(string&)


cdef extern from "cantera/extensions/PythonExtensionManager.h" namespace "Cantera":
    cdef cppclass CxxPythonExtensionManager "Cantera::PythonExtensionManager" (CxxExtensionManager):
        @staticmethod
        void registerSelf()


cdef extern from "cantera/base/ExtensionManagerFactory.h" namespace "Cantera":
    cdef cppclass CxxExtensionManagerFactory "Cantera::ExtensionManagerFactory":
        @staticmethod
        shared_ptr[CxxExtensionManager] build(string&)


ctypedef CxxDelegator* CxxDelegatorPtr

cdef int assign_delegates(object, CxxDelegator*) except -1
cdef void callback_v(PyFuncInfo& funcInfo) noexcept
