# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *
from .func1 cimport *

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

cdef extern from "cantera/base/Delegator.h" namespace "Cantera":
    cdef cppclass CxxDelegator "Cantera::Delegator":
        Delegator()

        void setDelegate(string&, function[void()], string&) except +translate_exception
        void setDelegate(string&, function[void(cbool)], string&) except +translate_exception
        void setDelegate(string&, function[void(double)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array1, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array1, double, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array2, double, double*, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array3, double*, double*, double*)], string&) except +translate_exception
        void setDelegate(string&, function[int(string&, size_t)], string&) except +translate_exception
        void setDelegate(string&, function[int(size_t&, string&)], string&) except +translate_exception

cdef extern from "cantera/cython/funcWrapper.h":
    # pyOverride is actually a templated function, but we have to specify the individual
    # instantiations because Cython doesn't understand variadic templates
    cdef function[void(double)] pyOverride(PyObject*, void(PyFuncInfo&, double))
    cdef function[void(cbool)] pyOverride(PyObject*, void(PyFuncInfo&, cbool))
    cdef function[void()] pyOverride(PyObject*, void(PyFuncInfo&))
    cdef function[void(size_array1, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array1, double*))
    cdef function[void(size_array1, double, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array1, double, double*))
    cdef function[void(size_array2, double, double*, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array2, double, double*, double*))
    cdef function[void(size_array3, double*, double*, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array3, double*, double*, double*))
    cdef function[int(string&, size_t)] pyOverride(PyObject*, int(PyFuncInfo&, string&, size_t))
    cdef function[int(size_t&, const string&)] pyOverride(
        PyObject*, int(PyFuncInfo&, size_t&, const string&))


ctypedef CxxDelegator* CxxDelegatorPtr

cdef int assign_delegates(object, CxxDelegator*) except -1
