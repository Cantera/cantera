# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *


cdef extern from "cantera/numerics/Func1.h":
    cdef cppclass CxxFunc1 "Cantera::Func1":
        double eval(double) except +translate_exception
        string type()
        string typeName()
        string write(string)


cdef extern from "cantera/cython/funcWrapper.h":
    ctypedef double (*callback_wrapper)(double, void*, void**) except? 0.0
    cdef int translate_exception()

    cdef cppclass CxxFunc1Py "Func1Py" (CxxFunc1):
        CxxFunc1Py(callback_wrapper, void*)

    cdef cppclass PyFuncInfo:
        PyFuncInfo()
        PyObject* func()
        void setFunc(PyObject*)
        PyObject* exceptionType()
        void setExceptionType(PyObject*)
        PyObject* exceptionValue()
        void setExceptionValue(PyObject*)


cdef extern from "cantera/numerics/Func1Factory.h":
    cdef shared_ptr[CxxFunc1] CxxNewFunc1 "Cantera::newFunc1" (
        string, double) except +translate_exception
    cdef shared_ptr[CxxFunc1] CxxNewFunc1 "Cantera::newFunc1" (
        string, vector[double]&) except +translate_exception
    cdef shared_ptr[CxxFunc1] CxxNewFunc1 "Cantera::newFunc1" (
        string, shared_ptr[CxxFunc1], shared_ptr[CxxFunc1]) except +translate_exception
    cdef shared_ptr[CxxFunc1] CxxNewFunc1 "Cantera::newFunc1" (
        string, shared_ptr[CxxFunc1], double) except +translate_exception


cdef class Func1:
    cdef shared_ptr[CxxFunc1] _func
    cdef CxxFunc1* func
    cdef object callable
    cdef object exception
    cpdef void _set_callback(self, object) except *
