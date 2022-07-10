# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

cdef extern from "cantera/cython/funcWrapper.h":
    ctypedef double (*callback_wrapper)(double, void*, void**)
    cdef int translate_exception()

    cdef cppclass CxxFunc1 "Func1Py":
        CxxFunc1(callback_wrapper, void*)
        double eval(double) except +translate_exception

    cdef cppclass PyFuncInfo:
        PyFuncInfo()
        PyObject* func()
        void setFunc(PyObject*)
        PyObject* exceptionType()
        void setExceptionType(PyObject*)
        PyObject* exceptionValue()
        void setExceptionValue(PyObject*)


cdef extern from "cantera/numerics/Func1.h":
    cdef cppclass CxxTabulated1 "Cantera::Tabulated1":
        CxxTabulated1(int, double*, double*, string) except +translate_exception
        double eval(double) except +translate_exception


cdef class Func1:
    cdef shared_ptr[CxxFunc1] _func
    cdef CxxFunc1* func
    cdef object callable
    cdef object exception
    cpdef void _set_callback(self, object) except *

cdef class TabulatedFunction(Func1):
    cpdef void _set_tables(self, object, object, string) except *
