#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

cdef class Func1:
    cdef shared_ptr[CxxFunc1] _func
    cdef CxxFunc1* func
    cdef object callable
    cdef object exception
    cpdef void _set_callback(self, object) except *

cdef class TabulatedFunction(Func1):
    cpdef void _set_tables(self, object, object, string) except *
