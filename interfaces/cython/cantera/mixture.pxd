#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

cdef class Mixture:
    cdef CxxMultiPhase* mix
    cdef list _phases
    cdef object _weakref_proxy
    cpdef int element_index(self, element) except *
