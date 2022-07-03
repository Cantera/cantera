#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

cdef int assign_delegates(object, CxxDelegator*) except -1
