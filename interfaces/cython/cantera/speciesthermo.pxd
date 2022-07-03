#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

cdef class SpeciesThermo:
    cdef shared_ptr[CxxSpeciesThermo] _spthermo
    cdef CxxSpeciesThermo* spthermo
    cdef _assign(self, shared_ptr[CxxSpeciesThermo] other)

cdef wrapSpeciesThermo(shared_ptr[CxxSpeciesThermo] spthermo)
