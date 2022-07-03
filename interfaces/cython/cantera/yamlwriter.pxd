#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .units cimport *

cdef class YamlWriter:
    cdef shared_ptr[CxxYamlWriter] _writer
    cdef CxxYamlWriter* writer
    @staticmethod
    cdef CxxUnitSystem _get_unitsystem(UnitSystem units)
