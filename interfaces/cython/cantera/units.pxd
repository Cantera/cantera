#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

import numpy as np
cimport numpy as np

cdef class Units:
    cdef CxxUnits units
    @staticmethod
    cdef copy(CxxUnits)

cdef class UnitSystem:
    cdef CxxUnitSystem unitsystem
