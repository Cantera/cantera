#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .base cimport *

cdef class Species:
    cdef shared_ptr[CxxSpecies] _species
    cdef CxxSpecies* species
    cdef _SolutionBase _phase

    cdef _assign(self, shared_ptr[CxxSpecies] other)

cdef class ThermoPhase(_SolutionBase):
    cdef double _mass_factor(self)
    cdef double _mole_factor(self)
    cpdef int element_index(self, element) except *
    cpdef int species_index(self, species) except *
    cdef np.ndarray _getArray1(self, thermoMethod1d method)
    cdef void _setArray1(self, thermoMethod1d method, values) except *
    cdef CxxPlasmaPhase* plasma
    cdef public object _enable_plasma

cdef class InterfacePhase(ThermoPhase):
    cdef CxxSurfPhase* surf
