#cython: language_level=3
#distutils: language = c++

import numpy as np
cimport numpy as np

from .ctcxx cimport *

cdef _assign_Solution(_SolutionBase soln, shared_ptr[CxxSolution] cxx_soln,
                      pybool reset_adjacent)
cdef object _wrap_Solution(shared_ptr[CxxSolution] cxx_soln)

cdef class _SolutionBase:
    cdef shared_ptr[CxxSolution] _base
    cdef CxxSolution* base
    cdef CxxThermoPhase* thermo
    cdef CxxKinetics* kinetics
    cdef CxxTransport* transport
    cdef int thermo_basis
    cdef np.ndarray _selected_species
    cdef object parent
    cdef object _adjacent
    cdef public object _references
