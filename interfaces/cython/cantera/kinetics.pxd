#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .base cimport *

cdef class Kinetics(_SolutionBase):
    pass

cdef class InterfaceKinetics(Kinetics):
    pass

cdef np.ndarray get_species_array(Kinetics kin, kineticsMethod1d method)
cdef np.ndarray get_reaction_array(Kinetics kin, kineticsMethod1d method)
cdef np.ndarray get_dense(Kinetics kin, kineticsMethodSparse method)
cdef tuple get_sparse(Kinetics kin, kineticsMethodSparse method)
