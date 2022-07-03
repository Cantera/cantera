#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .kinetics cimport *

cdef class ReactionPathDiagram:
    cdef CxxReactionPathDiagram diagram
    cdef CxxReactionPathBuilder builder
    cdef Kinetics kinetics
    cdef str element
    cdef pybool built
    cdef CxxStringStream* _log
