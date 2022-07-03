#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .func1 cimport *

cdef class ReactionRate:
    cdef shared_ptr[CxxReactionRate] _rate
    cdef CxxReactionRate* rate
    @staticmethod
    cdef wrap(shared_ptr[CxxReactionRate])
    cdef set_cxx_object(self)

cdef class ArrheniusRateBase(ReactionRate):
    cdef CxxArrheniusBase* base

cdef class FalloffRate(ReactionRate):
    cdef CxxFalloffRate* falloff
    cdef set_cxx_object(self)

cdef class CustomRate(ReactionRate):
    cdef CxxCustomFunc1Rate* cxx_object(self)
    cdef Func1 _rate_func

cdef class InterfaceRateBase(ArrheniusRateBase):
    cdef CxxInterfaceRateBase* interface

cdef class StickRateBase(InterfaceRateBase):
    cdef CxxStickingCoverage* stick

cdef class Reaction:
    cdef shared_ptr[CxxReaction] _reaction
    cdef CxxReaction* reaction
    @staticmethod
    cdef wrap(shared_ptr[CxxReaction])

cdef class CustomReaction(Reaction):
    cdef CustomRate _rate

cdef class Arrhenius:
    cdef CxxArrheniusRate* base
    cdef cbool own_rate
    cdef Reaction reaction # parent reaction, to prevent garbage collection
    @staticmethod
    cdef wrap(CxxArrheniusRate*)

cdef class Falloff:
    cdef shared_ptr[CxxFalloff] _falloff
    cdef CxxFalloff* falloff
