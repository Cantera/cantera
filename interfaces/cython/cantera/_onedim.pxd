#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .base cimport *
from .kinetics cimport *
from .func1 cimport *
from .thermo cimport *

import numpy as np
cimport numpy as np

cdef class Domain1D:
    cdef CxxDomain1D* domain
    cdef _SolutionBase gas
    cdef object _weakref_proxy
    cdef public pybool have_user_tolerances

cdef class Boundary1D(Domain1D):
    cdef CxxBoundary1D* boundary

cdef class Inlet1D(Boundary1D):
    cdef CxxInlet1D* inlet

cdef class Outlet1D(Boundary1D):
    cdef CxxOutlet1D* outlet

cdef class OutletReservoir1D(Boundary1D):
    cdef CxxOutletRes1D* outlet

cdef class SymmetryPlane1D(Boundary1D):
    cdef CxxSymm1D* symm

cdef class Surface1D(Boundary1D):
    cdef CxxSurf1D* surf

cdef class ReactingSurface1D(Boundary1D):
    cdef CxxReactingSurf1D* surf
    cdef public Kinetics surface

cdef class _FlowBase(Domain1D):
    cdef CxxStFlow* flow

cdef class IdealGasFlow(_FlowBase):
    pass

cdef class FreeFlow(IdealGasFlow):
    pass

cdef class IonFlow(_FlowBase):
    pass

cdef class AxisymmetricStagnationFlow(IdealGasFlow):
    pass

cdef class Sim1D:
    cdef CxxSim1D* sim
    cdef readonly object domains
    cdef object _initialized
    cdef object _initial_guess_args
    cdef object _initial_guess_kwargs
    cdef public Func1 _interrupt
    cdef public Func1 _time_step_callback
    cdef public Func1 _steady_callback

cdef CxxIdealGasPhase* getIdealGasPhase(ThermoPhase phase) except *
