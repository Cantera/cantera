#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .kinetics cimport *
from .func1 cimport *

cdef class ReactorBase:
    cdef CxxReactorBase* rbase
    cdef object _thermo
    cdef list _inlets
    cdef list _outlets
    cdef list _walls
    cdef list _surfaces
    cdef object _weakref_proxy

cdef class Reactor(ReactorBase):
    cdef CxxReactor* reactor
    cdef object _kinetics

cdef class Reservoir(ReactorBase):
    pass

cdef class ConstPressureReactor(Reactor):
    pass

cdef class IdealGasReactor(Reactor):
    pass

cdef class IdealGasConstPressureReactor(Reactor):
    pass

cdef class FlowReactor(Reactor):
    pass

cdef class ExtensibleReactor(Reactor):
    cdef CxxReactorAccessor* accessor

cdef class ReactorSurface:
    cdef CxxReactorSurface* surface
    cdef Kinetics _kinetics

cdef class WallBase:
    cdef CxxWallBase* wall
    cdef object _velocity_func
    cdef object _heat_flux_func
    cdef ReactorBase _left_reactor
    cdef ReactorBase _right_reactor
    cdef str name

cdef class Wall(WallBase):
    pass

cdef class FlowDevice:
    cdef CxxFlowDevice* dev
    cdef Func1 _rate_func
    cdef Func1 _time_func
    cdef str name
    cdef ReactorBase _upstream
    cdef ReactorBase _downstream

cdef class MassFlowController(FlowDevice):
    pass

cdef class Valve(FlowDevice):
    pass

cdef class PressureController(FlowDevice):
    pass

cdef class ReactorNet:
    cdef CxxReactorNet net
    cdef list _reactors
