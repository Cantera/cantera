# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

import numpy as np
cimport numpy as np

from .ctcxx cimport *

cdef extern from "cantera/thermo/ThermoFactory.h" namespace "Cantera":
    cdef cppclass CxxThermoPhase "Cantera::ThermoPhase"
    cdef shared_ptr[CxxThermoPhase] newThermo(string) except +translate_exception


cdef extern from "cantera/kinetics/KineticsFactory.h" namespace "Cantera":
    cdef cppclass CxxKinetics "Cantera::Kinetics"
    cdef shared_ptr[CxxKinetics] newKinetics (string) except +translate_exception


cdef extern from "cantera/transport/TransportFactory.h" namespace "Cantera":
    cdef cppclass CxxTransport "Cantera::Transport"
    cdef shared_ptr[CxxTransport] newTransport(CxxThermoPhase*, string) except +translate_exception


cdef extern from "cantera/base/Interface.h" namespace "Cantera":
    cdef shared_ptr[CxxSolution] newInterface(
        string, string, vector[string]) except +translate_exception


cdef extern from "cantera/base/Solution.h" namespace "Cantera":
    cdef cppclass CxxKinetics "Cantera::Kinetics"
    cdef cppclass CxxTransport "Cantera::Transport"
    cdef cppclass CxxAnyMap "Cantera::AnyMap"
    cdef cppclass CxxSolution "Cantera::Solution":
        CxxSolution()
        string name()
        void setName(string)
        string source()
        void setSource(string)
        CxxAnyMap& header() except +translate_exception
        shared_ptr[CxxThermoPhase] thermo()
        void setThermo(shared_ptr[CxxThermoPhase])
        shared_ptr[CxxKinetics] kinetics()
        void setKinetics(shared_ptr[CxxKinetics])
        shared_ptr[CxxTransport] transport()
        void setTransport(shared_ptr[CxxTransport])
        void setTransportModel(const string&) except +translate_exception
        CxxAnyMap parameters(cbool) except +translate_exception
        size_t nAdjacent()
        shared_ptr[CxxSolution] adjacent(size_t)

    cdef shared_ptr[CxxSolution] CxxNewSolution "Cantera::Solution::create" ()
    cdef shared_ptr[CxxSolution] newSolution (
        string, string, string, vector[shared_ptr[CxxSolution]]) except +translate_exception
    cdef shared_ptr[CxxSolution] newSolution (
        CxxAnyMap&, CxxAnyMap&, string, vector[shared_ptr[CxxSolution]]) except +translate_exception


ctypedef void (*transportMethod1d)(CxxTransport*, double*) except +translate_exception
ctypedef void (*transportMethod2d)(CxxTransport*, size_t, double*) except +translate_exception
ctypedef void (*transportPolyMethod1i)(CxxTransport*, size_t, double*) except +translate_exception
ctypedef void (*transportPolyMethod2i)(CxxTransport*, size_t, size_t, double*) except +translate_exception

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
