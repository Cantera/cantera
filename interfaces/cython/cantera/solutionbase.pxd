# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

import numpy as np
cimport numpy as np

from .ctcxx cimport *
from .delegator cimport CxxExternalHandle


cdef extern from "cantera/thermo/ThermoFactory.h" namespace "Cantera":
    cdef cppclass CxxThermoPhase "Cantera::ThermoPhase"
    cdef shared_ptr[CxxThermoPhase] newThermoModel(string) except +translate_exception


cdef extern from "cantera/kinetics/KineticsFactory.h" namespace "Cantera":
    cdef cppclass CxxKinetics "Cantera::Kinetics"
    cdef shared_ptr[CxxKinetics] newKinetics (string) except +translate_exception


cdef extern from "cantera/transport/TransportFactory.h" namespace "Cantera":
    cdef cppclass CxxTransport "Cantera::Transport"
    cdef shared_ptr[CxxTransport] newTransport(shared_ptr[CxxThermoPhase], string) except +translate_exception


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
        void holdExternalHandle(string&, shared_ptr[CxxExternalHandle])
        void registerChangedCallback(void*, function[void()])
        void removeChangedCallback(void*)

    cdef shared_ptr[CxxSolution] CxxNewSolution "Cantera::Solution::create" ()
    cdef shared_ptr[CxxSolution] newSolution (
        string, string, string, vector[shared_ptr[CxxSolution]]) except +translate_exception
    cdef shared_ptr[CxxSolution] newSolution (
        CxxAnyMap&, CxxAnyMap&, string, vector[shared_ptr[CxxSolution]]) except +translate_exception


cdef extern from "cantera/base/SolutionArray.h" namespace "Cantera":
    cdef cppclass CxxSolutionArray "Cantera::SolutionArray":
        shared_ptr[CxxSolutionArray] share(vector[int]&) except +translate_exception
        void reset() except +translate_exception
        int size()
        void resize(int) except +translate_exception
        vector[long int] apiShape() except +translate_exception
        void setApiShape(vector[long int]&) except +translate_exception
        int apiNdim()
        string info(vector[string]&, int, int) except +translate_exception
        CxxAnyMap meta()
        void setMeta(CxxAnyMap&)
        vector[string] componentNames() except +translate_exception
        cbool hasComponent(string&)
        CxxAnyValue getComponent(string&) except +translate_exception
        void setComponent(string&, CxxAnyValue&) except +translate_exception
        void setLoc(int) except +translate_exception
        void updateState(int) except +translate_exception
        vector[double] getState(int) except +translate_exception
        void setState(int, vector[double]&) except +translate_exception
        vector[string] listExtra()
        cbool hasExtra(string&)
        void addExtra(string&, cbool) except +translate_exception
        CxxAnyMap getAuxiliary(int) except +translate_exception
        void setAuxiliary(int, CxxAnyMap&) except +translate_exception
        void append(vector[double]&, CxxAnyMap&) except +translate_exception
        void save(string&, string&, string&, string&, cbool, int, string&) except +translate_exception
        CxxAnyMap restore(string&, string&, string&) except +translate_exception

    cdef shared_ptr[CxxSolutionArray] CxxNewSolutionArray "Cantera::SolutionArray::create" (
        shared_ptr[CxxSolution], int, CxxAnyMap&) except +translate_exception


ctypedef void (*transportMethod1d)(CxxTransport*, double*) except +translate_exception
ctypedef void (*transportMethod2d)(CxxTransport*, size_t, double*) except +translate_exception
ctypedef void (*transportPolyMethod1i)(CxxTransport*, size_t, double*) except +translate_exception
ctypedef void (*transportPolyMethod2i)(CxxTransport*, size_t, size_t, double*) except +translate_exception

cdef _assign_Solution(_SolutionBase soln, shared_ptr[CxxSolution] cxx_soln,
                      pybool reset_adjacent, pybool weak=?)
cdef object _wrap_Solution(shared_ptr[CxxSolution] cxx_soln)

cdef class _SolutionBase:
    cdef shared_ptr[CxxSolution] _base
    cdef weak_ptr[CxxSolution] weak_base
    cdef CxxSolution* base
    cdef CxxThermoPhase* thermo
    cdef CxxKinetics* kinetics
    cdef CxxTransport* transport
    cdef int thermo_basis
    cdef np.ndarray _selected_species
    cdef object parent
    cdef object _adjacent
    cdef object _soln_changed_callback
    cdef public object _references

cdef class SolutionArrayBase:
    cdef shared_ptr[CxxSolutionArray] _base
    cdef CxxSolutionArray* base
    cdef public object _weakref_proxy
