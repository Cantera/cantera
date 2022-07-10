# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

cdef extern from "cantera/equil/MultiPhase.h" namespace "Cantera":
    cdef cppclass CxxThermoPhase "Cantera::ThermoPhase"
    cdef cppclass CxxMultiPhase "Cantera::MultiPhase":
        CxxMultiPhase()
        void addPhase(CxxThermoPhase*, double) except +translate_exception
        void init() except +translate_exception
        void updatePhases() except +translate_exception

        void equilibrate(string, string, double, int, int, int, int) except +translate_exception

        size_t nSpecies()
        size_t nElements()
        size_t nPhases()
        size_t elementIndex(string) except +translate_exception
        size_t speciesIndex(size_t, size_t) except +translate_exception
        string speciesName(size_t) except +translate_exception
        double nAtoms(size_t, size_t) except +translate_exception

        double phaseMoles(size_t) except +translate_exception
        void setPhaseMoles(size_t, double) except +translate_exception
        void setMoles(double*) except +translate_exception
        void setMolesByName(string) except +translate_exception

        double speciesMoles(size_t) except +translate_exception
        double elementMoles(size_t) except +translate_exception

        void setTemperature(double) except +translate_exception
        double temperature()
        void setPressure(double) except +translate_exception
        double pressure()

        double minTemp() except +translate_exception
        double maxTemp() except +translate_exception
        double charge() except +translate_exception
        double phaseCharge(size_t) except +translate_exception
        void getChemPotentials(double*) except +translate_exception
        double enthalpy() except +translate_exception
        double entropy() except +translate_exception
        double gibbs() except +translate_exception
        double cp() except +translate_exception
        double volume() except +translate_exception


cdef class Mixture:
    cdef CxxMultiPhase* mix
    cdef list _phases
    cdef object _weakref_proxy
    cpdef int element_index(self, element) except *
