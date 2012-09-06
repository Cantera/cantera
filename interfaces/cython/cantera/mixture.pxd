from solution cimport *

cdef extern from "cantera/equil/MultiPhase.h" namespace "Cantera":
    cdef cppclass CxxMultiPhase "Cantera::MultiPhase":
        CxxMultiPhase()
        void addPhase(CxxThermoPhase*, double) except +
        void init() except +
        double nSpecies()
        void setTemperature(double)
        double temperature()
        void setPressure(double)
        double pressure()

cdef class Mixture:
    cdef CxxMultiPhase* mix
    cdef list _phases
