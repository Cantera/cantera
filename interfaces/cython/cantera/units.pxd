# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *

import numpy as np
cimport numpy as np

cdef extern from "cantera/base/Units.h" namespace "Cantera":
    cdef cppclass CxxUnits "Cantera::Units":
        CxxUnits()
        CxxUnits(CxxUnits&)
        CxxUnits(string, cbool) except +translate_exception
        string str()
        double factor()
        double dimension(string) except +translate_exception

    cdef cppclass CxxUnitSystem "Cantera::UnitSystem":
        CxxUnitSystem()
        stdmap[string, string] defaults()
        void setDefaults(stdmap[string, string]&) except +translate_exception
        double convert(CxxAnyValue&, string&) except +translate_exception
        double convert(CxxAnyValue&, CxxUnits&) except +translate_exception
        double convertTo(double, string&) except +translate_exception
        double convertTo(double, CxxUnits&) except +translate_exception
        double convertRateCoeff(CxxAnyValue&, CxxUnits&) except +translate_exception
        double convertActivationEnergy(CxxAnyValue&, string&) except +translate_exception
        double convertActivationEnergyTo(double, string&) except +translate_exception
        double convertActivationEnergyTo(double, CxxUnits&) except +translate_exception

    cdef cppclass CxxUnitStack "Cantera::UnitStack":
        CxxUnitStack()
        CxxUnitStack(CxxUnits&)
        CxxUnitStack(CxxUnitStack&)
        CxxUnits product()
        void join(double) except +translate_exception


cdef class Units:
    cdef CxxUnits units
    @staticmethod
    cdef Units copy(CxxUnits)

cdef class UnitStack:
    cdef CxxUnitStack stack
    @staticmethod
    cdef UnitStack copy(const CxxUnitStack&)


cdef class UnitSystem:
    cdef _set_unitSystem(self, shared_ptr[CxxUnitSystem] units)
    cdef shared_ptr[CxxUnitSystem] _unitsystem
    cdef CxxUnitSystem* unitsystem
