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
        CxxUnits(CxxUnits)
        CxxUnits(string, cbool) except +translate_exception
        string str()
        double factor()
        double dimension(string) except +translate_exception

    cdef cppclass CxxUnitSystem "Cantera::UnitSystem":
        CxxUnitSystem()
        stdmap[string, string] defaults()
        void setDefaults(stdmap[string, string]&) except +translate_exception

    cdef cppclass CxxUnitStack "Cantera::UnitStack":
        CxxUnitStack()
        CxxUnits product()


cdef class Units:
    cdef CxxUnits units
    @staticmethod
    cdef copy(CxxUnits)

cdef class UnitSystem:
    cdef CxxUnitSystem unitsystem
