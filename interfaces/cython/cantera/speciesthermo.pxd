# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *

cdef extern from "cantera/thermo/SpeciesThermoInterpType.h":
    cdef cppclass CxxSpeciesThermo "Cantera::SpeciesThermoInterpType":
        CxxSpeciesThermo()
        int reportType()
        void updatePropertiesTemp(double, double*, double*, double*) except +translate_exception
        double minTemp()
        double maxTemp()
        double refPressure()
        void reportParameters(size_t&, int&, double&, double&, double&, double* const) except +translate_exception
        int nCoeffs() except +translate_exception
        CxxAnyMap parameters(cbool) except +translate_exception
        CxxAnyMap& input()


cdef extern from "cantera/thermo/SpeciesThermoFactory.h":
    cdef CxxSpeciesThermo* CxxNewSpeciesThermo "Cantera::newSpeciesThermoInterpType"\
        (int, double, double, double, double*) except +translate_exception


cdef class SpeciesThermo:
    cdef shared_ptr[CxxSpeciesThermo] _spthermo
    cdef CxxSpeciesThermo* spthermo
    cdef _assign(self, shared_ptr[CxxSpeciesThermo] other)

cdef wrapSpeciesThermo(shared_ptr[CxxSpeciesThermo] spthermo)
