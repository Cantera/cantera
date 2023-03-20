# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *
from .units cimport *

cdef extern from "cantera/base/YamlWriter.h" namespace "Cantera":
    cdef cppclass CxxSolution "Cantera::Solution"
    cdef cppclass CxxYamlWriter "Cantera::YamlWriter":
        CxxYamlWriter()
        void setHeader(CxxAnyMap) except +translate_exception
        void addPhase(shared_ptr[CxxSolution]) except +translate_exception
        string toYamlString() except +translate_exception
        void toYamlFile(string&) except +translate_exception
        void setPrecision(int)
        void skipUserDefined(cbool)
        void setUnitSystem(CxxUnitSystem&) except +translate_exception


cdef class YamlWriter:
    cdef shared_ptr[CxxYamlWriter] _writer
    cdef CxxYamlWriter* writer
    @staticmethod
    cdef shared_ptr[CxxUnitSystem] _get_unitsystem(UnitSystem units)
