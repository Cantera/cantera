# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map as stdmap
from libcpp.cast cimport dynamic_cast
from libcpp.pair cimport pair
from libcpp cimport bool as cbool
from libcpp.functional cimport function
from libcpp.memory cimport shared_ptr, weak_ptr, dynamic_pointer_cast
from cpython cimport bool as pybool
from cpython.ref cimport PyObject
from cython.operator cimport dereference as deref, preincrement as inc

ctypedef stdmap[string,double] Composition

import numpy as np
cimport numpy as np

cdef extern from "cantera/cython/funcWrapper.h":
    cdef int translate_exception()

    cdef cppclass CxxAnyMap "Cantera::AnyMap"
    cdef cppclass CxxAnyValue "Cantera::AnyValue"
    cdef cppclass CxxUnits "Cantera::Units"
    cdef cppclass CxxUnitSystem "Cantera::UnitSystem"
