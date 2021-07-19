# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .func1 cimport *

cdef extern from "cantera/numerics/PreconditionerBase.h" namespace "Cantera":
    cdef cppclass CxxPreconditionerBase "Cantera::PreconditionerBase":
        CxxPreconditionerBase()

cdef extern from "cantera/numerics/AdaptivePreconditioner.h" namespace "Cantera":
    cdef cppclass CxxAdaptivePreconditioner "Cantera::AdaptivePreconditioner" (CxxPreconditionerBase):
        CxxAdaptivePreconditioner() except +
        void setThreshold(double threshold)
        double threshold()
        void setIlutFillFactor(int fillfactor)
        double ilutFillFactor()
        void setIlutDropTol(double droptol)
        double ilutDropTol()
        void printPreconditioner()

cdef extern from "cantera/numerics/PreconditionerFactory.h" namespace "Cantera":
    cdef CxxPreconditionerBase* newPreconditioner(string) except +translate_exception

cdef class PreconditionerBase:
    cdef CxxPreconditionerBase* pbase

cdef class AdaptivePreconditioner(PreconditionerBase):
    cdef CxxAdaptivePreconditioner* preconditioner
