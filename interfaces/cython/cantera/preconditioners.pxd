# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *
from .kinetics cimport CxxSparseMatrix

cdef extern from "cantera/numerics/PreconditionerBase.h" namespace "Cantera":
    cdef cppclass CxxPreconditionerBase "Cantera::PreconditionerBase":
        CxxPreconditionerBase()
        string preconditionerSide()
        string type()
        void setPreconditionerSide(string) except +translate_exception

cdef extern from "cantera/numerics/AdaptivePreconditioner.h" namespace "Cantera":
    cdef cppclass CxxAdaptivePreconditioner "Cantera::AdaptivePreconditioner" \
        (CxxPreconditionerBase):
        CxxAdaptivePreconditioner() except +translate_exception
        void setThreshold(double threshold)
        double threshold()
        void setIlutFillFactor(int fillfactor)
        double ilutFillFactor()
        void setIlutDropTol(double droptol)
        double ilutDropTol()
        void printPreconditioner()
        CxxSparseMatrix matrix() except +translate_exception

cdef extern from "cantera/oneD/MultiJac.h" namespace "Cantera":
    cdef cppclass CxxMultiJac "Cantera::MultiJac" (CxxPreconditionerBase):
        CxxMultiJac() except +translate_exception

cdef extern from "cantera/numerics/PreconditionerFactory.h" namespace "Cantera":
    cdef shared_ptr[CxxPreconditionerBase] newPreconditioner(string) except\
         +translate_exception

cdef class PreconditionerBase:
    @staticmethod
    cdef wrap(shared_ptr[CxxPreconditionerBase])
    cdef set_cxx_object(self)
    cdef shared_ptr[CxxPreconditionerBase] pbase

cdef class AdaptivePreconditioner(PreconditionerBase):
    cdef set_cxx_object(self)
    cdef CxxAdaptivePreconditioner* preconditioner

cdef class BandedJacobian(PreconditionerBase):
    cdef set_cxx_object(self)
    cdef CxxMultiJac* preconditioner

