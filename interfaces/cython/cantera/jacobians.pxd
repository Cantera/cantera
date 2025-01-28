# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *
from .kinetics cimport CxxSparseMatrix

cdef extern from "cantera/numerics/SystemJacobian.h" namespace "Cantera":
    cdef cppclass CxxSystemJacobian "Cantera::SystemJacobian":
        CxxSystemJacobian()
        string preconditionerSide()
        string type()
        void setPreconditionerSide(string) except +translate_exception

cdef extern from "cantera/numerics/EigenSparseJacobian.h" namespace "Cantera":
    cdef cppclass CxxEigenSparseJacobian "Cantera::EigenSparseJacobian" \
        (CxxSystemJacobian):
        CxxEigenSparseJacobian() except +translate_exception
        void printPreconditioner()
        CxxSparseMatrix matrix() except +translate_exception
        CxxSparseMatrix jacobian() except +translate_exception

cdef extern from "cantera/numerics/AdaptivePreconditioner.h" namespace "Cantera":
    cdef cppclass CxxAdaptivePreconditioner "Cantera::AdaptivePreconditioner" \
        (CxxEigenSparseJacobian):
        CxxAdaptivePreconditioner() except +translate_exception
        void setThreshold(double threshold)
        double threshold()
        void setIlutFillFactor(int fillfactor)
        double ilutFillFactor()
        void setIlutDropTol(double droptol)
        double ilutDropTol()

cdef extern from "cantera/oneD/MultiJac.h" namespace "Cantera":
    cdef cppclass CxxMultiJac "Cantera::MultiJac" (CxxSystemJacobian):
        CxxMultiJac() except +translate_exception

cdef extern from "cantera/numerics/SystemJacobianFactory.h" namespace "Cantera":
    cdef shared_ptr[CxxSystemJacobian] newSystemJacobian(string) except\
         +translate_exception

cdef class SystemJacobian:
    @staticmethod
    cdef wrap(shared_ptr[CxxSystemJacobian])
    cdef set_cxx_object(self)
    cdef shared_ptr[CxxSystemJacobian] _base

cdef class EigenSparseJacobian(SystemJacobian):
    cdef set_cxx_object(self)
    cdef CxxEigenSparseJacobian* sparse_jac

cdef class EigenSparseDirectJacobian(EigenSparseJacobian):
    pass

cdef class AdaptivePreconditioner(EigenSparseJacobian):
    cdef set_cxx_object(self)
    cdef CxxAdaptivePreconditioner* adaptive

cdef class BandedJacobian(SystemJacobian):
    cdef set_cxx_object(self)
    cdef CxxMultiJac* band_jac
