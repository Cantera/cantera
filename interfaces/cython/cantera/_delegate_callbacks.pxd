# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from .ctcxx cimport *
from .func1 cimport *
from .units cimport CxxUnitStack
from .delegator cimport CxxEigenTriplet


# Wrapper functions that call Python delegate methods from C++. Their signatures match
# the second argument of the ``pyOverride`` overloads declared in ``delegator.pxd``, and
# they are passed to ``pyOverride`` from ``delegator.assign_delegates``. They are declared
# here (rather than in ``delegator.pxd``) so that they can be defined in native Cython
# syntax: the C++ reference parameters they require cannot be expressed in the
# pure-Python ``delegator.py``.
cdef void callback_v(PyFuncInfo& funcInfo) noexcept
cdef void callback_v_d(PyFuncInfo& funcInfo, double arg) noexcept
cdef void callback_v_b(PyFuncInfo& funcInfo, cbool arg) noexcept
cdef void callback_v_AMr(PyFuncInfo& funcInfo, CxxAnyMap& arg) noexcept
cdef void callback_v_cAMr_cUSr(PyFuncInfo& funcInfo, const CxxAnyMap& arg1,
                               const CxxUnitStack& arg2) noexcept
cdef void callback_v_csr_vp(PyFuncInfo& funcInfo,
                            const string& arg1, void* obj) noexcept
cdef void callback_v_dp(PyFuncInfo& funcInfo, span[double] arg) noexcept
cdef void callback_v_d_dp(PyFuncInfo& funcInfo, double arg1,
                          span[double] arg2) noexcept
cdef void callback_v_dp_dp_dp(PyFuncInfo& funcInfo,
        span[double] arg1, span[double] arg2, span[double] arg3) noexcept
cdef void callback_v_vETr(PyFuncInfo& funcInfo,
                          vector[CxxEigenTriplet]& trips) noexcept
cdef int callback_d_vp(PyFuncInfo& funcInfo, double& out, void* obj) noexcept
cdef int callback_s_sz(PyFuncInfo& funcInfo, string& out, size_t arg) noexcept
cdef int callback_sz_csr(PyFuncInfo& funcInfo, size_t& out, const string& arg) noexcept
cdef void callback_v_d_dp_dp(PyFuncInfo& funcInfo, double arg1,
                             span[double] arg2, span[double] arg3) noexcept
