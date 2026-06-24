# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# This module holds the parts of the delegate machinery that must be implemented in
# native Cython syntax because they rely on C++ reference parameters, which have no
# spelling in Cython's pure-Python syntax:
#
# - The ``callback_*`` wrapper functions, which call Python delegate methods from C++.
#   Their signatures match the function-pointer arguments of the ``pyOverride`` overloads
#   (see ``funcWrapper.h`` and ``delegator.pxd``) and are passed to ``pyOverride`` from
#   ``delegator.assign_delegates``.
# - The ``ct_*`` functions, declared ``public`` so that they can be called from
#   ``PythonExtensionManager.cpp`` (via the generated ``_delegate_callbacks.h`` header).
#
# The higher-level Python logic (the ``extension`` decorator and ``assign_delegates``)
# lives in the pure-Python ``delegator.py`` module.
#
# ## Implementation for each delegated function type
#
# Besides the C++ functions implemented in the `Delegator` class, each delegated
# function type requires several additional components in Cython:
# - A declaration for the overload of `Delegator::setDelegate` in `delegator.pxd`.
# - A declaration of the `pyOverride` function in `delegator.pxd`. Although the
#   `pyOverride` function is templated and does not need to be modified for each type,
#   each version of this function must be separately declared because Cython doesn't
#   understand variadic templates. The signature that needs to be declared has:
#   - as its return value, a `std::function` object where the function type matches the
#     function type needed by `Delegator::setDelegate` (which is different from the
#     original member function for functions that take array arguments or have return
#     values)
#   - as its first argument, a `PyObject*`
#   - as its second argument, a function pointer corresponding to the return value but
#     with a `PyFuncInfo&` added as the first argument
# - A `cdef ` "callback" wrapper function implemented in this module, whose signature
#   matches the second argument to `pyOverride` as described above, and implements the
#   behavior described below.
# - A case in the `elif callback == ...` tree in the `assign_delegates` function (in
#   `delegator.py`), which creates the C++ function object wrapper for a Python method
#   and calls `Delegator::setDelegate`

# ## Wrapper functions for calling Python functions from C++
#
# The following functions are used along with the `pyOverride` function to wrap a Python
# function inside a C++ `std::function` object that can be used with the C++ `Delegator`
# class (see `Delegator.h`). The functions here are responsible for the following
# behaviors:
# - Mapping array arguments from C pointers to Cython "memoryviewslice" objects
# - Catching and stashing any exceptions that might be thrown by the wrapped Python
#   function. These will be translated to C++ exceptions and be re-thrown by the
#   `pyOverride`-generated wrapper
# - Translating the return value semantics of the Python function into the form
#   required by the `Delegator` class, where the Python function's return value is
#   an "output" argument of the function, and the actual return value (int) is used
#   to indicate whether the Python function returned a value or None.
# - Converting between C++ and Python strings
#
# The callback function names use a naming scheme based on the function signature of
# the corresponding C++ member function. After the prefix `callback_` is a notation
# of the C++ member function's return type, followed by an underscore, then the
# notations for each argument, separated by underscores. The following shorthand is
# used for different return / argument types:
# - `v` for `void`
# - `b` for `bool`
# - `d` for `double`
# - `s` for `std::string`
# - `sz` for `size_t`
# - prefix `c` for `const` arguments
# - suffix `r` for reference arguments
# - suffix `p` for pointer arguments
#
# See `funcWrapper.h` for the definition of the `PyFuncInfo` class and the `pyOverride`
# function.

import sys as _sys
import importlib as _importlib

from libc.stdlib cimport malloc
from libc.string cimport strcpy

from ._utils cimport stringify, pystr, anymap_to_py, py_to_anymap
from .units cimport Units, UnitStack
from .reaction cimport (ExtensibleRate, ExtensibleRateData, CxxReaction,
    CxxReactionRateDelegator, CxxReactionDataDelegator)
from .solutionbase cimport CxxSolution, _assign_Solution
from .delegator cimport CxxExtensionManager, CxxExtensionManagerFactory

from cpython.object cimport PyTypeObject, traverseproc, visitproc, inquiry

cdef extern from "cantera/extensions/PythonExtensionManager.h" namespace "Cantera":
    cdef cppclass CxxPythonExtensionManager "Cantera::PythonExtensionManager":
        @staticmethod
        void registerSelf()


# Wrapper for functions of type void()
cdef void callback_v(PyFuncInfo& funcInfo) noexcept:
    try:
        (<object>funcInfo.func())()
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double)
cdef void callback_v_d(PyFuncInfo& funcInfo, double arg) noexcept:
    try:
        (<object>funcInfo.func())(arg)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(bool)
cdef void callback_v_b(PyFuncInfo& funcInfo, cbool arg) noexcept:
    try:
        (<object>funcInfo.func())(arg)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(AnyMap&)
cdef void callback_v_AMr(PyFuncInfo& funcInfo, CxxAnyMap& arg) noexcept:
    pyArg = anymap_to_py(<CxxAnyMap&>arg)  # cast away constness
    try:
        (<object>funcInfo.func())(pyArg)
        # return updated AnyMap to C++. Odd syntax is a workaround for Cython's
        # unwillingness to assign to a reference
        (&arg)[0] = py_to_anymap(pyArg)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(const AnyMap&, const UnitStack&)
cdef void callback_v_cAMr_cUSr(PyFuncInfo& funcInfo, const CxxAnyMap& arg1,
                               const CxxUnitStack& arg2) noexcept:

    pyArg1 = anymap_to_py(<CxxAnyMap&>arg1)  # cast away constness
    pyArg2 = UnitStack.copy(arg2)
    try:
        (<object>funcInfo.func())(pyArg1, pyArg2)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(const string&, void*)
cdef void callback_v_csr_vp(PyFuncInfo& funcInfo,
                            const string& arg1, void* obj) noexcept:
    try:
        (<object>funcInfo.func())(pystr(arg1), <object>obj)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(span<double>)
cdef void callback_v_dp(PyFuncInfo& funcInfo, span[double] arg) noexcept:
    cdef double[:] view = <double[:arg.size()]>arg.data() if arg.size() else None

    try:
        (<object>funcInfo.func())(view)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double, span<double>)
cdef void callback_v_d_dp(PyFuncInfo& funcInfo, double arg1,
                          span[double] arg2) noexcept:
    cdef double[:] view = <double[:arg2.size()]>arg2.data() if arg2.size() else None

    try:
        (<object>funcInfo.func())(arg1, view)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(span<double>, span<double>, span<double>)
cdef void callback_v_dp_dp_dp(PyFuncInfo& funcInfo,
        span[double] arg1, span[double] arg2, span[double] arg3) noexcept:
    cdef double[:] view1 = <double[:arg1.size()]>arg1.data() if arg1.size() else None
    cdef double[:] view2 = <double[:arg2.size()]>arg2.data() if arg2.size() else None
    cdef double[:] view3 = <double[:arg3.size()]>arg3.data() if arg3.size() else None
    try:
        (<object>funcInfo.func())(view1, view2, view3)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(SparseTriplets&)
cdef void callback_v_vETr(PyFuncInfo& funcInfo,
                          vector[CxxEigenTriplet]& trips) noexcept:
    try:
        pyTrips = []
        ret = (<object>funcInfo.func())(pyTrips)
        if ret is not None:
            pyTrips.extend(ret)
        for row, col, value in pyTrips:
            trips.push_back(CxxEigenTriplet(<size_t>row, <size_t>col, <double>value))
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type double(void*)
cdef int callback_d_vp(PyFuncInfo& funcInfo, double& out, void* obj) noexcept:
    try:
        ret = (<object>funcInfo.func())(<object>obj)
        if ret is None:
            return 0
        else:
            (&out)[0] = ret
            return 1
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1

# Wrapper for functions of type string(size_t)
cdef int callback_s_sz(PyFuncInfo& funcInfo, string& out, size_t arg) noexcept:
    try:
        ret = (<object>funcInfo.func())(arg)
        if ret is None:
            return 0
        else:
            (&out)[0] = stringify(ret)
            return 1
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1

# Wrapper for functions of type size_t(string&)
cdef int callback_sz_csr(PyFuncInfo& funcInfo, size_t& out, const string& arg) noexcept:
    try:
        ret = (<object>funcInfo.func())(pystr(arg))
        if ret is None:
            return 0
        else:
            (&out)[0] = ret
            return 1
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1

# Wrapper for functions of type void(double, span<double>, span<double>)
cdef void callback_v_d_dp_dp(PyFuncInfo& funcInfo, double arg1,
                             span[double] arg2, span[double] arg3) noexcept:
    cdef double[:] view1 = <double[:arg2.size()]>arg2.data() if arg2.size() else None
    cdef double[:] view2 = <double[:arg3.size()]>arg3.data() if arg3.size() else None

    try:
        (<object>funcInfo.func())(arg1, view1, view2)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)


# Specifications for ReactionRate delegators that have not yet been registered with
# ReactionRateFactory. This list is read by PythonExtensionManager::registerRateBuilders
# and then cleared.
_rate_delegators = []
_rate_data_delegators = []


cdef public char* ct_getExceptionString(object exType, object exValue, object exTraceback):
    import traceback
    result = str(exValue) + "\n\n"
    result += "".join(traceback.format_exception(exType, exValue, exTraceback))
    tmp = bytes(result.encode())
    cdef char* c_string = <char*> malloc((len(tmp) + 1) * sizeof(char))
    strcpy(c_string, tmp)
    return c_string


cdef public object ct_newPythonExtensibleRate(CxxReactionRateDelegator* delegator,
                                              const string& module_name,
                                              const string& class_name):

    mod = _importlib.import_module(module_name.decode())
    cdef ExtensibleRate rate = getattr(mod, class_name.decode())(init=False)
    rate.set_cxx_object(delegator)
    return rate


cdef public object ct_newPythonExtensibleRateData(CxxReactionDataDelegator* delegator,
                                                  const string& module_name,
                                                  const string& class_name):

    mod = _importlib.import_module(module_name.decode())
    cdef ExtensibleRateData data = getattr(mod, class_name.decode())()
    data.set_cxx_object(delegator)
    return data


cdef public ct_registerReactionDelegators():
    cdef shared_ptr[CxxExtensionManager] mgr = (
        CxxExtensionManagerFactory.build(stringify("python")))

    for module, cls, name in _rate_delegators:
        mgr.get().registerRateBuilder(stringify(module), stringify(cls), stringify(name))

    _rate_delegators.clear()

    for module, cls, name in _rate_data_delegators:
        mgr.get().registerRateDataBuilder(stringify(module), stringify(cls), stringify(name))

    _rate_data_delegators.clear()


cdef public object ct_wrapSolution(shared_ptr[CxxSolution] soln):
    from .composite import Solution
    pySoln = Solution(init=False)
    _assign_Solution(pySoln, soln, False, weak=True)
    return pySoln

# The pair of (ExtensibleRate, ReactionRateDelegator) objects both hold owned references
# to one another. To allow the Python garbage collector to break this reference cycle,
# we need to implement custom behavior to detect when the only object referring the
# ReactionRateDelegator is the ExtensibleRate.
# Implementation roughly follows from https://github.com/mdavidsaver/cython-c--demo

# Capture the original implementations of the tp_traverse and tp_clear methods
cdef PyTypeObject* extensibleRate_t = <PyTypeObject*>ExtensibleRate
cdef traverseproc extensibleRate_base_traverse = extensibleRate_t.tp_traverse
cdef inquiry extensibleRate_base_clear = extensibleRate_t.tp_clear
assert extensibleRate_base_traverse != NULL
assert extensibleRate_base_clear != NULL

cdef int traverse_ExtensibleRate(PyObject* raw, visitproc visit, void* arg) except -1:
    cdef ExtensibleRate self = <ExtensibleRate>raw
    cdef int ret = 0
    # If self._rate.use_count() is 1, this ExtensibleRate rate is the only object
    # referencing self._rate. To let the GC see the cycle where self._rate references
    # self, we tell it to visit self. If self._rate.use_count() is more than one, there
    # are other C++ objects referring to self._rate, and we don't want the GC to see a
    # cycle, so we skip visiting self.
    if self._rate.use_count() == 1:
        ret = visit(<PyObject*>self, arg)
    if ret:
        return ret
    # Call the original traverser to deal with all other members
    ret = extensibleRate_base_traverse(raw, visit, arg)
    return ret

cdef int clear_ExtensibleRate(object obj) except -1:
    cdef ExtensibleRate self = <ExtensibleRate>obj
    # If the GC has called this method, this ExtensibleRate is the only object holding
    # a reference to the ReactionRateDelegator, and resetting the shared_ptr will delete
    # both that object and its reference to the ExtensibleRate, allowing the ref count
    # for the ExtensibleRate to go to zero.
    self._rate.reset()
    return extensibleRate_base_clear(obj)

# Assign the augmented garbage collector functions
extensibleRate_t.tp_traverse = traverse_ExtensibleRate
extensibleRate_t.tp_clear = clear_ExtensibleRate


CxxPythonExtensionManager.registerSelf()
