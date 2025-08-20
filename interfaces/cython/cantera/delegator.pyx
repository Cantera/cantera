# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import inspect as _inspect
import sys as _sys
import importlib as _importlib

from libc.stdlib cimport malloc
from libc.string cimport strcpy

from ._utils import CanteraError
from ._utils cimport stringify, pystr, anymap_to_py, py_to_anymap
from .units cimport Units, UnitStack
# from .reaction import ExtensibleRate, ExtensibleRateData
from .reaction cimport (ExtensibleRate, ExtensibleRateData, CxxReaction,
    CxxReactionRateDelegator, CxxReactionDataDelegator)
from .solutionbase cimport CxxSolution, _assign_Solution
from cython.operator import dereference as deref

from cpython.object cimport PyTypeObject, traverseproc, visitproc, inquiry

# ## Implementation for each delegated function type
#
# Besides the C++ functions implemented in the `Delegator` class, each delegated
# function type requires several additional components in Cython:
# - A declaration for the overload of `Delegator::setDelegate` in `_cantera.pxd`.
# - A declaration of the `pyOverride` function in `_cantera.pxd`. Although the
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
# - A `cdef ` "callback" wrapper function implemented in `delegator.pyx`, whose
#   signature matches the second argument to `pyOverride` as described above, and
#   implements the behavior described below.
# - A case in the `elif callback == ...` tree in the `assign_delegates` function, which
#   creates the C++ function object wrapper for a Python method and calls
#   `Delegator::setDelegate`
#
# ## Implementation for specific delegated functions
#
# Beyond the C++ implementation required in a class derived from `Delegator`, each
# delegated function needs only to have an entry in the corresponding Python class's
# `delegatable_methods` class variable. This variable is a mapping where the keys are
# the base names (without the `before_` / `replace_` / `after_` prefixes) of the
# delegate functions, using Python naming conventions, and the values are tuples of the
# corresponding C++ member function names and the matching signature of the C++ member
# function. These signatures should match the ones checked in the `assign_delegates`
# method.

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

# Wrapper for functions of type void(double*)
cdef void callback_v_dp(PyFuncInfo& funcInfo, size_array1 sizes, double* arg) noexcept:
    cdef double[:] view = <double[:sizes[0]]>arg if sizes[0] else None

    try:
        (<object>funcInfo.func())(view)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double, double*)
cdef void callback_v_d_dp(PyFuncInfo& funcInfo, size_array1 sizes, double arg1,
                          double* arg2) noexcept:
    cdef double[:] view = <double[:sizes[0]]>arg2 if sizes[0] else None

    try:
        (<object>funcInfo.func())(arg1, view)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double*, double*, double*)
cdef void callback_v_dp_dp_dp(PyFuncInfo& funcInfo,
        size_array3 sizes, double* arg1, double* arg2, double* arg3) noexcept:

    cdef double[:] view1 = <double[:sizes[0]]>arg1 if sizes[0] else None
    cdef double[:] view2 = <double[:sizes[1]]>arg2 if sizes[1] else None
    cdef double[:] view3 = <double[:sizes[2]]>arg3 if sizes[2] else None
    try:
        (<object>funcInfo.func())(view1, view2, view3)
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

# Wrapper for functions of type void(double, double*, double*)
cdef void callback_v_d_dp_dp(PyFuncInfo& funcInfo, size_array2 sizes, double arg1,
                             double* arg2, double* arg3) noexcept:
    cdef double[:] view1 = <double[:sizes[0]]>arg2 if sizes[0] else None
    cdef double[:] view2 = <double[:sizes[1]]>arg3 if sizes[1] else None

    try:
        (<object>funcInfo.func())(arg1, view1, view2)
    except BaseException as e:
        exc_type, exc_value = _sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

cdef int assign_delegates(obj, CxxDelegator* delegator) except -1:
    """
    Use methods defined in the Python class ``obj`` as delegates for the C++
    object ``delegator``. This function should be called in the ``__init__``
    method of classes where the wrapped C++ type is derived from the C++
    ``Delegator`` class.

    Methods that can be delegated are described by the ``delegatable_methods``
    dict of ``obj``.

    * The keys are the base names of the Python delegate methods. For methods where
      delegation is _optional_, the name is prefixed with ``before_``, ``after_``, or
      ``replace_`` in a specific implementation of a delegated class. For example, for
      the base name ``eval``, the delegate class can define one of these methods:
      ``before_eval``, ``after_eval``, or ``replace_eval``. For methods where delegation
      is _required_, no prefix is used.

    * The values are tuples of two or three elements, where the first element is the
      name of the corresponding C++ method, and the second element indicates the
      signature of the delegate function, such as ``void(double*)``. The third element,
      if present, indicates that the delegate is required, how it is executed with
      respect to the base class method (that is, ``before``, ``after``, or ``replace``).
    """
    delegator.setDelegatorName(stringify(obj.__class__.__name__))

    # Find all delegate methods, and make sure there aren't multiple
    # conflicting implementations
    cdef string cxx_name
    cdef string cxx_when
    obj._delegates = []
    for name, options in obj.delegatable_methods.items():
        if len(options) == 3:
            # Delegate with pre-selected mode, without using prefix on method name
            when = options[2]
            method = getattr(obj, name)
        else:
            when = None

        replace = 'replace_{}'.format(name)
        if hasattr(obj, replace):
            when = 'replace'
            method = getattr(obj, replace)

        before = 'before_{}'.format(name)
        if hasattr(obj, before):
            if when is not None:
                raise CanteraError(
                    "Only one delegate supported for '{}'".format(name))
            when = 'before'
            method = getattr(obj, before)

        after = 'after_{}'.format(name)
        if hasattr(obj, after):
            if when is not None:
                raise CanteraError(
                    "Only one delegate supported for '{}'".format(name))
            when = 'after'
            method = getattr(obj, after)

        if when is None:
            continue

        cxx_name = stringify(options[0])
        callback = options[1].replace(' ', '')

        # Make sure that the number of arguments needed by the C++ function
        # corresponds to the number of arguments accepted by the Python delegate
        if callback.endswith("()"):
            callback_args = 0
        else:
            callback_args = callback.count(",") + 1

        signature = _inspect.signature(method)
        params = signature.parameters.values()
        min_args = len([p for p in params if p.default is p.empty])
        max_args = len([p for p in params if p.kind is not p.KEYWORD_ONLY])

        if not (min_args <= callback_args <= max_args):
            raise ValueError(f"Function with signature {name}{signature}\n"
                "does not have the right number of arguments to be used as a delegate "
                f"for a function with the signature\n{callback}")

        cxx_when = stringify(when)
        if callback == 'void()':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v), cxx_when)
        elif callback == 'void(double)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_d), cxx_when)
        elif callback == 'void(AnyMap&)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_AMr), cxx_when)
        elif callback == 'void(AnyMap&,UnitStack&)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_cAMr_cUSr), cxx_when)
        elif callback == 'void(string,void*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_csr_vp), cxx_when)
        elif callback == 'void(double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_dp), cxx_when)
        elif callback == 'void(bool)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_b), cxx_when)
        elif callback == 'void(double,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_d_dp), cxx_when)
        elif callback == 'void(double*,double*,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_dp_dp_dp), cxx_when)
        elif callback == 'double(void*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_d_vp), cxx_when)
        elif callback == 'string(size_t)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_s_sz), cxx_when)
        elif callback == 'size_t(string)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_sz_csr), cxx_when)
        elif callback == 'void(double,double*,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(<PyObject*>method, callback_v_d_dp_dp), cxx_when)
        else:
            raise ValueError("Don't know how to set delegates for functions "
                f"with signature '{callback}'")

        # A Python object needs to hold references to the bound methods to prevent them
        # from being deleted, while still being eventually reachable by the garbage
        # collector
        obj._delegates.append(method)

    return 0


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

def extension(*, name, data=None):
    """
    A decorator for declaring Cantera extensions that should be registered with
    the corresponding factory classes to create objects with the specified *name*.

    This decorator can be used in combination with an ``extensions`` section in a YAML
    input file to trigger registration of extensions marked with this decorator,
    For example, consider an input file containing top level ``extensions`` and
    ``reactions`` sections such as:

    .. code:: yaml

        extensions:
        - type: python
          name: my_cool_module

        ...  # phases and species sections

        reactions:
        - equation: O + H2 <=> H + OH  # Reaction 3
          type: cool-rate
          A: 3.87e+04
          b: 2.7
          Ea: 6260.0

    and a Python module ``my_cool_module.py``::

        import cantera as ct

        class CoolRateData(ct.ExtensibleRateData):
            def update(self, soln):
                ...

        @ct.extension(name="cool-rate", data=CoolRateData)
        class CoolRate(ct.ExtensibleRate):
            def set_parameters(self, params, units):
                ...
            def eval(self, data):
                ...

    Loading this input file from any Cantera user interface would cause Cantera to load
    the ``my_cool_module.py`` module and register the ``CoolRate`` and ``CoolRateData``
    classes to handle reactions whose ``type`` in the YAML file is set to ``cool-rate``.

    .. versionadded:: 3.0
    """
    def decorator(cls):
        cdef shared_ptr[CxxExtensionManager] mgr = (
            CxxExtensionManagerFactory.build(stringify("python")))

        if issubclass(cls, ExtensibleRate):
            cls._reaction_rate_type = name
            # Registering immediately supports the case where the main
            # application is Python
            mgr.get().registerRateBuilder(
                stringify(cls.__module__), stringify(cls.__name__), stringify(name))

            # Deferred registration supports the case where the main application
            # is not Python
            _rate_delegators.append((cls.__module__, cls.__name__, name))

            # Register the ReactionData delegator
            if not issubclass(data, ExtensibleRateData):
                raise ValueError("'data' must inherit from 'ExtensibleRateData'")
            mgr.get().registerRateDataBuilder(
                stringify(data.__module__), stringify(data.__name__), stringify(name))
            _rate_data_delegators.append((data.__module__, data.__name__, name))
        else:
            raise TypeError(f"{cls} is not extensible")
        return cls

    return decorator
