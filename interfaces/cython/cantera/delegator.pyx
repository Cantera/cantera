# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import inspect

# ## Implementation for each delegated function type
#
# Besides the C++ functions implemented in the `Delegator` class, each delegated
# function type requires several additional components in Cython:
# - A delcaration for the overload of `Delegator::setDelegate` in `_cantera.pxd`.
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
cdef void callback_v(PyFuncInfo& funcInfo):
    try:
        (<object>funcInfo.func())()
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double)
cdef void callback_v_d(PyFuncInfo& funcInfo, double arg):
    try:
        (<object>funcInfo.func())(arg)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(bool)
cdef void callback_v_b(PyFuncInfo& funcInfo, cbool arg):
    try:
        (<object>funcInfo.func())(arg)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double*)
cdef void callback_v_dp(PyFuncInfo& funcInfo, size_array1 sizes, double* arg):
    cdef double[:] view = <double[:sizes[0]]>arg if sizes[0] else None

    try:
        (<object>funcInfo.func())(view)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double, double*)
cdef void callback_v_d_dp(PyFuncInfo& funcInfo, size_array1 sizes, double arg1,
                                double* arg2):
    cdef double[:] view = <double[:sizes[0]]>arg2 if sizes[0] else None

    try:
        (<object>funcInfo.func())(arg1, view)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type void(double*, double*, double*)
cdef void callback_v_dp_dp_dp(PyFuncInfo& funcInfo,
        size_array3 sizes, double* arg1, double* arg2, double* arg3):

    cdef double[:] view1 = <double[:sizes[0]]>arg1 if sizes[0] else None
    cdef double[:] view2 = <double[:sizes[1]]>arg2 if sizes[1] else None
    cdef double[:] view3 = <double[:sizes[2]]>arg3 if sizes[2] else None
    try:
        (<object>funcInfo.func())(view1, view2, view3)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)

# Wrapper for functions of type string(size_t)
cdef int callback_s_sz(PyFuncInfo& funcInfo, string& out, size_t arg):
    try:
        ret = (<object>funcInfo.func())(arg)
        if ret is None:
            return 0
        else:
            (&out)[0] = stringify(ret)
            return 1
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1

# Wrapper for functions of type size_t(string&)
cdef int callback_sz_csr(PyFuncInfo& funcInfo, size_t& out, const string& arg):
    try:
        ret = (<object>funcInfo.func())(pystr(arg))
        if ret is None:
            return 0
        else:
            (&out)[0] = ret
            return 1
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
        funcInfo.setExceptionType(<PyObject*>exc_type)
        funcInfo.setExceptionValue(<PyObject*>exc_value)
    return -1

# Wrapper for functions of type void(double, double*, double*)
cdef void callback_v_d_dp_dp(PyFuncInfo& funcInfo, size_array2 sizes, double arg1,
                                double* arg2, double* arg3):
    cdef double[:] view1 = <double[:sizes[0]]>arg2 if sizes[0] else None
    cdef double[:] view2 = <double[:sizes[1]]>arg3 if sizes[1] else None

    try:
        (<object>funcInfo.func())(arg1, view1, view2)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]
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

    * The keys are the base names of the Python delegate methods,
      which may be prefixed with ``before_``, ``after_``, or ``replace_`` in
      a specific implementation of a delegated class. For example, for the base
      name ``eval``, the delegate class can define one of these methods:
      ``before_eval``, ``after_eval``, or ``replace_eval``.

    * The values are tuples of two elements, where the first element is the name
      of the corresponding C++ method, and the second element indicates the
      signature of the delegate function, such as ``void(double*)``.
    """
    # Find all delegate methods, and make sure there aren't multiple
    # conflicting implementations
    cdef string cxx_name
    cdef string cxx_when
    obj._delegates = []
    for name in obj.delegatable_methods:
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

        cxx_name = stringify(obj.delegatable_methods[name][0])
        callback = obj.delegatable_methods[name][1].replace(' ', '')

        # Make sure that the number of arguments needed by the C++ function
        # corresponds to the number of arguments accepted by the Python delegate
        if callback.endswith("()"):
            callback_args = 0
        else:
            callback_args = callback.count(",") + 1

        signature = inspect.signature(method)
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
