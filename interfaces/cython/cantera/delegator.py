# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3
# pyright: reportMissingImports=false, reportAttributeAccessIssue=false
# pyright: reportUndefinedVariable=false, reportUnboundVariable=false
# pyright: reportInvalidTypeArguments=false, reportAssignmentType=false
# pyright: reportIndexIssue=false, reportInvalidTypeForm=false

import cython
import inspect as _inspect
from typing import TYPE_CHECKING

from cython.cimports.cantera._utils import stringify
from ._utils import CanteraError
from cython.cimports.cantera.reaction import ExtensibleRate, ExtensibleRateData
from cython.cimports.cantera._delegate_callbacks import (
    callback_v, callback_v_d, callback_v_b, callback_v_AMr, callback_v_cAMr_cUSr,
    callback_v_csr_vp, callback_v_dp, callback_v_d_dp, callback_v_dp_dp_dp,
    callback_v_vETr, callback_d_vp, callback_s_sz, callback_sz_csr, callback_v_d_dp_dp)

if TYPE_CHECKING:
    # ``ExtensibleRate`` / ``ExtensibleRateData`` are cimported above for the runtime
    # ``issubclass`` checks; the type checkers cannot follow the ``cython.cimports``
    # path, so re-import them here under aliases (used only in the string annotations
    # below). A plain top-level Python import would create a cycle: ``delegator`` is
    # initialized before ``reaction`` in ``_cantera``.
    from collections.abc import Callable as _Callable
    from .reaction import (
        ExtensibleRate as _ExtensibleRate,
        ExtensibleRateData as _ExtensibleRateData,
    )

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
#
# The native-Cython machinery that these functions depend on (the `callback_*` wrapper
# functions and the `ct_*` functions called from C++) lives in the `_delegate_callbacks`
# module, because it relies on C++ reference parameters that cannot be expressed in
# Cython's pure-Python syntax.


@cython.cfunc
@cython.exceptval(-1, check=False)
def assign_delegates(obj, delegator: cython.pointer(CxxDelegator)) -> cython.int:
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
    cxx_name: string
    cxx_when: string
    obj._delegates = []
    valid_names = set()
    for name, options in obj.delegatable_methods.items():
        valid_names.add(name)
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
                pyOverride(cython.cast(cython.pointer(PyObject), method), callback_v),
                cxx_when)
        elif callback == 'void(double)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method), callback_v_d),
                cxx_when)
        elif callback == 'void(AnyMap&)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method), callback_v_AMr),
                cxx_when)
        elif callback == 'void(AnyMap&,UnitStack&)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method),
                           callback_v_cAMr_cUSr),
                cxx_when)
        elif callback == 'void(string,void*)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method),
                           callback_v_csr_vp),
                cxx_when)
        elif callback == 'void(double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method), callback_v_dp),
                cxx_when)
        elif callback == 'void(bool)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method), callback_v_b),
                cxx_when)
        elif callback == 'void(double,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method),
                           callback_v_d_dp),
                cxx_when)
        elif callback == 'void(double*,double*,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method),
                           callback_v_dp_dp_dp),
                cxx_when)
        elif callback == 'void(SparseTriplets&)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method),
                           callback_v_vETr),
                cxx_when)
        elif callback == 'double(void*)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method), callback_d_vp),
                cxx_when)
        elif callback == 'string(size_t)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method), callback_s_sz),
                cxx_when)
        elif callback == 'size_t(string)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method),
                           callback_sz_csr),
                cxx_when)
        elif callback == 'void(double,double*,double*)':
            delegator.setDelegate(cxx_name,
                pyOverride(cython.cast(cython.pointer(PyObject), method),
                           callback_v_d_dp_dp),
                cxx_when)
        else:
            raise ValueError("Don't know how to set delegates for functions "
                f"with signature '{callback}'")

        # A Python object needs to hold references to the bound methods to prevent them
        # from being deleted, while still being eventually reachable by the garbage
        # collector
        obj._delegates.append(method)

    # Check for methods on the base object that have no matching delegate
    for name in dir(obj):
        if (name.startswith("before_") or name.startswith("after_") or
            name.startswith("replace_")):
            base_name = name.split("_", 1)[1]
            if base_name not in valid_names:
                raise ValueError(f"'{base_name}' is not a delegatable method")

    return 0


def extension(
    *, name: str, data: type[_ExtensibleRateData] | None = None
) -> _Callable[[type[_ExtensibleRate]], type[_ExtensibleRate]]:
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
        from ._delegate_callbacks import _rate_delegators, _rate_data_delegators
        mgr: shared_ptr[CxxExtensionManager] = (
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
