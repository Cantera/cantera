# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import os
import warnings
from cpython.ref cimport PyObject
from libcpp.utility cimport move
import numbers as _numbers
import importlib as _importlib
from collections import namedtuple as _namedtuple
import numpy as np
from .units cimport Units


_scipy_sparse = None
def _import_scipy_sparse():
    # defer scipy import
    global _scipy_sparse
    if _scipy_sparse is not None:
        return

    try:
        _importlib.metadata.version('scipy')
    except _importlib.metadata.PackageNotFoundError:
        raise ImportError('Method requires a working scipy installation.')
    else:
        from scipy import sparse as _scipy_sparse


cdef unique_ptr[CxxPythonLogger] _logger = make_unique[CxxPythonLogger]()
CxxSetLogger(move(_logger))

cdef string stringify(x) except *:
    """ Converts Python strings to std::string. """
    if x is None:
        return stringify("")
    if isinstance(x, bytes):
        return string(<bytes>x)
    else:
        tmp = bytes(x.encode())
        return string(tmp)

cdef pystr(string x):
    return x.decode()

def add_directory(directory):
    """ Add a directory to search for Cantera data files. """
    CxxAddDirectory(stringify(str(directory)))

def get_data_directories():
    """ Get a list of the directories Cantera searches for data files. """
    return pystr(CxxGetDataDirectories(stringify(os.pathsep))).split(os.pathsep)

__sundials_version__ = pystr(get_sundials_version())

__version__ = pystr(CxxVersion())

if __version__ != pystr(get_cantera_version_py()):
    raise ImportError("Mismatch between Cantera Python module version "
        f"({pystr(get_cantera_version_py())}) and Cantera shared library "
        f"version ({__version__})")

__git_commit__ = pystr(CxxGitCommit())

if __git_commit__ != pystr(get_cantera_git_commit_py()):
    raise ImportError("Mismatch between Cantera Python module Git commit "
        f"({pystr(get_cantera_git_commit_py())}) and Cantera shared library "
        f"git commit ({__git_commit__})")

_USE_SPARSE = False

def debug_mode_enabled():
    return CxxDebugModeEnabled()

def print_stack_trace_on_segfault():
    """
    Enable printing a stack trace if a segfault occurs. Not recommended for general
    use as it is possible for this to deadlock.

    .. versionadded:: 3.0
    """
    CxxPrintStackTraceOnSegfault()

def appdelete():
    """ Delete all global Cantera C++ objects """
    CxxAppdelete()

def use_sparse(sparse=True):
    """
    Enable sparse output using `scipy.sparse`. Sparse output requires a working
    *SciPy* installation. Use pip or conda to install ``scipy`` to enable this method.
    """
    global _USE_SPARSE
    if sparse:
        try:
            _import_scipy_sparse()
        except ImportError:
            raise
    _USE_SPARSE = sparse

def make_deprecation_warnings_fatal():
    warnings.filterwarnings('error', category=DeprecationWarning,
                            module='cantera')  # for warnings in Python code
    warnings.filterwarnings('error', category=DeprecationWarning,
                            message='.*Cantera.*')  # for warnings in Cython code
    Cxx_make_deprecation_warnings_fatal()

def suppress_deprecation_warnings():
    warnings.filterwarnings('ignore', category=DeprecationWarning,
                            module='cantera')  # for warnings in Python code
    warnings.filterwarnings('ignore', category=DeprecationWarning,
                            message='.*Cantera.*')  # for warnings in Cython code
    Cxx_suppress_deprecation_warnings()

def suppress_thermo_warnings(pybool suppress=True):
    Cxx_suppress_thermo_warnings(suppress)

def use_legacy_rate_constants(pybool legacy):
    """
    Set definition used for rate constant calculation.

    If set to `False` (default value), rate constants of three-body reactions are
    consistent with conventional definitions (for example Eq. 9.75 in
    :cite:t:`kee2003`). If set to `True`, output for rate constants of three-body
    reactions is multiplied by third-body concentrations, consistent with Cantera's
    behavior prior to version 3.0.
    """
    Cxx_use_legacy_rate_constants(legacy)

def hdf_support():
    """
    Returns list of libraries that include HDF support:
    - 'native': if Cantera was compiled with C++ HighFive HDF5 support.

    .. versionadded:: 3.0
    """
    out = []
    if CxxUsesHDF5():
        out.append("native")
    return set(out)

cdef Composition comp_map(X) except *:
    if isinstance(X, (str, bytes)):
        return parseCompString(stringify(X))

    # assume X is dict-like
    cdef Composition m
    for species,value in (<object>X).items():
        m[stringify(species)] = value
    return m

cdef comp_map_to_dict(Composition m):
    return {pystr(species):value for species,value in (<object>m).items()}

class CanteraError(RuntimeError):
    @staticmethod
    def set_stack_trace_depth(depth):
        """
        Set the number of stack frames to include when a `CanteraError` is displayed. By
        default, or if the depth is set to 0, no stack information will be shown.
        """
        CxxCanteraError.setStackTraceDepth(depth)

_DimensionalValue = _namedtuple('_DimensionalValue',
                                ('value', 'units', 'activation_energy'),
                                defaults=[False])

cdef public PyObject* pyCanteraError = <PyObject*>CanteraError


cdef class AnyMap(dict):
    """
    A key-value store representing objects defined in Cantera's YAML input format.

    Extends the capabilities of a normal `dict` object by providing functions for
    converting values between different unit systems. See :ref:`sec-yaml-units` for
    details on how units are handled in YAML input files.
    """
    def __cinit__(self, *args, **kwawrgs):
        self.unitsystem = UnitSystem()

    cdef _set_CxxUnitSystem(self, shared_ptr[CxxUnitSystem] units):
        self.unitsystem._set_unitSystem(units)

    def default_units(self):
        return self.unitsystem.defaults()

    @property
    def units(self):
        """Get the `UnitSystem` applicable to this `AnyMap`."""
        return self.unitsystem

    def convert(self, str key, dest):
        """
        Convert the value corresponding to the specified *key* to the units defined by
        *dest*. *dest* may be a string or a `Units` object.
        """
        return self.unitsystem.convert_to(self[key], dest)

    def convert_activation_energy(self, key, dest):
        """
        Convert the value corresponding to the specified *key* to the units defined by
        *dest*. *dest* may be a string or a `Units` object defining units that are
        interpretable as an activation energy.
        """
        return self.unitsystem.convert_activation_energy_to(self[key], dest)

    def convert_rate_coeff(self, str key, dest):
        """
        Convert the value corresponding to the specified *key* to the units defined by
        *dest*, with special handling for `UnitStack` input and potentially-undefined
        rate coefficient units.
        """
        return self.unitsystem.convert_rate_coeff_to(self[key], dest)

    def set_quantity(self, str key, value, src):
        """
        Set the element *key* of this map to the specified value, converting from the
        units defined by *src* to the correct unit system for this map when serializing
        to YAML.
        """
        self[key] = _DimensionalValue(value, src)

    def set_activation_energy(self, str key, value, src):
        """
        Set the element *key* of this map to the specified value, converting from the
        activation energy units defined by *src* to the correct unit system for this map
        when serializing to YAML.
        """
        self[key] = _DimensionalValue(value, src, True)


cdef anyvalue_to_python(string name, CxxAnyValue& v):
    cdef CxxAnyMap a
    cdef CxxAnyValue b
    if v.empty():
        # It is not possible to determine the associated type; return None
        return None
    if v.isScalar():
        if v.isType[string]():
            return pystr(v.asType[string]())
        elif v.isType[double]():
            return v.asType[double]()
        elif v.isType[long]():
            # 'long' is equivalent to 'long int'
            # Cython requires the former in this context
            return v.asType[long]()
        elif v.isType[cbool]():
            return v.asType[cbool]()
        else:
            raise TypeError("Unable to convert value with key '{}' "
                            "from AnyValue of held type '{}'".format(
                                pystr(name), v.type_str()))
    elif v.isType[CxxAnyMap]():
        return anymap_to_py(v.asType[CxxAnyMap]())
    elif v.isType[vector[CxxAnyMap]]():
        return [anymap_to_py(a) for a in v.asType[vector[CxxAnyMap]]()]
    elif v.isType[vector[double]]():
        return v.asType[vector[double]]()
    elif v.isType[vector[string]]():
        return [pystr(s) for s in v.asType[vector[string]]()]
    elif v.isType[vector[long]]():
        return v.asType[vector[long]]()
    elif v.isType[vector[cbool]]():
        return v.asType[vector[cbool]]()
    elif v.isType[vector[CxxAnyValue]]():
        return [anyvalue_to_python(name, b)
                for b in v.asType[vector[CxxAnyValue]]()]
    elif v.isType[vector[vector[double]]]():
        return v.asType[vector[vector[double]]]()
    elif v.isType[vector[vector[string]]]():
        return [[pystr(s) for s in row]
                for row in v.asType[vector[vector[string]]]()]
    elif v.isType[vector[vector[long]]]():
        return v.asType[vector[vector[long]]]()
    elif v.isType[vector[vector[cbool]]]():
        return v.asType[vector[vector[cbool]]]()
    else:
        raise TypeError("Unable to convert value with key '{}' "
                        "from AnyValue of held type '{}'".format(
                            pystr(name), v.type_str()))


cdef anymap_to_py(CxxAnyMap& m):
    cdef pair[string,CxxAnyValue] item
    m.applyUnits()
    cdef AnyMap out = AnyMap()
    out._set_CxxUnitSystem(m.unitsShared())
    for item in m.ordered():
        out[pystr(item.first)] = anyvalue_to_python(item.first, item.second)
    return out


cdef void setQuantity(CxxAnyMap& m, str k, v: _DimensionalValue) except *:
    cdef CxxAnyValue testval = python_to_anyvalue(v.value)
    cdef CxxAnyValue target
    if isinstance(v.units, str):
        if testval.isScalar():
            target.setQuantity(testval.asType[double](), stringify(v.units),
                               <cbool?>v.activation_energy)
        else:
            target.setQuantity(testval.asVector[double](), stringify(v.units))
    elif isinstance(v.units, Units):
        target.setQuantity(testval.asType[double](), (<Units>v.units).units)
    else:
        raise TypeError(f'Expected a string or Units object. Got {type(v.units)}')
    m[stringify(k)] = target


cdef CxxAnyMap py_to_anymap(data, cbool hyphenize=False) except *:
    cdef CxxAnyMap m
    if hyphenize:
        # replace "_" by "-": while Python dictionaries typically use "_" in key names,
        # the YAML convention uses "-" in field names
        def _hyphenize(data):
            if isinstance(data, dict):
                return {k.replace("_", "-"): _hyphenize(v) for k, v in data.items()}
            return data

        data = _hyphenize(data)

    for k, v in data.items():
        if isinstance(v, _DimensionalValue):
            setQuantity(m, k, v)
        else:
            m[stringify(k)] = python_to_anyvalue(v, k)
    return m

cdef get_types(item):
    """ Helper function used by python_to_anyvalue """
    if not len(item):
        # Empty list, so no specific type can be inferred
        return None, 1

    elif isinstance(item, np.ndarray):
        if isinstance(item.flat[0], _numbers.Integral):
            itype = _numbers.Integral
        elif isinstance(item.flat[0], _numbers.Real):
            itype = _numbers.Real
        elif isinstance(item.flat[0], (str, np.str_)):
            itype = str
        else:
            itype = item.dtype.type
        return itype, item.ndim

    else:
        itype = set()
        for i in item:
            if isinstance(i, dict):
                # Treat all classes derived from dict as equivalent, and
                # convertible to AnyMap
                itype.add(dict)
            elif isinstance(i, (list, tuple)):
                # Treat all classes derived from list or tuple as equivalent,
                # and convertible to std::vector
                itype.add(list)
            elif type(i) == bool:
                itype.add(bool)  # otherwise bools will get counted as integers
            elif isinstance(i, _numbers.Integral):
                itype.add(_numbers.Integral)
            elif isinstance(i, _numbers.Real):
                itype.add(_numbers.Real)
            else:
                itype.add(type(i))

        if itype == {_numbers.Real, _numbers.Integral}:
            itype = {_numbers.Real}

        if itype == {list}:
            inner_types = set()
            ndim_inner = set()
            for j in item:
                type_j, ndim_j = get_types(j)
                inner_types.add(type_j)
                ndim_inner.add(ndim_j)
            if len(inner_types) == 1 and len(ndim_inner) == 1:
                # Inner types and dimensions match up to second nesting level.
                # Checks for higher levels are skipped as the items are converted
                # to vector<AnyValue> rather than vector<vector<type>>.
                return inner_types.pop(), ndim_inner.pop() + 1
            else:
                return None, 1
        elif len(itype) == 1:
            return itype.pop(), 1
        else:
            return None, 1

cdef CxxAnyValue python_to_anyvalue(item, name=None) except *:
    cdef CxxAnyValue v
    if name is not None:
        v.setKey(stringify(name))
    if isinstance(item, dict):
        v = py_to_anymap(item)
    elif isinstance(item, (list, tuple, set, np.ndarray)):
        itype, ndim = get_types(item)
        if ndim == 1:
            if itype == str:
                v = list_string_to_anyvalue(item)
            elif itype == bool:
                v = list_bool_to_anyvalue(item)
            elif itype == _numbers.Integral:
                v = list_int_to_anyvalue(item)
            elif itype == _numbers.Real:
                v = list_double_to_anyvalue(item)
            elif itype == dict:
                v = list_dict_to_anyvalue(item)
            else:
                v = list_to_anyvalue(item)
        elif ndim == 2:
            if itype == str:
                v = list2_string_to_anyvalue(item)
            elif itype == bool:
                v = list2_bool_to_anyvalue(item)
            elif itype == _numbers.Integral:
                v = list2_int_to_anyvalue(item)
            elif itype == _numbers.Real:
                v = list2_double_to_anyvalue(item)
            else:
                v = list_to_anyvalue(item)
        else:
            v = list_to_anyvalue(item)
    elif isinstance(item, (str, np.str_, np.bytes_)):
        v = stringify(item)
    elif isinstance(item, (bool, np.bool_)):
        v = <cbool>(item)
    elif isinstance(item, (int, np.int32, np.int64)):
        v = <long int>(item)
    elif isinstance(item, (float, np.float32, np.float64)):
        v = <double>(item)
    elif item is None:
        pass  # None corresponds to "empty" AnyValue
    elif name is not None:
        raise CanteraError("Unable to convert item of type {!r}"
            " with key {!r} to AnyValue".format(type(item), name))
    else:
        raise CanteraError("Unable to convert item of type {!r}"
            " to AnyValue".format(type(item)))
    return v

# Helper functions for converting specific types to AnyValue

cdef vector[CxxAnyValue] list_to_anyvalue(data) except *:
    cdef vector[CxxAnyValue] v
    v.resize(len(data))
    cdef size_t i
    for i, item in enumerate(data):
        v[i] = python_to_anyvalue(item)
    return v

cdef vector[double] list_double_to_anyvalue(data):
    cdef vector[double] v
    v.resize(len(data))
    cdef size_t i
    for i, item in enumerate(data):
        v[i] = <double>item
    return v

cdef vector[long] list_int_to_anyvalue(data):
    cdef vector[long] v
    v.resize(len(data))
    cdef size_t i
    for i, item in enumerate(data):
        v[i] = <long int>item
    return v

cdef vector[cbool] list_bool_to_anyvalue(data):
    cdef vector[cbool] v
    v.resize(len(data))
    cdef size_t i
    for i, item in enumerate(data):
        v[i] = <cbool>item
    return v

cdef vector[string] list_string_to_anyvalue(data):
    cdef vector[string] v
    v.resize(len(data))
    cdef size_t i
    for i, item in enumerate(data):
        v[i] = stringify(item)
    return v

cdef vector[CxxAnyMap] list_dict_to_anyvalue(data) except *:
    cdef vector[CxxAnyMap] v
    v.resize(len(data))
    cdef size_t i
    for i, item in enumerate(data):
        v[i] = py_to_anymap(item)
    return v

cdef vector[vector[double]] list2_double_to_anyvalue(data):
    cdef vector[vector[double]] v
    v.resize(len(data))
    cdef size_t i, j
    for i, item in enumerate(data):
        v[i].resize(len(item)) # allows for ragged nested lists
        for j, jtem in enumerate(item):
            v[i][j] = <double>jtem
    return v

cdef vector[vector[long]] list2_int_to_anyvalue(data):
    cdef vector[vector[long]] v
    v.resize(len(data))
    cdef size_t i, j
    for i, item in enumerate(data):
        v[i].resize(len(item)) # allows for ragged nested lists
        for j, jtem in enumerate(item):
            v[i][j] = <long int>jtem
    return v

cdef vector[vector[cbool]] list2_bool_to_anyvalue(data):
    cdef vector[vector[cbool]] v
    v.resize(len(data))
    cdef size_t i, j
    for i, item in enumerate(data):
        v[i].resize(len(item)) # allows for ragged nested lists
        for j, jtem in enumerate(item):
            v[i][j] = <cbool>jtem
    return v

cdef vector[vector[string]] list2_string_to_anyvalue(data):
    cdef vector[vector[string]] v
    v.resize(len(data))
    cdef size_t i, j
    for i, item in enumerate(data):
        v[i].resize(len(item)) # allows for ragged nested lists
        for j, jtem in enumerate(item):
            v[i][j] = stringify(jtem)
    return v

def _py_to_any_to_py(dd):
    # used for internal testing purposes only
    cdef string name = stringify("test")
    cdef CxxAnyValue vv = python_to_anyvalue(dd)
    return anyvalue_to_python(name, vv), pystr(vv.type_str())

def _py_to_anymap_to_py(pp):
    # used for internal testing purposes only
    cdef CxxAnyMap m = py_to_anymap(pp)
    return anymap_to_py(m)
