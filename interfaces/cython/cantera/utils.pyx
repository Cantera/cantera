# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import os
import warnings
from cpython.ref cimport PyObject
import numbers
import pkg_resources

# avoid explicit dependence of cantera on scipy
try:
    pkg_resources.get_distribution('scipy')
except pkg_resources.DistributionNotFound:
    _scipy_sparse = ImportError('Method requires a working scipy installation.')
else:
    from scipy import sparse as _scipy_sparse

cdef CxxPythonLogger* _logger = new CxxPythonLogger()
CxxSetLogger(_logger)

cdef string stringify(x) except *:
    """ Converts Python strings to std::string. """
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

__sundials_version__ = '.'.join(str(get_sundials_version()))

__version__ = pystr(get_cantera_version())

__git_commit__ = pystr(CxxGitCommit())

_USE_SPARSE = False

def debug_mode_enabled():
    return CxxDebugModeEnabled()

def appdelete():
    """ Delete all global Cantera C++ objects """
    CxxAppdelete()

def use_sparse(sparse=True):
    """
    Enable sparse output using `scipy.sparse`. Sparse output requires a working
    *SciPy* installation. Use pip or conda to install ``scipy`` to enable this method.
    """
    global _USE_SPARSE
    if sparse and isinstance(_scipy_sparse, ImportError):
        raise _scipy_sparse
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

    If set to 'False', rate constants of three-body reactions are consistent with
    conventional definitions. If set to 'True', output for rate constants of
    three-body reactions is multipied by third-body concentrations (legacy behavior).
    For the pre-compiled Cantera 2.6 distribution, the default value is set to 'True',
    which implies no change compared to previous behavior. For user-compiled Cantera,
    the default behavior can be changed by the SCons flag 'legacy_rate_constants'.

    .. deprecated:: 2.6

        Behavior to change after Cantera 2.6; for Cantera 2.6, rate constants of
        three-body reactions are multiplied with third-body concentrations
        (no change to legacy behavior). After Cantera 2.6, results will no longer
        include third-body concentrations and be consistent with conventional
        definitions (see Eq. 9.75 in Kee, Coltrin and Glarborg, 'Chemically
        Reacting Flow', Wiley Interscience, 2003).
    """
    Cxx_use_legacy_rate_constants(legacy)

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
    pass

cdef public PyObject* pyCanteraError = <PyObject*>CanteraError

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
        return anymap_to_dict(v.asType[CxxAnyMap]())
    elif v.isType[vector[CxxAnyMap]]():
        return [anymap_to_dict(a) for a in v.asType[vector[CxxAnyMap]]()]
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


cdef anymap_to_dict(CxxAnyMap& m):
    cdef pair[string,CxxAnyValue] item
    m.applyUnits()
    if m.empty():
        return {}
    return {pystr(item.first): anyvalue_to_python(item.first, item.second)
            for item in m.ordered()}

cdef CxxAnyMap dict_to_anymap(data) except *:
    cdef CxxAnyMap m
    for k, v in data.items():
        m[stringify(k)] = python_to_anyvalue(v, k)
    return m

cdef get_types(item):
    """ Helper function used by python_to_anyvalue """
    if not len(item):
        # Empty list, so no specific type can be inferred
        return None, 1

    elif isinstance(item, np.ndarray):
        if isinstance(item.flat[0], numbers.Integral):
            itype = numbers.Integral
        elif isinstance(item.flat[0], numbers.Real):
            itype = numbers.Real
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
            elif isinstance(i, numbers.Integral):
                itype.add(numbers.Integral)
            elif isinstance(i, numbers.Real):
                itype.add(numbers.Real)
            else:
                itype.add(type(i))

        if itype == {numbers.Real, numbers.Integral}:
            itype = {numbers.Real}

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
    if isinstance(item, dict):
        v = dict_to_anymap(item)
    elif isinstance(item, (list, tuple, set, np.ndarray)):
        itype, ndim = get_types(item)
        if ndim == 1:
            if itype == str:
                v = list_string_to_anyvalue(item)
            elif itype == bool:
                v = list_bool_to_anyvalue(item)
            elif itype == numbers.Integral:
                v = list_int_to_anyvalue(item)
            elif itype == numbers.Real:
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
            elif itype == numbers.Integral:
                v = list2_int_to_anyvalue(item)
            elif itype == numbers.Real:
                v = list2_double_to_anyvalue(item)
            else:
                v = list_to_anyvalue(item)
        else:
            v = list_to_anyvalue(item)
    elif isinstance(item, str):
        v = stringify(item)
    elif isinstance(item, bool):
        v = <cbool>(item)
    elif isinstance(item, int):
        v = <long int>(item)
    elif isinstance(item, float):
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
        v[i] = dict_to_anymap(item)
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
    # @internal  used for testing purposes only
    cdef string name = stringify("test")
    cdef CxxAnyValue vv = python_to_anyvalue(dd)
    return anyvalue_to_python(name, vv), pystr(vv.type_str())
