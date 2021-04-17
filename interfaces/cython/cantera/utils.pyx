# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
import os
import warnings
from cpython.ref cimport PyObject

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
    CxxAddDirectory(stringify(directory))

def get_data_directories():
    """ Get a list of the directories Cantera searches for data files. """
    return pystr(CxxGetDataDirectories(stringify(os.pathsep))).split(os.pathsep)

__sundials_version__ = '.'.join(str(get_sundials_version()))

__version__ = pystr(get_cantera_version())

__git_commit__ = pystr(CxxGitCommit())

def appdelete():
    """ Delete all global Cantera C++ objects """
    CxxAppdelete()

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
    if v.isEmpty():
        # It is not possible to determine the associated type; return empty list
        return []
    if v.isScalar():
        if v.isType[string]():
            return pystr(v.asType[string]())
        elif v.isType[double]():
            return v.asType[double]()
        elif v.isType[long]():
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
    m.applyUnits()
    if m.isEmpty():
        return {}
    return {pystr(item.first): anyvalue_to_python(item.first, item.second)
            for item in m.ordered()}


cdef CxxAnyValue python_to_anyvalue(v, name=None) except *:
    # convert Python values to C++ AnyValue
    cdef CxxAnyMap a
    cdef CxxAnyValue b

    if name is None:
        msg = "Unable to convert value to AnyValue:\n"
    else:
        msg = "Unable to convert value with key '{}' to AnyValue:\n".format(name)

    if v is None:
        return b
    elif isinstance(v, set):
        raise NotImplementedError("{}cannot process Python set.".format(msg))
    elif isinstance(v, dict):
        a = dict_to_anymap(v)
        b = a
        return b

    def is_scalar(v):
        return not hasattr(v, '__len__') or isinstance(v, str)

    all_same = True
    multi_dict = False
    if not is_scalar(v):
        # ensure that data are homogeneous
        # (np.array converts inhomogeneous sequences to string arrays)
        if not len(v):
            # np.array converts empty arrays to float
            pass
        elif isinstance(v, np.ndarray):
            pass
        elif is_scalar(v[0]):
            all_same = all([type(val) == type(v[0]) for val in v])
        elif isinstance(v[0], dict):
            multi_dict = all([isinstance(val, dict) for val in v])
        else:
            all_same = []
            sizes = []
            for row in range(len(v)):
                sizes.append(len(v[row]))
                all_same.append(all([type(val) == type(v[0][0]) for val in v[row]]))
            if not all(all_same):
                raise NotImplementedError(
                        "{}cannot process nested sequences with inhomogeneous data "
                        "types.".format(msg)
                    )
            elif np.unique(sizes).size != 1:
                raise NotImplementedError(
                        "{}cannot process arrays represented by ragged nested "
                        "sequences.".format(msg)
                    )
            all_same = all(all_same)

    cdef vector[CxxAnyValue] vv_any
    if not all_same:
        # inhomogeneous sequence
        for val in v:
            vv_any.push_back(python_to_anyvalue(val, name))
        b = vv_any
        return b

    cdef vector[CxxAnyMap] vv_map
    if multi_dict:
        for val in v:
            vv_map.push_back(dict_to_anymap(val))
        b = vv_map
        return b

    # data are homogeneous: convert to np.array for convenience
    vv = np.array(v)
    ndim = vv.ndim

    cdef double v_double
    cdef vector[double] vv_double
    cdef vector[vector[double]] vvv_double
    if vv.dtype == float:
        if vv.size == 0:
            # empty list
            return b
        if ndim == 0:
            v_double = v
            b = v_double
            return b
        if ndim == 1:
            for val in vv:
                vv_double.push_back(val)
            b = vv_double
            return b
        if ndim == 2:
            vvv_double.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_double[row].push_back(val)
            b = vvv_double
            return b
        raise NotImplementedError(
                "{}cannot process float array with {} dimensions".format(msg, ndim)
            )

    cdef long v_int
    cdef vector[long] vv_int
    cdef vector[vector[long]] vvv_int
    if vv.dtype == int:
        if ndim == 0:
            v_int = v
            b = v_int
            return b
        if ndim == 1:
            for val in vv:
                vv_int.push_back(val)
            b = vv_int
            return b
        if ndim == 2:
            vvv_int.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_int[row].push_back(val)
            b = vvv_int
            return b
        raise NotImplementedError(
                "{}cannot process integer array with {} dimensions".format(msg, ndim)
            )

    cdef cbool v_bool
    cdef vector[cbool] vv_bool
    cdef vector[vector[cbool]] vvv_bool
    if vv.dtype == bool:
        if ndim == 0:
            v_bool = v
            b = v_bool
            return b
        if ndim == 1:
            for val in vv:
                vv_bool.push_back(val)
            b = vv_bool
            return b
        if ndim == 2:
            vvv_bool.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_bool[row].push_back(val)
            b = vvv_bool
            return b
        raise NotImplementedError(
                "{}cannot process boolean array with {} dimensions".format(msg, ndim)
            )

    cdef string v_string
    cdef vector[string] vv_string
    cdef vector[vector[string]] vvv_string
    if isinstance(v, str) or (ndim and isinstance(vv.ravel()[0], str)):
        if ndim == 0:
            v_string = stringify(v)
            b = v_string
            return b
        if ndim == 1:
            for val in vv:
                vv_string.push_back(stringify(val))
            b = vv_string
            return b
        if ndim == 2:
            vvv_string.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_string[row].push_back(stringify(val))
            b = vvv_string
            return b
        raise NotImplementedError(
                "{}cannot process string array with {} dimensions".format(msg, ndim)
            )

    raise NotImplementedError(
            "{}unknown conversion for variable with value\n".format(msg,v)
        )


cdef CxxAnyMap dict_to_anymap(dict dd) except *:
    # convert Python dictionary to C++ AnyMap
    cdef CxxAnyMap mm
    cdef string kk
    for key, val in dd.items():
        kk = stringify(key)
        mm[kk] = python_to_anyvalue(val, key)
    return mm


def _py_to_any_to_py(dd):
    # @internal  used for testing purposes only
    cdef string name = stringify("test")
    cdef CxxAnyValue vv = python_to_anyvalue(dd)
    return anyvalue_to_python(name, vv)
