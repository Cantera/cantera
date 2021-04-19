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


cdef CxxAnyValue python_to_anyvalue(item, name=None) except *:
    # convert Python values to C++ AnyValue
    cdef CxxAnyMap a
    cdef CxxAnyValue b

    if name is None:
        msg = "Unable to convert item of type '{}' to AnyValue:\n".format(type(item))
    else:
        msg = ("Unable to convert item of type '{}' with key '{}' to AnyValue:\n"
               "".format(type(item), name))

    if item is None:
        return b
    elif isinstance(item, set):
        raise NotImplementedError("{}cannot process Python set.".format(msg))
    elif isinstance(item, dict):
        return <CxxAnyValue>dict_to_anymap(item)

    def is_scalar(item):
        return not hasattr(item, '__len__') or isinstance(item, str)

    itype = set()
    if not (is_scalar(item) or isinstance(item, np.ndarray)):
        itype = set(type(i) for i in item)
        # ensure that data are homogeneous
        # (np.array converts inhomogeneous sequences to string arrays)
        if not len(item):
            # np.array converts empty arrays to float
            pass
        elif is_scalar(item[0]):
            # the content may be inhomogeneous or 1-D homogeneous, but it will
            # not be a 2-D array and thus does not have to be checked
            pass
        elif itype == {dict}:
            # list of dictionaries
            pass
        else:
            itype = set()
            sizes = []
            for row in item:
                sizes.append(len(row))
                itype = itype.union(set(type(val) for val in row))
            if len(itype) != 1:
                raise NotImplementedError(
                    "{}cannot process nested sequences with inhomogeneous data "
                    "types.".format(msg))
            elif np.unique(sizes).size != 1:
                raise NotImplementedError(
                    "{}cannot process arrays represented by ragged nested "
                    "sequences.".format(msg))

    cdef vector[CxxAnyValue] vv_any
    if len(itype) > 1:
        # inhomogeneous sequence
        for val in item:
            vv_any.push_back(python_to_anyvalue(val, name))
        b = vv_any
        return b

    cdef vector[CxxAnyMap] vv_map
    if itype == {dict}:
        for val in item:
            vv_map.push_back(dict_to_anymap(val))
        b = vv_map
        return b

    # data are homogeneous: convert to np.array for convenience
    vv = np.array(item)
    ndim = vv.ndim

    cdef vector[double] vv_double
    cdef vector[vector[double]] vvv_double
    if vv.dtype == float:
        if vv.size == 0:
            # empty list
            pass
        elif ndim == 0:
            b = <double>item
        elif ndim == 1:
            for val in vv:
                vv_double.push_back(<double>val)
            b = vv_double
        elif ndim == 2:
            vvv_double.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_double[row].push_back(<double>val)
            b = vvv_double
        else:
            raise NotImplementedError(
                "{}cannot process float array with {} dimensions".format(msg, ndim))
        return b

    cdef vector[long] vv_int
    cdef vector[vector[long]] vvv_int
    if vv.dtype == int:
        if ndim == 0:
            b = <long>item
        elif ndim == 1:
            for val in vv:
                vv_int.push_back(<long>val)
            b = vv_int
        elif ndim == 2:
            vvv_int.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_int[row].push_back(<long>val)
            b = vvv_int
        else:
            raise NotImplementedError(
                "{}cannot process integer array with {} dimensions".format(msg, ndim))
        return b

    cdef vector[cbool] vv_bool
    cdef vector[vector[cbool]] vvv_bool
    if vv.dtype == bool:
        if ndim == 0:
            b = <cbool>item
        elif ndim == 1:
            for val in vv:
                vv_bool.push_back(<cbool>val)
            b = vv_bool
        elif ndim == 2:
            vvv_bool.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_bool[row].push_back(<cbool>val)
            b = vvv_bool
        else:
            raise NotImplementedError(
                "{}cannot process boolean array with {} dimensions".format(msg, ndim))
        return b

    cdef vector[string] vv_string
    cdef vector[vector[string]] vvv_string
    if isinstance(item, str) or (ndim and isinstance(vv.ravel()[0], str)):
        if ndim == 0:
            b = stringify(item)
        elif ndim == 1:
            for val in vv:
                vv_string.push_back(stringify(val))
            b = vv_string
        elif ndim == 2:
            vvv_string.resize(vv.shape[0])
            for row in range(vv.shape[0]):
                for val in vv[row, :]:
                    vvv_string[row].push_back(stringify(val))
            b = vvv_string
        else:
            raise NotImplementedError(
                "{}cannot process string array with {} dimensions".format(msg, ndim))
        return b

    raise NotImplementedError(
        "{}unknown conversion for variable with value\n".format(msg, item))


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
