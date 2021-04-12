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
    if v.isScalar():
        if v.isType[string]():
            return pystr(v.asType[string]())
        elif v.isType[double]():
            return v.asType[double]()
        elif v.isType[long]():
            return v.asType[long]()
        elif v.isType[cbool]():
            return v.asType[cbool]()
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
    return {pystr(item.first): anyvalue_to_python(item.first, item.second)
            for item in m.ordered()}
