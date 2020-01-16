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
