# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

import sys
import os
from cpython.ref cimport PyObject

cdef int _pythonMajorVersion = sys.version_info[0]

cdef CxxPythonLogger* _logger = new CxxPythonLogger()
CxxSetLogger(_logger)

cdef string stringify(x) except *:
    """ Converts Python strings to std::string. """
    # This method works with both Python 2.x and 3.x.
    if isinstance(x, bytes):
        return string(<bytes>x)
    else:
        tmp = bytes(x.encode())
        return string(tmp)

cdef pystr(string x):
    cdef bytes s = x.c_str()
    if _pythonMajorVersion == 2:
        # Python 2.x
        return s
    else:
        # Python 3.x
        return s.decode()

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
    Cxx_make_deprecation_warnings_fatal()

def suppress_thermo_warnings(pybool suppress=True):
    Cxx_suppress_thermo_warnings(suppress)

cdef Composition comp_map(X) except *:
    if isinstance(X, (str, unicode, bytes)):
        return parseCompString(stringify(X))

    # assume X is dict-like
    cdef Composition m
    for species,value in X.items():
        m[stringify(species)] = value
    return m

cdef comp_map_to_dict(Composition m):
    return {pystr(species):value for species,value in m.items()}

class CanteraError(RuntimeError):
    pass

cdef public PyObject* pyCanteraError = <PyObject*>CanteraError
