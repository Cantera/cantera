# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from libcpp.string cimport string
from libc.stdlib cimport malloc
from libc.string cimport strcpy

import importlib
import inspect

import cantera as ct
from cantera.reaction cimport ExtensibleRate

cdef public char* ct_getPythonExtensibleRateTypes(const string& module_name):
    """
    Load the named module and find classes derived from ExtensibleRate.

    Returns a string where each line contains the class name and the corresponding
    rate name, separated by a space
    """
    mod = importlib.import_module(module_name.decode())
    names = "\n".join(
        f"{name} {cls._reaction_rate_type}"
        for name, cls in inspect.getmembers(mod)
        if inspect.isclass(cls) and issubclass(cls, ct.ExtensibleRate))
    tmp = bytes(names.encode())
    cdef char* c_string = <char*> malloc((len(tmp) + 1) * sizeof(char))
    strcpy(c_string, tmp)
    return c_string
