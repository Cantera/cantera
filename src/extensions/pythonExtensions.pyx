# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from libcpp.string cimport string
from libc.stdlib cimport malloc
from libc.string cimport strcpy
from cpython.ref cimport Py_INCREF

import importlib
import inspect

import cantera as ct
from cantera.reaction cimport ExtensibleRate, CxxReactionRate
from cantera.delegator cimport CxxDelegator, assign_delegates


cdef extern from "cantera/kinetics/ReactionRateDelegator.h" namespace "Cantera":
    cdef cppclass CxxReactionRateDelegator "Cantera::ReactionRateDelegator" (CxxDelegator, CxxReactionRate):
        CxxReactionRateDelegator()


cdef public char* ct_getPythonExtensibleRateTypes(const string& module_name) except NULL:
    """
    Load the named module and find classes derived from ExtensibleRate.

    Returns a string where each line contains the class name and the corresponding
    rate name, separated by a space
    """
    mod = importlib.import_module(module_name.decode())
    names = "\n".join(
        f"{name}\t{cls._reaction_rate_type}"
        for name, cls in inspect.getmembers(mod)
        if inspect.isclass(cls) and issubclass(cls, ct.ExtensibleRate))
    tmp = bytes(names.encode())
    cdef char* c_string = <char*> malloc((len(tmp) + 1) * sizeof(char))
    strcpy(c_string, tmp)
    return c_string


cdef public object ct_newPythonExtensibleRate(CxxReactionRateDelegator* delegator,
                                              const string& module_name,
                                              const string& class_name):

    mod = importlib.import_module(module_name.decode())
    cdef ExtensibleRate rate = getattr(mod, class_name.decode())(init=False)
    rate.set_cxx_object(delegator)
    assign_delegates(rate, delegator)
    return rate
