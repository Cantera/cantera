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
from cantera._utils cimport stringify
from cantera.ctcxx cimport shared_ptr
from cantera.solutionbase cimport CxxSolution, _assign_Solution
from cantera.reaction cimport ExtensibleRate, CxxReactionRate, ExtensibleRateData
from cantera.delegator cimport CxxDelegator, assign_delegates


cdef extern from "cantera/kinetics/ReactionRateDelegator.h" namespace "Cantera":
    cdef cppclass CxxReactionDataDelegator "Cantera::ReactionDataDelegator" (CxxDelegator):
        CxxReactionDataDelegator()

    cdef cppclass CxxReactionRateDelegator "Cantera::ReactionRateDelegator" (CxxDelegator, CxxReactionRate):
        CxxReactionRateDelegator()


cdef extern from "cantera/base/ExtensionManager.h" namespace "Cantera":
    cdef cppclass CxxExtensionManager "Cantera::ExtensionManager":
        void registerRateBuilder(string&, string&, string&)
        void registerRateDataBuilder(string&, string&, string&)

cdef extern from "cantera/base/ExtensionManagerFactory.h" namespace "Cantera":
    cdef cppclass CxxExtensionManagerFactory "Cantera::ExtensionManagerFactory":
        @staticmethod
        shared_ptr[CxxExtensionManager] build(string&)

cdef public char* ct_getExceptionString(object exType, object exValue, object exTraceback):
    import traceback
    result = str(exValue) + "\n\n"
    result += "".join(traceback.format_exception(exType, exValue, exTraceback))
    tmp = bytes(result.encode())
    cdef char* c_string = <char*> malloc((len(tmp) + 1) * sizeof(char))
    strcpy(c_string, tmp)
    return c_string


cdef public object ct_newPythonExtensibleRate(CxxReactionRateDelegator* delegator,
                                              const string& module_name,
                                              const string& class_name):

    mod = importlib.import_module(module_name.decode())
    cdef ExtensibleRate rate = getattr(mod, class_name.decode())(init=False)
    rate.set_cxx_object(delegator)
    return rate


cdef public object ct_newPythonExtensibleRateData(CxxReactionDataDelegator* delegator,
                                                  const string& module_name,
                                                  const string& class_name):

    mod = importlib.import_module(module_name.decode())
    cdef ExtensibleRateData data = getattr(mod, class_name.decode())()
    data.set_cxx_object(delegator)
    return data


cdef public ct_registerReactionDelegators():
    cdef shared_ptr[CxxExtensionManager] mgr = (
        CxxExtensionManagerFactory.build(stringify("python")))

    for module, cls, name in ct.delegator._rate_delegators:
        mgr.get().registerRateBuilder(stringify(module), stringify(cls), stringify(name))

    ct.delegator._rate_delegators.clear()

    for module, cls, name in ct.delegator._rate_data_delegators:
        mgr.get().registerRateDataBuilder(stringify(module), stringify(cls), stringify(name))

    ct.delegator._rate_data_delegators.clear()


cdef public object ct_wrapSolution(shared_ptr[CxxSolution] soln):
    pySoln = ct.Solution(init=False)
    _assign_Solution(pySoln, soln, False, weak=True)
    return pySoln
