# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from libcpp.unordered_map cimport unordered_map

from .ctcxx cimport *

cdef extern from "cantera/base/AnyMap.h" namespace "Cantera":
    cdef cppclass CxxAnyValue "Cantera::AnyValue"

    cdef cppclass CxxAnyMap "Cantera::AnyMap":
        cppclass Iterator:
            pair[string, CxxAnyValue]& operator*()
            Iterator& operator++()
            cbool operator!=(Iterator&)

        cppclass OrderedIterator:
            pair[string, CxxAnyValue]& operator*()
            OrderedIterator& operator++()
            cbool operator!=(OrderedIterator&)

        cppclass OrderedProxy:
            OrderedIterator begin()
            OrderedIterator end()

        CxxAnyMap()
        Iterator begin()
        Iterator end()
        OrderedProxy ordered() except +translate_exception
        CxxAnyValue& operator[](string) except +translate_exception
        cbool empty()
        cbool hasKey(string)
        void clear()
        void update(CxxAnyMap& other, cbool)
        string keys_str()
        void applyUnits()

    cdef cppclass CxxAnyValue "Cantera::AnyValue":
        CxxAnyValue()
        CxxAnyValue& operator=(string) except +translate_exception
        CxxAnyValue& operator=(double) except +translate_exception
        CxxAnyValue& operator=(cbool) except +translate_exception
        CxxAnyValue& operator=(int) except +translate_exception
        CxxAnyValue& operator=(long) except +translate_exception
        CxxAnyValue& operator=(CxxAnyMap) except +translate_exception
        CxxAnyValue& operator=[T](vector[T]) except +translate_exception
        unordered_map[string, CxxAnyMap*] asMap(string) except +translate_exception
        CxxAnyMap& getMapWhere(string, string) except +translate_exception
        T& asType "as" [T]() except +translate_exception
        string type_str()
        cbool empty()
        cbool isType "is" [T]()
        cbool isScalar()

    CxxAnyMap AnyMapFromYamlFile "Cantera::AnyMap::fromYamlFile" (string) except +translate_exception
    CxxAnyMap AnyMapFromYamlString "Cantera::AnyMap::fromYamlString" (string) except +translate_exception

cdef extern from "cantera/base/stringUtils.h" namespace "Cantera":
    cdef Composition parseCompString(string) except +translate_exception

cdef extern from "cantera/base/global.h" namespace "Cantera":
    cdef void CxxAddDirectory "Cantera::addDirectory" (string)
    cdef string CxxGetDataDirectories "Cantera::getDataDirectories" (string)
    cdef void CxxAppdelete "Cantera::appdelete" ()
    cdef void Cxx_make_deprecation_warnings_fatal "Cantera::make_deprecation_warnings_fatal" ()
    cdef void Cxx_suppress_deprecation_warnings "Cantera::suppress_deprecation_warnings" ()
    cdef void Cxx_suppress_thermo_warnings "Cantera::suppress_thermo_warnings" (cbool)
    cdef void Cxx_use_legacy_rate_constants "Cantera::use_legacy_rate_constants" (cbool)
    cdef string CxxGitCommit "Cantera::gitCommit" ()
    cdef cbool CxxDebugModeEnabled "Cantera::debugModeEnabled" ()

cdef extern from "cantera/cython/utils_utils.h":
    cdef string get_cantera_version()
    cdef int get_sundials_version()
    cdef cppclass CxxPythonLogger "PythonLogger":
        pass

    cdef void CxxSetLogger "setLogger" (CxxPythonLogger*)

cdef string stringify(x) except *
cdef pystr(string x)

cdef comp_map_to_dict(Composition m)
cdef Composition comp_map(X) except *

cdef CxxAnyMap dict_to_anymap(data) except *
cdef anymap_to_dict(CxxAnyMap& m)
