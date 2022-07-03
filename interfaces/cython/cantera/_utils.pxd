#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from cython.operator cimport dereference as deref, preincrement as inc

cdef string stringify(x) except *
cdef pystr(string x)

cdef comp_map_to_dict(Composition m)
cdef Composition comp_map(X) except *

cdef CxxAnyMap dict_to_anymap(data) except *
cdef anymap_to_dict(CxxAnyMap& m)
