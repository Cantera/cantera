# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language=c++

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map as stdmap
from libcpp.cast cimport dynamic_cast
from libcpp.pair cimport pair
from libcpp cimport bool as cbool
from libcpp.functional cimport function
from libcpp.memory cimport shared_ptr, weak_ptr, dynamic_pointer_cast
from cpython cimport bool as pybool
from cpython.ref cimport PyObject
from cython.operator cimport dereference as deref, preincrement as inc

ctypedef stdmap[string,double] Composition

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector

# See https://github.com/cython/cython/pull/6539
# TODO: Replace once Cython >= 3.1.0 is required
cdef extern from "<span>" namespace "std" nogil:
    # Only Extent = std::dynamic_extent is supported until Cython can also
    # support integer templates. See https://github.com/cython/cython/pull/426
    cdef cppclass span[T]:
        ctypedef T value_type
        ctypedef size_t size_type
        ctypedef ptrdiff_t difference_type

        size_t extent

        cppclass iterator:
            iterator() except +
            iterator(iterator&) except +
            T& operator*()
            iterator operator++()
            iterator operator--()
            iterator operator++(int)
            iterator operator--(int)
            iterator operator+(size_type)
            iterator operator-(size_type)
            difference_type operator-(iterator)
            difference_type operator-(const_iterator)
            bint operator==(iterator)
            bint operator==(const_iterator)
            bint operator!=(iterator)
            bint operator!=(const_iterator)
            bint operator<(iterator)
            bint operator<(const_iterator)
            bint operator>(iterator)
            bint operator>(const_iterator)
            bint operator<=(iterator)
            bint operator<=(const_iterator)
            bint operator>=(iterator)
            bint operator>=(const_iterator)

        cppclass reverse_iterator:
            reverse_iterator() except +
            reverse_iterator(reverse_iterator&) except +
            T& operator*()
            reverse_iterator operator++()
            reverse_iterator operator--()
            reverse_iterator operator++(int)
            reverse_iterator operator--(int)
            reverse_iterator operator+(size_type)
            reverse_iterator operator-(size_type)
            difference_type operator-(iterator)
            difference_type operator-(const_iterator)
            bint operator==(reverse_iterator)
            bint operator==(const_reverse_iterator)
            bint operator!=(reverse_iterator)
            bint operator!=(const_reverse_iterator)
            bint operator<(reverse_iterator)
            bint operator<(const_reverse_iterator)
            bint operator>(reverse_iterator)
            bint operator>(const_reverse_iterator)
            bint operator<=(reverse_iterator)
            bint operator<=(const_reverse_iterator)
            bint operator>=(reverse_iterator)
            bint operator>=(const_reverse_iterator)

        span()
        span(T*, size_type) except +  # span[It](It, size_type)
        span(T*, T*) except +  # span[It, End](It, End)
        span(vector&)  # span[U, N](array[T, N]& arr)
        span(span&)

        T& operator[](ssize_t)

        T& back()
        iterator begin()
        T* data()
        bint empty()
        iterator end()
        span[T] first(size_type)
        T& front()
        span[T] last(size_type)
        reverse_iterator rbegin()
        reverse_iterator rend()
        size_type size()
        span[T] subspan(size_type)
        span[T] subspan(size_type, size_type)

    cdef size_t dynamic_extent


cdef extern from "cantera/cython/funcWrapper.h":
    cdef int translate_exception()

    cdef cppclass CxxAnyMap "Cantera::AnyMap"
    cdef cppclass CxxAnyValue "Cantera::AnyValue"
    cdef cppclass CxxUnits "Cantera::Units"
    cdef cppclass CxxUnitSystem "Cantera::UnitSystem"
