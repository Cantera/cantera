# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys

cdef double func_callback(double t, void* obj, void** err):
    """
    This function is called from C/C++ to evaluate a `Func1` object *obj*,
    returning the value of the function at *t*. If an exception occurs while
    evaluating the function, the Python exception info is saved in the
    two-element array *err*.
    """
    try:
        return (<Func1>obj).callable(t)
    except BaseException as e:
        exc_type, exc_value = sys.exc_info()[:2]

        # Stash the exception info to prevent it from being garbage collected
        (<Func1>obj).exception = exc_type, exc_value
        err[0] = <void*>exc_type
        err[1] = <void*>exc_value
        return 0.0


cdef class Func1:
    """
    This class is used as a wrapper for a function of one variable, i.e.
    :math:`y = f(t)`, that is defined in Python and can be called by the
    Cantera C++ core. `Func1` objects are constructed from callable Python
    objects, e.g. functions or classes which implement the `__call__` method::

        >>> f1 = Func1(math.sin)
        >>> f1(math.pi/4)
        0.7071067811865475

        >>> f2 = Func1(lambda t: t**2 + 1)
        >>> f2(3)
        10

        >>> class Multiplier:
        ...     def __init__(self, factor):
        ...         self.factor = factor
        ...     def __call__(self, t):
        ...         return self.factor * t
        >>> f3 = Func1(Multiplier(5))
        >>> f3(6)
        30.0

    For simplicity, constant functions can be defined by passing a constant
    value directly::

        >>> f4 = Func1(2.5)
        >>> f4(0.1)
        2.5

    Note that all methods which accept `Func1` objects will also accept the
    callable object and create the wrapper on their own, so it is often not
    necessary to explicitly create a `Func1` object.
    """
    def __cinit__(self, *args, **kwargs):
        self.exception = None
        self.callable = None

    def __init__(self, c):
        if hasattr(c, '__call__'):
            # callback function
            self._set_callback(c)
        else:
            arr = np.array(c)
            try:
                if arr.ndim == 0:
                    # handle constants or unsized numpy arrays
                    k = float(c)
                    self._set_callback(lambda t: k)
                elif arr.size == 1:
                    # handle lists, tuples or numpy arrays with a single element
                    k = float(c[0])
                    self._set_callback(lambda t: k)
                else:
                    raise TypeError

            except TypeError:
                raise TypeError(
                    "'Func1' objects must be constructed from a number or "
                    "a callable object") from None

    cpdef void _set_callback(self, c) except *:
        self.callable = c
        self.func = new CxxFunc1(func_callback, <void*>self)

    def __dealloc__(self):
        del self.func

    def __call__(self, t):
        return self.func.eval(t)

    def __reduce__(self):
        msg = "'{}' objects are not picklable".format(type(self).__name__)
        raise NotImplementedError(msg)

    def __copy__(self):
        msg = "'{}' objects are not copyable".format(type(self).__name__)
        raise NotImplementedError(msg)


cdef class TabulatedFunction(Func1):
    """
    A `TabulatedFunction` object representing a tabulated function is defined by
    sample points and corresponding function values. Inputs are specified by
    two iterable objects containing sample point location and function values.
    Between sample points, values are evaluated based on the optional argument
    ``method``; options are ``'linear'`` (linear interpolation, default) or
    ``'previous'`` (nearest previous value). Outside the sample interval, the
    value at the closest end point is returned.

    Examples for `TabulatedFunction` objects are::

        >>> t1 = TabulatedFunction([0, 1, 2], [2, 1, 0])
        >>> [t1(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 1.5, 0.5, 0.0, 0.0]

        >>> t2 = TabulatedFunction(np.array([0, 1, 2]), np.array([2, 1, 0]))
        >>> [t2(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 1.5, 0.5, 0.0, 0.0]

    The optional ``method`` keyword argument changes the type of interpolation
    from the ``'linear'`` default to ``'previous'``::

        >>> t3 = TabulatedFunction([0, 1, 2], [2, 1, 0], method='previous')
        >>> [t3(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 2.0, 1.0, 0.0, 0.0]
    """

    def __init__(self, time, fval, method='linear'):
        self._set_tables(time, fval, stringify(method))

    cpdef void _set_tables(self, time, fval, string method) except *:
        tt = np.asarray(time, dtype=np.double)
        ff = np.asarray(fval, dtype=np.double)
        if tt.size != ff.size:
            raise ValueError("Sizes of arrays do not match "
                             "({} vs {})".format(tt.size, ff.size))
        elif tt.size == 0:
            raise ValueError("Arrays must not be empty.")
        cdef np.ndarray[np.double_t, ndim=1] tvec = tt
        cdef np.ndarray[np.double_t, ndim=1] fvec = ff
        self.func = <CxxFunc1*>(new CxxTabulated1(tt.size, &tvec[0], &fvec[0],
                                                  method))
