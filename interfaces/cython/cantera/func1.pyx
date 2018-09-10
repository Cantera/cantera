# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

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

    For simplicity, constant functions can be defined by passing the constant
    value directly::

        >>> f4 = Func1(2.5)
        >>> f4(0.1)
        2.5

    Note that all methods which accept `Func1` objects will also accept the
    callable object and create the wrapper on their own, so it is generally
    unnecessary to explicitly create a `Func1` object.
    """
    def __cinit__(self, c):
        self.exception = None
        if hasattr(c, '__call__'):
            self.callable = c
        else:
            try:
                # calling float() converts numpy arrays of size 1 to scalars
                k = float(c)
            except TypeError:
                if hasattr(c, '__len__') and len(c) == 1:
                    # Handle lists or tuples with a single element
                    k = float(c[0])
                else:
                    raise TypeError('Func1 must be constructed from a number or'
                                    ' a callable object')
            self.callable = lambda t: k

        self.func = new CxxFunc1(func_callback, <void*>self)

    def __dealloc__(self):
        del self.func

    def __call__(self, t):
        return self.func.eval(t)

    def __reduce__(self):
        raise NotImplementedError('Func1 object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Func1 object is not copyable')
