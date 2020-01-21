# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
from numbers import Real as _real

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

    A `Func1` object representing a tabulated function is defined by sample
    points and corresponding function values. Inputs are specified either by
    two iterable objects containing sample point location and function values,
    or a single array that concatenates those inputs in two columns. Between
    sample points, values are evaluated using a linear interpolation, whereas
    outside the sample interval, the value at the closest end point is returned.
    Examples for tabulated `Func1` object are::

        >>> t1 = Func1([[0, 2], [1, 1], [2, 0]])
        >>> [t1(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 1.5, 0.5, 0.0, 0.0]

        >>> t2 = Func1([0, 1, 2], [2, 1, 0])
        >>> [t2(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 1.5, 0.5, 0.0, 0.0]

        >>> t3 = Func1(np.array([0, 1, 2]), (2, 1, 0))
        >>> [t3(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 1.5, 0.5, 0.0, 0.0]

    Note that all methods which accept `Func1` objects will also accept the
    callable object and create the wrapper on their own, so it is generally
    unnecessary to explicitly create a `Func1` object.
    """
    def __cinit__(self, *args, **kwargs):
        self.exception = None
        self.callable = None

    def __init__(self, *args):
        if len(args) == 1:
            c = args[0]
            if hasattr(c, '__call__'):
                # callback function
                self.__set_callback(c)
            else:
                arr = np.array(c)
                try:
                    if arr.ndim == 0:
                        # handle constants or unsized numpy arrays
                        k = float(c)
                        self.__set_callback(lambda t: k)
                    elif arr.size == 1:
                        # handle lists, tuples or numpy arrays with a single element
                        k = float(c[0])
                        self.__set_callback(lambda t: k)
                    elif arr.ndim == 2:
                        # tabulated function (single argument)
                        if arr.shape[1] == 2:
                            time = arr[:, 0]
                            fval = arr[:, 1]
                            self.__set_tables(time, fval)
                        else:
                            raise ValueError(
                                "Invalid dimensions: specification of "
                                "tabulated function with a single array "
                                "requires two columns")
                    else:
                        raise TypeError

                except TypeError:
                    raise TypeError(
                        "Func1 must be constructed from a callable object, "
                        "a number or a numeric array with two columns")

        elif len(args) == 2:
            # tabulated function (two arguments mimic C++ interface)
            time, fval = args
            self.__set_tables(time, fval)

        else:
            raise ValueError("Invalid number of arguments")


    cpdef void __set_callback(self, c) except *:
        self.callable = c
        self.func = new CxxFunc1(func_callback, <void*>self)

    cpdef void __set_tables(self, time, fval) except *:
        cdef vector[double] tvec, fvec
        for t in time:
            tvec.push_back(t)
        for f in fval:
            fvec.push_back(f)
        self.func = <CxxFunc1*>(new CxxTabulated1(tvec, fvec))

    def __dealloc__(self):
        del self.func

    def __call__(self, t):
        return self.func.eval(t)

    def __reduce__(self):
        raise NotImplementedError('Func1 object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Func1 object is not copyable')
