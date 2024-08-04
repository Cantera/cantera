# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import sys
cimport numpy as np
import numpy as np
import warnings

from ._utils cimport *


cdef double func_callback(double t, void* obj, void** err) except? 0.0:
    """
    This function is called from C/C++ to evaluate a `Func1` object ``obj``,
    returning the value of the function at ``t``. If an exception occurs while
    evaluating the function, the Python exception info is saved in the
    two-element array ``err``.
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
    This class is used as a wrapper for a function of one variable,
    :math:`y = f(t)`, that is defined in Python and can be called by the
    Cantera C++ core. `Func1` objects are constructed from callable Python
    objects, for example functions or classes which implement the `__call__` method::

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

    Note that `Func1` objects also allow for direct access to functor objects
    implemented in C++ based on associated type specifiers::

        >>> f5 = Func1("exp", 3.)  # C++ 'Exp1' functor
        >>> f5.cxx_type
        'Cantera::Exp1'
        >>> f5.write()
        '\\exp(3x)'
        >>> f5(2.)
        403.4287934927351
        >>> f6 = Func1("Arrhenius", [9630.0, 2.0, 2012.878])  # C++ 'Arrhenius1' functor
        >>> f6(1500)
        5662665826.195515

    For implemented C++ functor types, see the Cantera C++ :ct:`Func1` documentation.

    `Func1` objects support operator overloading which facilitates the construction of
    compound functions, where some standard simplifications are implemented::

        >>> f7 = 2 * f5 + 3
        >>> f7.write()
        '2\\exp(3x) + 3'
        >>> f7(2.)
        809.8575869854702
        >>> f8 = f5 * f5
        >>> f8.write()
        '\\exp(6x)'
        >>> f8(2.)
        162754.79141900392

    .. versionchanged:: 3.1

        Implementations for operator overloading and direct support for C++ functor
        construction are added.
    """
    def __cinit__(self, *args, **kwargs):
        self.exception = None
        self.callable = None

    def __init__(self, c, *args, init=True):
        if init is False:
            # used by 'create' classmethod
            return
        if hasattr(c, '__call__'):
            # callback function
            self._set_callback(c)
            return

        cdef shared_ptr[CxxFunc1] cxx_func
        if isinstance(c, str):
            cxx_func = Func1._make_cxx_func1(stringify(c), args)
            self._func = cxx_func
            self.func = cxx_func.get()
            return

        try:
            arr = np.array(c)
            if arr.ndim == 0:
                # handle constants or unsized numpy arrays
                k = float(c)
            elif arr.size == 1:
                # handle lists, tuples or numpy arrays with a single element
                k = float(arr.flat[0])
            else:
                raise TypeError
            cxx_func = Func1._make_cxx_func1(stringify("constant"), (k,))
            self._func = cxx_func
            self.func = cxx_func.get()

        except TypeError:
            raise TypeError(
                "'Func1' objects must be constructed from a number or "
                "a callable object") from None

    cpdef void _set_callback(self, c) except *:
        self.callable = c
        self._func.reset(new CxxFunc1Py(func_callback, <void*>self))
        self.func = self._func.get()

    @property
    def type(self):
        """
        Return the type of the underlying C++ functor object.

        .. versionadded:: 3.0
        """
        return pystr(self.func.type())

    @staticmethod
    cdef shared_ptr[CxxFunc1] _make_cxx_func1(string cxx_string, tuple args):
        """Create C++ functor from type specifier and arguments."""
        cdef shared_ptr[CxxFunc1] func
        cdef Func1 f0
        cdef Func1 f1
        cdef vector[double] arr
        func1_type = pystr(CxxCheckFunc1(cxx_string))
        if func1_type == "undefined":
            functor_type = pystr(cxx_string)
            raise NotImplementedError(f"Functor '{functor_type}' is not implemented.")
        if len(args) == 0 and func1_type == "standard":
            # basic functor with no parameter
            func = CxxNewFunc1(cxx_string, 1.)
        elif len(args) == 1 and func1_type == "standard":
            if hasattr(args[0], "__len__"):
                # advanced functor with array and no parameter
                for v in args[0]:
                    arr.push_back(v)
                func = CxxNewFunc1(cxx_string, arr)
            else:
                # basic functor with scalar parameter
                func = CxxNewFunc1(cxx_string, float(args[0]))
        elif len(args) == 2:
            if func1_type == "compound":
                # compounding functor
                if isinstance(args[0], Func1) and isinstance(args[1], Func1):
                    f0 = args[0]
                    f1 = args[1]
                elif isinstance(args[0], Func1) and not hasattr(args[1], "__len__"):
                    f0 = args[0]
                    f1 = Func1(float(args[1]))
                elif isinstance(args[1], Func1) and not hasattr(args[0], "__len__"):
                    f0 = Func1(float(args[0]))
                    f1 = args[1]
                else:
                    raise ValueError(f"Invalid arguments for compounding functor.")
                func = CxxNewFunc1(cxx_string, f0._func, f1._func)
            elif func1_type == "modified":
                # modifying functor
                if isinstance(args[0], Func1) and isinstance(args[1], (float, int)):
                    f0 = args[0]
                else:
                    raise ValueError(f"Invalid arguments for modifying functor.")
                func = CxxNewFunc1(cxx_string, f0._func, float(args[1]))
            elif func1_type == "standard":
                # tabulating functor
                if hasattr(args[0], "__len__") and hasattr(args[1], "__len__"):
                    for v in args[0]:
                        arr.push_back(v)
                    for v in args[1]:
                        arr.push_back(v)
                else:
                    raise ValueError(f"Invalid arguments for tabulating functor.")
                func = CxxNewFunc1(cxx_string, arr)
            else:
                raise ValueError("Invalid arguments")
        else:
            raise ValueError("Invalid arguments")
        return func

    @staticmethod
    cdef Func1 _make_func1(shared_ptr[CxxFunc1] func):
        """Create Python Func1 from C++ functor."""
        cdef Func1 out = Func1(None, init=False)
        out._func = func
        out.func = out._func.get()
        return out

    @staticmethod
    def cxx_functor(functor_type, *args):
        """
        Retrieve a C++ `Func1` functor (advanced feature).

        For implemented functor types, see the Cantera C++ ``Func1`` documentation.

        .. versionadded:: 3.0

        .. deprecated:: 3.1

            To be removed after Cantera 3.1; replaced by alternative constructor.
        """
        warnings.warn(
            "To be removed after Cantera 3.1; use alternative constructor instead.",
            DeprecationWarning)
        cdef shared_ptr[CxxFunc1] func
        func = Func1._make_cxx_func1(stringify(functor_type), args)
        return Func1._make_func1(func)

    def __add__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewSumFunction(self._func, f1._func))

    def __radd__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewSumFunction(f1._func, self._func))

    def __sub__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewDiffFunction(self._func, f1._func))

    def __rsub__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewDiffFunction(f1._func, self._func))

    def __mul__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewProdFunction(self._func, f1._func))

    def __rmul__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewProdFunction(f1._func, self._func))

    def __truediv__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewRatioFunction(self._func, f1._func))

    def __rtruediv__(self, other):
        if not isinstance(other, Func1):
            other = Func1(other)
        cdef Func1 f1 = other
        return Func1._make_func1(CxxNewRatioFunction(f1._func, self._func))

    @property
    def cxx_type(self):
        """
        Return the type of the underlying C++ functor object.

        .. versionadded:: 3.1
        """
        return pystr(self.func.typeName())

    def write(self, name="x"):
        """
        Write a :math:`LaTeX` expression representing a functor.

        :param name:
            Name of the variable to be used.

        .. versionadded:: 3.0
        """
        return pystr(self.func.write(stringify(name)))

    def __call__(self, t):
        return self.func.eval(t)

    def __reduce__(self):
        msg = "'{}' objects are not picklable".format(type(self).__name__)
        raise NotImplementedError(msg)

    def __copy__(self):
        msg = "'{}' objects are not copyable".format(type(self).__name__)
        raise NotImplementedError(msg)


cdef class Tabulated1(Func1):
    """
    A `Tabulated1` object representing a tabulated function is defined by
    sample points and corresponding function values. Inputs are specified by
    two iterable objects containing sample point location and function values.
    Between sample points, values are evaluated based on the optional argument
    ``method``; options are ``'linear'`` (linear interpolation, default) or
    ``'previous'`` (nearest previous value). Outside the sample interval, the
    value at the closest end point is returned.

    Examples for `Tabulated1` objects are::

        >>> t1 = Tabulated1([0, 1, 2], [2, 1, 0])
        >>> [t1(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 1.5, 0.5, 0.0, 0.0]

        >>> t2 = Tabulated1(np.array([0, 1, 2]), np.array([2, 1, 0]))
        >>> [t2(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 1.5, 0.5, 0.0, 0.0]

    The optional ``method`` keyword argument changes the type of interpolation
    from the ``'linear'`` default to ``'previous'``::

        >>> t3 = Tabulated1([0, 1, 2], [2, 1, 0], method='previous')
        >>> [t3(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]]
        [2.0, 2.0, 2.0, 1.0, 0.0, 0.0]

    .. versionadded:: 3.0
    """

    def __init__(self, time, fval, method='linear'):
        cdef vector[double] arr
        for v in time:
            arr.push_back(v)
        for v in fval:
            arr.push_back(v)
        cdef string cxx_string = stringify(f"tabulated-{method}")
        self._func = CxxNewFunc1(cxx_string, arr)
        self.func = self._func.get()
