"""
The classes in this module are designed to allow constructing
user-defined functions of one variable in Python that can be used with the
Cantera C++ kernel. These classes are mostly shadow classes for
corresponding classes in the C++ kernel.
"""

from Cantera.num import array, asarray, ravel, shape, transpose
import _cantera
import types


class Func1:
    """
    Functors of one variable.

    A Functor is an object that behaves like a function. :class:`Func1`
    is the base class from which several functor classes derive. These
    classes are designed to allow specifying functions of time from Python
    that can be used by the C++ kernel.

    Functors can be added, multiplied, and divided to yield new functors.

    >>> f1 = Polynomial([1.0, 0.0, 3.0])  # 3*t*t + 1
    >>> f1(2.0)
    13
    >>> f2 = Polynomial([-1.0, 2.0])      # 2*t - 1
    >>> f2(2.0)
    5
    >>> f3 = f1/f2                        # (3*t*t + 1)/(2*t - 1)
    >>> f3(2.0)
    4.3333333
    """

    def __init__(self, typ, n, coeffs=[]):
        """
        The constructor is meant to be called from constructors of subclasses
        of Func1: :class:`Polynomial`, :class:`Gaussian`, :class:`Arrhenius`,
        :class:`Fourier`, :class:`Const`, :class:`PeriodicFunction`.
        """
        self.n = n
        self._own = 1
        self._func_id = 0
        self._typ = typ
        if _cantera.nummod == 'numpy':
            self.coeffs = array(coeffs, dtype=float, ndmin=1)
        else:
            self.coeffs = asarray(coeffs,'d')
        self._func_id = _cantera.func_new(typ, n, self.coeffs)

    def __del__(self):
        if self._func_id and self._own:
            _cantera.func_del(self._func_id)

    def __repr__(self):
        return self.write()

    def __call__(self, t):
        """Implements function syntax, so that F(t) is equivalent to
        F.value(t)."""
        if type(t) == types.NoneType:
            return self
        if type(t) == types.InstanceType:
            return CompositeFunction(self, t)
        else:
            return _cantera.func_value(self._func_id, t)

    def __add__(self, other):
        """Overloads operator '+'

        Returns a new function self(t) + other(t)"""

        # if 'other' is a number, then create a 'Const' functor for
        # it.
        if type(other) == types.FloatType:
            return SumFunction(self, Const(other))

        return SumFunction(self, other)

    def __radd__(self, other):
        """Overloads operator '+'

        Returns a new function other(t) + self(t)"""
        # if 'other' is a number, then create a 'Const' functor for
        # it.
        if type(other) == types.FloatType:
            return SumFunction(Const(other),self)
        return SumFunction(other, self)

    def __sub__(self, other):
        """Overloads operator '-'

        Returns a new function self(t) - other(t)"""

        # if 'other' is a number, then create a 'Const' functor for
        # it.
        if type(other) != types.InstanceType:
            return DiffFunction(self, Const(other))

        return DiffFunction(self, other)

    def __rsub__(self, other):
        """Overloads operator '-'

        Returns a new function other(t) - self(t)"""

        # if 'other' is a number, then create a 'Const' functor for
        # it.
        if type(other) != types.InstanceType:
            return DiffFunction(Const(other), self)
        return DiffFunction(other, self)

    def __mul__(self, other):
        """Overloads operator '*'

           Return a new function self(t)*other(t)"""
        if type(other) != types.InstanceType:
            return ProdFunction(self, Const(other))
        return ProdFunction(self, other)

    def __rmul__(self, other):
        """Overloads operator '*'

        Returns a new function other(t)*self(t)"""
        if type(other) != types.InstanceType:
            return ProdFunction(Const(other), self)
        return ProdFunction(other, self)

    def __div__(self, other):
        """Overloads operator '/'

        Returns a new function self(t)/other(t)"""
        if type(other) != types.InstanceType:
            return RatioFunction(self, Const(other))
        return RatioFunction(self, other)

    def __rdiv__(self, other):
        """Overloads operator '/'

        Returns a new function other(t)/self(t)"""
        if type(other) != types.InstanceType:
            return RatioFunction(Const(other), self)
        return RatioFunction(other, self)

    def func_id(self):
        """Internal. Return the integer index used internally to access the
        kernel-level object."""
        return self._func_id

    def write(self, arg = 'x', length = 1000):
        return _cantera.func_write(self._func_id, length, arg)


class Sin(Func1):
    def __init__(self,omega=1.0):
        Func1.__init__(self,100,1,omega)
class Cos(Func1):
    def __init__(self, omega=1.0):
        Func1.__init__(self,102,1,omega)
class Exp(Func1):
    def __init__(self,A=1.0):
        Func1.__init__(self,104,1,A)
class Pow(Func1):
    def __init__(self, n):
        Func1.__init__(self,106,1,n)

class Polynomial(Func1):
    r"""
    A polynomial.
    Instances of class 'Polynomial' evaluate

    .. math:: f(t) = \sum_{n = 0}^N a_n t^n .

    The coefficients are supplied as a list, beginning with :math:`a_N` and
    ending with :math:`a_0`.

    >>> p1 = Polynomial([1.0, -2.0, 3.0])   #    3t^2 - 2t + 1
    >>> p2 = Polynomial([6.0, 8.0])         #    8t + 6

    """
    def __init__(self, coeffs=[]):
        """
        coeffs - polynomial coefficients
        """
        Func1.__init__(self, 2, len(coeffs)-1, coeffs)



class Gaussian(Func1):
    r"""A Gaussian pulse. Instances of class 'Gaussian' evaluate

    .. math::  f(t) = A \exp[-(t - t_0) / \tau]

    where

    .. math:: \tau = \frac{\mbox{FWHM}}{2.0\sqrt{\ln(2.0)}}

    'FWHM' denotes the full width at half maximum.

    As an example, here is how to create a Gaussian pulse with peak amplitude
    10.0, centered at time 2.0, with full-width at half max = 0.2:

    >>> f = Gaussian(A = 10.0, t0 = 2.0, FWHM = 0.2)
    >>> f(2.0)
    10
    >>> f(1.9)
    5
    >>> f(2.1)
    5
    """
    def __init__(self, A, t0, FWHM):
        coeffs = array([A, t0, FWHM], 'd')
        Func1.__init__(self, 4, 0, coeffs)


class Fourier(Func1):
    r"""
    Fourier series. Instances of class 'Fourier' evaluate the Fourier series

    .. math::

        f(t) = \frac{a_0}{2} +
            \sum_{n=1}^N [a_n \cos(n\omega t) + b_n \sin(n \omega t)]

    where

    .. math::

        a_n = \frac{\omega}{\pi}
        \int_{-\pi/\omega}^{\pi/\omega} f(t) \cos(n \omega t) dt

        b_n = \frac{\omega}{\pi}
        \int_{-\pi/\omega}^{\pi/\omega} f(t) \sin(n \omega t) dt.

    The function :math:`f(t)` is periodic, with period :math:`T = 2\pi/\omega`.

    As an example, a function with Fourier components up to the second harmonic
    is constructed as follows:

    >>> coeffs = [(a0, b0), (a1, b1), (a2, b2)]
    >>> f = Fourier(omega, coeffs)

    Note that ``b0`` must be specified, but is not used. The value of ``b0``
    is arbitrary.
    """
    def __init__(self, omega, coefficients):
        """
        :param omega:
            fundamental frequency [radians/sec].
        :param coefficients:
            List of (a,b) pairs, beginning with n = 0.
        """
        cc = asarray(coefficients,'d')
        n, m = cc.shape
        if m <> 2:
            raise CanteraError('provide (a, b) for each term')
        cc[0,1] = omega
        Func1.__init__(self, 1, n-1, ravel(transpose(cc)))


class Arrhenius(Func1):
    r"""Sum of modified Arrhenius terms. Instances of class 'Arrhenius' evaluate

    .. math::  f(T) = \sum_{n=1}^N A_n T^{b_n}\exp(-E_n/T)

    Example:

    >>> f = Arrhenius([(a0, b0, e0), (a1, b1, e1)])
    """
    def __init__(self, coefficients):
        """
        :param coefficients:
            sequence of (*A*, *b*, *E*) triplets.
        """
        cc = asarray(coefficients,'d')
        n, m = cc.shape
        if m <> 3:
            raise CanteraError('Three Arrhenius parameters (A, b, E) required.')
        Func1.__init__(self, 3, n, ravel(cc))


class Const(Func1):
    """Constant function.
    Objects created by function Const act as functions that have a constant
    value. These are used internally whenever a statement like

    >>> f = Gausian(2.0, 1.0, 0.1) + 4.0

    is encountered. The addition operator of class Func1 is defined so that
    this is equivalent to

    >>> f = SumFunction(Gaussian(2.0, 1.0, 0.1), Const(4.0))

    Function Const returns instances of class Polynomial that have
    degree zero, with the constant term set to the desired value.
    """
    def __init__(self, value):
        Func1.__init__(self,110,1,value)
    #return Polynomial([value])


class PeriodicFunction(Func1):
    """Converts a function into a periodic function with period T."""
    def __init__(self, func, T):
        """
        :param func:
            initial non-periodic function
        :param T:
            period [s]
        """
        Func1.__init__(self, 50, func.func_id(), array([T],'d'))
        func._own = 0


# functions that combine two functions

class ComboFunc1(Func1):
    """
    Combines two functions.
    This class is the base class for functors that combine two
    other functors in a binary operation.
    """
    def __init__(self, typ, f1, f2):
        self._own = 1
        self._func_id = 0
        self._typ = typ
        if type(f1) == types.IntType:
            f1 = Const(f1)
        if type(f2) == types.IntType:
            f2 = Const(f2)
        self.f1 = f1
        self.f2 = f2
        self.f1._own = 0
        self.f2._own = 0
        self._func_id = _cantera.func_newcombo(typ, f1.func_id(), f2.func_id())


class SumFunction(ComboFunc1):
    """Sum of two functions.
    Instances of class SumFunction evaluate the sum of two supplied functors.
    It is not necessary to explicitly create an instance of SumFunction, since
    the addition operator of the base class is overloaded to return a SumFunction
    instance.

    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1 + f2     # functor to evaluate (2t + 1) + (3t - 5)

    In this example, object 'f3' is a functor of class'SumFunction' that calls
    f1 and f2 and returns their sum.
    """

    def __init__(self, f1, f2):
        """
        :param f1:
            first functor.
        :param f2:
            second functor.
        """
        ComboFunc1.__init__(self, 20, f1, f2)


class DiffFunction(ComboFunc1):
    """Difference of two functions.
    Instances of class DiffFunction evaluate the difference of two supplied
    functors. It is not necessary to explicitly create an instance of
    DiffFunction, since the subtraction operator of the base class is
    overloaded to return a DiffFunction instance.

    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1 - f2     # functor to evaluate (2t + 1) - (3t - 5)

    In this example, object 'f3' is a functor of class'DiffFunction' that
    calls f1 and f2 and returns their difference.
    """

    def __init__(self, f1, f2):
        """
        :param f1:
            first functor.
        :param f2:
            second functor.
        """
        ComboFunc1.__init__(self, 25, f1, f2)

class ProdFunction(ComboFunc1):
    """Product of two functions.  Instances of class ProdFunction
    evaluate the product of two supplied functors.  It is not
    necessary to explicitly create an instance of 'ProdFunction',
    since the multiplication operator of the base class is overloaded
    to return a 'ProdFunction' instance.

    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1 * f2     # functor to evaluate (2t + 1)*(3t - 5)

    In this example, object 'f3' is a functor of class'ProdFunction'
    that calls f1 and f2 and returns their product.
    """
    def __init__(self, f1, f2):
        """
        :param f1:
            first functor.
        :param f2:
            second functor.
        """
        ComboFunc1.__init__(self, 30, f1, f2)


class RatioFunction(ComboFunc1):
    """Ratio of two functions.
    Instances of class RatioFunction evaluate the ratio of two supplied functors.
    It  is not necessary to explicitly create an instance of 'RatioFunction', since
    the division operator of the base class is overloaded to return a RatioFunction
    instance.

    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1 / f2     # functor to evaluate (2t + 1)/(3t - 5)

    In this example, object 'f3' is a functor of class'RatioFunction' that
    calls f1 and f2 and returns their ratio.
    """
    def __init__(self, f1, f2):
        """
        :param f1:
            first functor.
        :param f2:
            second functor.
        """
        ComboFunc1.__init__(self, 40, f1, f2)

class CompositeFunction(ComboFunc1):
    """
    Function of a function.
    Instances of class CompositeFunction evaluate f(g(t)) for two supplied
    functors f and g. It  is not necessary to explicitly create an instance
    of 'CompositeFunction', since the () operator of the base class is
    overloaded to return a CompositeFunction when called with a functor
    argument.

    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1(f2)     # functor to evaluate 2(3t - 5) + 1

    In this example, object 'f3' is a functor of class'CompositeFunction'
    that calls f1 and f2 and returns f1(f2(t)).
    """
    def __init__(self, f1, f2):
        """
        :param f1:
            first functor.
        :param f2:
            second functor.
        """
        ComboFunc1.__init__(self, 60, f1, f2)


class DerivativeFunction(Func1):
    def __init__(self, f):
        self.f = f
        #f._own = 0
        self._own = 1
        self._func_id = _cantera.func_derivative(f.func_id())


def derivative(f):
    """
    Take the derivative of a functor *f*
    """
    return DerivativeFunction(f)
