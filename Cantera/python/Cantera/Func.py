"""

The classes in this module are designed to allow constructing
user-defined functions of one variable in Python that can be used with the
Cantera C++ kernel. These classes are mostly shadow classes for
corresponding classes in the C++ kernel.

"""

from Numeric import array, asarray, ravel, shape, transpose
import _cantera
import types


class Func1:
    """A class for functors of one variable.

    A Functor is an object that behaves like a function. Class 'Func1'
    is the base class from which several functor classes derive. These
    classes are designed to be used with the Cantera kernel.  """
    
    def __init__(self, typ, n, coeffs=[]):
        """
        typ - functor type

        n - order

        coeffs - coefficient array
        """
        self.n = n
        self.coeffs = asarray(coeffs,'d')
        self._func_id = _cantera.func_new(typ, n, self.coeffs)
        
    def __del__(self):
        _cantera.func_del(self._func_id)

    def __call__(self, t):
        """Implements function syntax, so that F(t) is equivalent to
        F.value(t)."""
        return _cantera.func_value(self._func_id, t)
    
    def __add__(self, other):
        """Overloads operator '+'

        Returns a new function self(t) + other(t)"""
        if type(other) == types.FloatType:
            return SumFunction(self, Const(other))
        return SumFunction(self, other)
    
    def __radd__(self, other):
        """Overloads operator '+'

        Returns a new function other(t) + self(t)"""        
        if type(other) == types.FloatType:
            return SumFunction(Const(other),self)        
        return SumFunction(other, self)

    def __mul__(self, other):
        """Overloads operator '*'
        
           Return a new function self(t)*other(t)"""        
        return ProdFunction(self, other)

    def __rmul__(self, other):
        """Overloads operator '*'

        Returns a new function other(t)*self(t)"""        
        return ProdFunction(other, self)
    
    def __div__(self, other):
        """Overloads operator '/'

        Returns a new function self(t)/other(t)"""        
        return RatioFunction(self, other)

    def __rdiv__(self, other):
        """Overloads operator '/'

        Returns a new function other(t)/self(t)"""        
        return RatioFunction(other, self)            

    def func_id(self):
        """Return the integer index used internally to access the
        kernel-level object."""
        return self._func_id


class Polynomial(Func1):
    """A polynomial.
    Instances of class 'Polynomial' evaluate
    \f[
    f(t) = \sum_{n = 0}^N a_n t^n.
    \f]
    The coefficients are supplied as a list, beginning with \f$a_N\f$ and ending with \f$a_0\f$.
    >>> p1 = Polynomial([1.0, -2.0, 3.0])   #    3t^2 - 2t + 1
    >>> p2 = Polynomial([6.0, 8.0])         #    8t + 6
    """
    def __init__(self, coeffs=[]):
        """
        coeffs - polynomial coefficients
        """
        Func1.__init__(self, 2, len(coeffs)-1, coeffs)


class Gaussian(Func1):
    """A Gaussian pulse. Instances of class 'Gaussian' evaluate 
    \f[
    f(t) = A \exp[-(t - t_0) / \tau]
    \f]
    where
    \f[
    \tau = \frac{\mbox{FWHM}}{2.0\sqrt{\log(2.0)}}
    \f]
    Here FWHM denotes the full width at half maximum. 
    """
    def __init__(self, A = 0.0, t0 = 0.0, FWHM = 1.0):
        """

        A - Peak value.

        t0 - time at which pulse is centered.

        FWHM - full width at half-maximum.

        """
        coeffs = array([A, t0, FWHM], 'd')
        Func1.__init__(self, 4, 0, coeffs)

        
class Fourier(Func1):
    """Fourier series. Instances of class 'Fourier' evaluate the Fourier series
    \f[
    f(t) = \frac{a_0}{2} + \sum_{n=1}^N [a_n \cos(n\omega t) + b_n \sin(n \omega t)]
    \f]
    where
    \f[
    a_n = \int_{-\pi/\omega}^{\pi/\omega} f(t) \cos(n \omega t) dt
    \f]
    and
    \f[
    b_n = \int_{-\pi/\omega}^{\pi/\omega} f(t) \sin(n \omega t) dt.
    \f]
    The function \f$ f(t) \f$ must be periodic, with period \f$ T = 2\pi/\omega \f$.
    >>> coeffs = [(a0, b0), (a1, b1), (a2, b2)]
    >>> f = Fourier(omega, coeffs)
    Note that b0 must be specified, but is not
    used. The value of b0 is arbitrary.
    """
    def __init__(self, omega, coefficients):
        """
        omega - fundamental frequency [radians/sec].
        
        coefficients - List of (a,b) pairs, beginning with \f$n = 0\f$.
        
        """
        cc = asarray(coefficients,'d')       
        n, m = cc.shape
        if m <> 2:
            raise CanteraError('provide (a, b) for each term')
        cc[0,1] = omega
        Func1.__init__(self, 1, n-1, ravel(transpose(cc)))


class Arrhenius(Func1):
    """Sum of modified Arrhenius terms. Instances of class 'Arrhenius' evaluate
    \f[
    f(T) = \sum_{i=1}^n A_n T^{b_n}\exp(-E_n/T)
    \f]
    Example:

    >>> f = Arrhenius([(a0, b0, e0), (a1, b1, e1)])
    
    """
    def __init__(self, coefficients):
        """
        coefficients - sequence of \f$(A, b, E)\f$ triplets.
        
        """
        cc = asarray(coefficients,'d')
        n, m = cc.shape
        if m <> 3:
            raise CanteraError('Three Arrhenius parameters (A, b, E) required.')
        Func1.__init__(self, 3, n, ravel(cc))        



def Const(value):
    """Constant function.
    >>> f = Const(4.0)  # evaluates f(t) = 4.0.
    """
    return Polynomial([value])


class PeriodicFunction(Func1):
    """Converts a function into a periodic function with period T."""
    def __init__(self, func, T):
        """
        func - initial non-periodic function

        T - period [s]
        """
        Func1.__init__(self, 50, func.func_id(), array([T],'d'))
        

# functions that combine two functions

class SumFunction(Func1):
    """Sum of two functions.
    Instances of class SumFunction evaluate the sum of two supplied functors.
    It is not necessary to explicitly create an instance of SumFunction, since
    the addition operator of the base class is overloaded to return a SumFunction
    instance.
    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1 + f2     # functor to evaluate (2t + 1) + (3t - 5)
    In this example, object 'f3' is a functor of class'SumFunction' that calls f1 and f2
    and returns their sum.
    """
    
    def __init__(self, f1, f2):
        """
        f1 - first functor.
        
        f2 - second functor.
        """
        self.f1 = f1
        self.f2 = f2
        self.n = -1
        self._func_id = _cantera.func_newcombo(20, f1.func_id(), f2.func_id())

class ProdFunction(Func1):
    """Product of two functions.
    Instances of class ProdFunction evaluate the product of two supplied functors.
    It is not necessary to explicitly create an instance of 'ProdFunction', since 
    the multiplication operator of the base class is overloaded to return a 'ProdFunction'
    instance.
    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1 * f2     # functor to evaluate (2t + 1)*(3t - 5)
    In this example, object 'f3' is a functor of class'ProdFunction' that calls f1 and f2
    and returns their product.
    """    
    def __init__(self, f1, f2):
        """
        f1 - first functor.
        
        f2 - second functor.
        """        
        self.f1 = f1
        self.f2 = f2        
        self.n = -1
        self._func_id = _cantera.func_newcombo(30, f1.func_id(), f2.func_id())

class RatioFunction(Func1):
    """Ratio of two functions.
    Instances of class RatioFunction evaluate the ratio of two supplied functors.
    It  is not necessary to explicitly create an instance of 'RatioFunction', since
    the division operator of the base class is overloaded to return a RatioFunction
    instance.
    >>> f1 = Polynomial([2.0, 1.0])
    >>> f2 = Polynomial([3.0, -5.0])
    >>> f3 = f1 / f2     # functor to evaluate (2t + 1)/(3t - 5)
    In this example, object 'f3' is a functor of class'RatioFunction' that calls f1 and f2
    and returns their ratio.
    """        
    def __init__(self, f1, f2):
        """
        f1 - first functor.
        
        f2 - second functor.
        """                
        self.f1 = f1
        self.f2 = f2        
        self.n = -1
        self._func_id = _cantera.func_newcombo(40, f1.func_id(), f2.func_id())
    
        
