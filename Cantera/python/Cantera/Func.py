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
    """Base class for functions of one variable."""
    def __init__(self, typ, n, coeffs=[]):
        self.n = n
        self.coeffs = asarray(coeffs,'d')
        self._func_id = ctfunc.func_new(typ, n, self.coeffs)
        
    def __del__(self):
        ctfunc.func_del(self._func_id)

    def __call__(self, t):
        return ctfunc.func_value(self._func_id, t)
    
    def __add__(self, other):
        if type(other) == types.FloatType:
            return SumFunction(self, Const(other))
        return SumFunction(self, other)
    
    def __radd__(self, other):
        if type(other) == types.FloatType:
            return SumFunction(Const(other),self)        
        return SumFunction(other, self)

    def __mul__(self, other):
        return ProdFunction(self, other)

    def __rmul__(self, other):
        return ProdFunction(other, self)
    
    def __div__(self, other):
        return RatioFunction(self, other)

    def __rdiv__(self, other):
        return RatioFunction(other, self)            

    def func_id(self):
        return self._func_id


class Polynomial(Func1):
    """A polynomial. The degree is determined by the number of coefficients
    supplied. Example:
    
    p = Polynomial([1.0, -2.0, 3.0])   # 3t^2 - 2t + 1
    """
    def __init__(self, coeffs=[]):
        Func1.__init__(self, 2, len(coeffs)-1, coeffs)



class Fourier(Func1):
    """
    Fourier series.
    
    f(t) = a[0]/2 + sum_{i=1}^n [a[i]*cos(n*omega*t) + b[i]*sin(n*omega*t)]

    Note that b[0] must be specified for symmetry with 'a', but is not
    used.
    """
    
    def __init__(self, omega, c):
        cc = asarray(c,'d')       
        n, m = cc.shape
        if m <> 2:
            raise CanteraError('provide (a, b) for each term')
        cc[0,1] = omega
        Func1.__init__(self, 1, n-1, ravel(transpose(cc)))


class Arrhenius(Func1):
    def __init__(self, c):
        cc = asarray(c,'d')
        n, m = cc.shape
        if m <> 3:
            raise CanteraError('Three Arrhenius parameters (A, b, E) required.')
        Func1.__init__(self, 3, n, ravel(cc))        



def Const(value):
    return Polynomial([value])



# functions that combine two functions

class SumFunction(Func1):
    def __init__(self, f1, f2):
        self.f1 = f1
        self.f2 = f2
        self.n = -1
        self._func_id = ctfunc.func_newcombo(20, f1.func_id(), f2.func_id())

class ProdFunction(Func1):
    def __init__(self, f1, f2):
        self.f1 = f1
        self.f2 = f2        
        self.n = -1
        self._func_id = ctfunc.func_newcombo(30, f1.func_id(), f2.func_id())

class RatioFunction(Func1):
    def __init__(self, f1, f2):
        self.f1 = f1
        self.f2 = f2        
        self.n = -1
        self._func_id = ctfunc.func_newcombo(40, f1.func_id(), f2.func_id())
    
        
