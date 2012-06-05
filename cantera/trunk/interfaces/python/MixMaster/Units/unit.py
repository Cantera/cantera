import operator

class unit:

    _zero = (0,) * 7
    _negativeOne = (-1, ) * 7

    _labels = ('m', 'kg', 's', 'A', 'K', 'mol', 'cd')


    def __init__(self, value, derivation):
        self.value = value
        self.derivation = derivation
        return


    def __add__(self, other):
        if not self.derivation == other.derivation:
            raise ImcompatibleUnits(self, other)

        return unit(self.value + other.value, self.derivation)


    def __sub__(self, other):
        if not self.derivation == other.derivation:
            raise ImcompatibleUnits(self, other)

        return unit(self.value - other.value, self.derivation)


    def __mul__(self, other):
        if type(other) == type(0) or type(other) == type(0.0):
            return unit(other*self.value, self.derivation)

        value = self.value * other.value
        derivation = tuple(map(operator.add, self.derivation, other.derivation))

        return unit(value, derivation)


    def __div__(self, other):
        if type(other) == type(0) or type(other) == type(0.0):
            return unit(self.value/other, self.derivation)

        value = self.value / other.value
        derivation = tuple(map(operator.sub, self.derivation, other.derivation))

        return unit(value, derivation)


    def __pow__(self, other):
        if type(other) != type(0) and type(other) != type(0.0):
            raise BadOperation

        value = self.value ** other
        derivation = tuple(map(operator.mul, [other]*7, self.derivation))

        return unit(value, derivation)


    def __pos__(self): return self


    def __neg__(self): return unit(-self.value, self.derivation)


    def __abs__(self): return unit(abs(self.value), self.derivation)


    def __invert__(self):
        value = 1./self.value
        derivation = tuple(map(operator.mul, self._negativeOne, self.derivation))
        return unit(value, derivation)


    def __rmul__(self, other):
        return unit.__mul__(self, other)

    def __rdiv__(self, other):
        if type(other) != type(0) and type(other) != type(0.0):
            raise BadOperation(self, other)

        value = other/self.value
        derivation = tuple(map(operator.mul, self._negativeOne, self.derivation))

        return unit(value, derivation)


    def __float__(self):
        return self.value
        #if self.derivation == self._zero: return self.value
        #raise BadConversion(self)


    def __str__(self):
        str = "%g" % self.value
        for i in range(0, 7):
            exponent = self.derivation[i]
            if exponent == 0: continue
            if exponent == 1:
                str = str + " %s" % (self._labels[i])
            else:
                str = str + " %s^%d" % (self._labels[i], exponent)

        return str

dimensionless = unit(1, unit._zero)
