from ._cantera import *

class Quantity(object):
    """
    A class representing a specific quantity of a `Solution`.
    """
    def __init__(self, phase, mass=None, moles=None, constant='UV'):
        self.state = phase.TDY
        self._phase = phase

        if mass is not None:
            self.mass = mass
        elif moles is not None:
            self.moles = moles
        else:
            self.mass = 1.0

        assert constant in ('TP','TV','HP','SP','SV','UV')
        self.constant = constant

    @property
    def phase(self):
        self._phase.TDY = self.state
        return self._phase

    @property
    def moles(self):
        return self.mass / self.phase.mean_molecular_weight

    @moles.setter
    def moles(self, n):
        self.mass = n * self.phase.mean_molecular_weight

    @property
    def volume(self):
        return self.mass * self.phase.volume_mass

    @property
    def int_energy(self):
        return self.mass * self.phase.int_energy_mass

    @property
    def enthalpy(self):
        return self.mass * self.phase.enthalpy_mass

    @property
    def entropy(self):
        return self.mass * self.phase.entropy_mass

    @property
    def gibbs(self):
        return self.mass * self.phase.gibbs_mass

    def __imul__(self, other):
        self.mass *= other
        return self

    def __mul__(self, other):
        return Quantity(self.phase, mass=self.mass * other)

    def __rmul__(self, other):
        return Quantity(self.phase, mass=self.mass * other)

    def __iadd__(self, other):
        assert(self.constant == other.constant)
        a1,b1 = getattr(self.phase, self.constant)
        a2,b2 = getattr(other.phase, self.constant)
        m = self.mass + other.mass
        a = (a1 * self.mass + a2 * other.mass) / m
        b = (b1 * self.mass + b2 * other.mass) / m
        self._phase.Y = (self.Y * self.mass + other.Y * other.mass) / m
        setattr(self._phase, self.constant, (a,b))
        self.state = self._phase.TDY
        self.mass = m
        return self

    def __add__(self, other):
        newquantity = Quantity(self.phase, mass=self.mass, constant=self.constant)
        newquantity += other
        return newquantity

# Synonyms for total properties
Quantity.V = Quantity.volume
Quantity.U = Quantity.int_energy
Quantity.H = Quantity.enthalpy
Quantity.S = Quantity.entropy
Quantity.G = Quantity.gibbs

# Add properties to act as pass-throughs for attributes of class Solution
def _prop(attr):
    def getter(self):
        return getattr(self.phase, attr)

    def setter(self, value):
        setattr(self.phase, attr, value)
        self.state = self._phase.TDY

    return property(getter, setter, doc=getattr(Solution, attr).__doc__)

for _attr in dir(Solution):
    if _attr.startswith('_') or _attr in Quantity.__dict__:
        continue
    else:
        setattr(Quantity, _attr, _prop(_attr))
