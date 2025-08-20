# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections.abc import Sequence as _Sequence
import numbers as _numbers
import numpy as np

from ._utils cimport *

cdef class Units:
    """
    A representation of the units associated with a dimensional quantity.

    This class is a light-weight interface to internal Cantera capabilities that are
    used for converting quantities between unit systems and checking for dimensional
    consistency. Internally, `Units` objects are mainly used within the `UnitSystem`
    class to convert values from a user-specified Unit system to Cantera's base units
    (SI + kmol).

    The Python API handles display of `Units` that do not require a conversion factor,
    with other functions not enabled.
    """
    def __cinit__(self, name=None):
        if name:
            self.units = CxxUnits(stringify(name), True)

    def __repr__(self) -> str:
        return f"<Units({str(self)}) at {id(self):0x}>"

    def __str__(self) -> str:
        return pystr(self.units.str())

    def dimension(self, primary: str) -> float:
        """The dimension of the given unit component.

        .. versionadded:: 3.0

        :param primary:
            A string with the desired unit component. One of ``"mass"``,
            ``"length"``, ``"time"``, ``"temperature"``, ``"current"``, or
            ``"quantity"``.
        """
        return self.units.dimension(stringify(primary))

    @property
    def dimensions(self) -> dict[str, float]:
        """A dictionary of the primary unit components to their dimensions.

        .. versionadded:: 3.0
        """
        dimensions = ("mass", "length", "time", "temperature", "current", "quantity")
        return {d: self.dimension(d) for d in dimensions}

    @property
    def factor(self) -> float:
        """The factor required to convert from this unit to Cantera's base units.

        .. versionadded:: 3.0
        """
        return self.units.factor()

    @staticmethod
    cdef Units copy(CxxUnits other):
        """Copy a C++ Units object to a Python object."""
        cdef Units units = Units()
        units.units = CxxUnits(other)
        return units


cdef class UnitStack:
    def __cinit__(self):
        self.stack = CxxUnitStack(CxxUnits())

    @staticmethod
    cdef UnitStack copy(const CxxUnitStack& other):
        """Copy a C++ UnitStack object to a Python object."""
        cdef UnitStack stack = UnitStack()
        stack.stack = CxxUnitStack(other)
        return stack

    def product(self):
        cdef CxxUnits units = self.stack.product()
        return Units.copy(units)

    def join(self, exponent):
        self.stack.join(exponent)


cdef class UnitSystem:
    """
    Unit system used for YAML input and output.

    The `UnitSystem` class is used to specify dimensional values for a given unit
    system. The main use is for converting values specified in input files to Cantera's
    native unit system, which is SI units except for the use of kmol as the base
    unit of quantity, that is, kilogram, meter, second, kelvin, ampere, and kmol.

    Generally, this class is used indirectly, through methods interacting with `AnyMap`
    objects such as `ExtensibleRate.set_parameters` and `ExtensibleRate.get_parameters`.

    The default unit system used by Cantera is SI+kmol::

        ct.UnitSystem({
            "length": "m", "mass": "kg", "time": "s",
            "quantity": "kmol", "pressure": "Pa", "energy": "J",
            "temperature": "K", "current": "A", "activation-energy": "J / kmol"})

    A CGS+mol unit system with activation energy units of cal/mol is created as::

        ct.UnitSystem({
            "length": "cm", "mass": "g", "time": "s",
            "quantity": "mol", "pressure": "dyn / cm^2", "energy": "erg",
            "temperature": "K", "current": "A", "activation-energy": "cal / mol"})

    Defaults for dimensions not specified will be left unchanged. Accordingly,
    the default unit system is retrieved as::

        ct.UnitSystem()

    """
    def __cinit__(self, units=None):
        self._unitsystem.reset(new CxxUnitSystem())
        self.unitsystem = self._unitsystem.get()
        if units:
            self.units = units

    def __repr__(self):
        return f"<UnitSystem at {id(self):0x}>"

    cdef _set_unitSystem(self, shared_ptr[CxxUnitSystem] units):
        self._unitsystem = units
        self.unitsystem = self._unitsystem.get()

    def defaults(self):
        return self.unitsystem.defaults()

    property units:
        """
        Units used by the unit system
        """
        def __get__(self):
            cdef stdmap[string, string] cxxunits = self.unitsystem.defaults()
            cdef pair[string, string] item
            return {pystr(item.first): pystr(item.second) for item in cxxunits}
        def __set__(self, units):
            cdef stdmap[string, string] cxxunits
            for dimension, unit in units.items():
                cxxunits[stringify(dimension)] = stringify(unit)
            self.unitsystem.setDefaults(cxxunits)

    def convert_to(self, quantity, dest):
        """
        Convert *quantity* to the units defined by *dest*, using this `UnitSystem` to
        define the default units of *quantity*. *quantity* can be one of the following:

        - A number, for example ``3.14``
        - A "quantity string" containing a number and a dimension string, separated by
          a space. For example, ``"3.14 kmol/m^3"``
        - A NumPy array of either numeric values or quantity strings as described above
        - A list, tuple, or other sequence of any shape containing numeric values or
          quantity strings, For example ``("3000 mm", 3.14, "12 cm")``

        *dest* can be a string or `Units` object specifying the destination units.
        """
        cdef double value
        cdef CxxAnyValue val_units
        if isinstance(quantity, str):
            val_units = python_to_anyvalue(quantity)
            if isinstance(dest, str):
                return self.unitsystem.convert(val_units, stringify(dest))
            elif isinstance(dest, Units):
                return self.unitsystem.convert(val_units, (<Units>dest).units)
            else:
                raise TypeError("'dest' must be a string or 'Units' object")
        elif isinstance(quantity, _numbers.Real):
            value = quantity
            if isinstance(dest, str):
                return self.unitsystem.convertTo(value, stringify(dest))
            elif isinstance(dest, Units):
                return self.unitsystem.convertTo(value, (<Units>dest).units)
            else:
                raise TypeError("'dest' must be a string or 'Units' object")
        elif isinstance(quantity, np.ndarray):
            return np.vectorize(lambda item: self.convert_to(item, dest))(quantity)
        elif isinstance(quantity, _Sequence):
            return [self.convert_to(item, dest) for item in quantity]
        else:
            raise TypeError("'quantity' must be either a string or a number")

    def convert_activation_energy_to(self, quantity, str dest):
        """
        Convert *quantity* to the activation energy units defined by *dest*, using this
        `UnitSystem` to define the default units of *quantity*. *quantity* can be one of
        the following:

        - A number, for example ``3.14``
        - A "quantity string" containing a number and a dimension string, separated by
          a space. For example, ``"3.14 J/kmol"``
        - A NumPy array of either numeric values or quantity strings as described above
        - A list, tuple, or other sequence of any shape containing numeric values or
          quantity strings, For example ``("30 kcal/mol", 3.14, "12000 K")``

        *dest* can be a string or `Units` object specifying the destination units, which
        must be interpretable as unit of energy per unit quantity (for example, J/kmol),
        energy (for example, eV) or temperature (K).
        """
        cdef double value
        cdef CxxAnyValue val_units
        if isinstance(quantity, str):
            val_units = python_to_anyvalue(quantity)
            return self.unitsystem.convertActivationEnergy(val_units, stringify(dest))
        elif isinstance(quantity, _numbers.Real):
            value = quantity
            return self.unitsystem.convertActivationEnergyTo(value, stringify(dest))
        elif isinstance(quantity, np.ndarray):
            return np.vectorize(
                lambda item: self.convert_activation_energy_to(item, dest))(quantity)
        elif isinstance(quantity, _Sequence):
            return [self.convert_activation_energy_to(item, dest) for item in quantity]
        else:
            raise TypeError("'quantity' must be either a string or a number")

    def convert_rate_coeff_to(self, quantity, dest):
        """
        Convert a *quantity* representing a rate coefficient to the units defined by
        *dest*, using this `UnitSystem` to define the default units of *quantity*.
        Behaves similar to`convert_to` but with special handling for the case of
        standalone rate constants, where the destination units are not actually known,
        and where the units may be specified using `Units` or `UnitStack` objects.
        """
        cdef Units dest_units
        if isinstance(dest, Units):
            dest_units = dest
        elif isinstance(dest, UnitStack):
            dest_units = dest.product()
        elif isinstance(dest, str):
            dest_units = Units(dest)
        else:
            raise TypeError("'dest' must be a string, Units, or UnitStack object")

        cdef CxxAnyValue val_units
        if isinstance(quantity, (str, _numbers.Real)):
            val_units = python_to_anyvalue(quantity)
            return self.unitsystem.convertRateCoeff(val_units, dest_units.units)
        elif isinstance(quantity, np.ndarray):
            return np.vectorize(
                lambda q: self.convert_rate_coeff_to(q, dest_units))(quantity)
        elif isinstance(quantity, _Sequence):
            return [self.convert_rate_coeff_to(q, dest_units) for q in quantity]
        else:
            raise TypeError("'quantity' must be either a string or a number")
