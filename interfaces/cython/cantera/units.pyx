# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

cdef class Units:
    """
    A representation of the units associated with a dimensional quantity.

    Used for converting quantities between unit systems and checking for dimensional
    consistency. `Units` objects are mainly used within the `UnitSystem` class to
    convert values from a user-specified Unit system to Cantera's base units
    (SI + kmol).

    String representations of units can be written using multiplication, division,
    and exponentiation. Spaces are ignored. Positive, negative, and decimal exponents
    are permitted. Examples are::

        ct.Units("kg*m/s^2")
        ct.Units("J/kmol")
        ct.Units("m*s^-2")
        ct.Units("J/kg/K")
        ct.Units("4184.0 J/kmol")

    Metric prefixes are recognized for all units, for example ``nm``, ``hPa``, ``mg``,
    ``EJ``, ``mL``, ``kcal``. In addition, a numeric scaling factor can be provided.
    """
    def __cinit__(self, name=None, init=True):
        if init and name:
            self.units = CxxUnits(stringify(name))
        elif init:
            self.units = CxxUnits()

    def __repr__(self):
        return f"<Units({pystr(self.units.str())}) at {id(self):0x}>"

    property factor:
        """
        Return the factor for converting from this unit to Cantera's base units.
        """
        def __get__(self):
            return self.units.factor()

    @staticmethod
    cdef copy(CxxUnits other):
        """Copy a C++ Units object to a Python object."""
        cdef Units units = Units(init=False)
        units.units = CxxUnits(other)
        return units

cdef class UnitSystem:
    """
    Unit conversion utility

    Provides functions for converting dimensional values from a given unit system.
    The main use is for converting values specified in input files to Cantera's
    native unit system, which is SI units except for the use of kmol as the base
    unit of quantity, i.e. kilogram, meter, second, kelvin, ampere, and kmol.

    Special functions for converting activation energies allow these values to be
    expressed as either energy per quantity, energy (e.g. eV), or temperature by
    applying a factor of the Avogadro number or the gas constant where needed.

    The default unit system used by Cantera is SI+kmol::

        ct.UnitSystem({
            "length": "m", "mass": "kg", "time": "s",
            "quantity": "kmol", "pressure": "Pa", "energy": "J",
            "temperature": "K", "current": "A", "activation-energy": "J/kmol"})

    A CGS+mol unit system with activation energy units of cal/mol is created as::

        ct.UnitSystem({
            "length": "cm", "mass": "g", "time": "s",
            "quantity": "mol", "pressure": "dyn/cm^2", "energy": "erg",
            "temperature": "K", "current": "A", "activation-energy": "cal/mol"})

    Defaults for dimensions not specified will be left unchanged. Accordingly,
    the default unit system is retrieved as::

        ct.UnitSystem()
    """
    def __cinit__(self, units=None):
        self.unitsystem = CxxUnitSystem()
        if units:
            self.units = units

    def __repr__(self):
        units = f"{self.units}".replace(",", ",\n")
        return f"<UnitSystem at {id(self):0x}> with\n{units}"

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
