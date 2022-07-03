# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.
from typing import Dict

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
    def dimensions(self) -> Dict[str, str]:
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
    cdef copy(CxxUnits other):
        """Copy a C++ Units object to a Python object."""
        cdef Units units = Units()
        units.units = CxxUnits(other)
        return units


cdef class UnitSystem:
    """
    Unit system used for YAML input and output.

    The `UnitSystem` class is used to specify dimensional values for a given unit
    system. The main use is for converting values specified in input files to Cantera's
    native unit system, which is SI units except for the use of kmol as the base
    unit of quantity, that is, kilogram, meter, second, kelvin, ampere, and kmol.

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
        self.unitsystem = CxxUnitSystem()
        if units:
            self.units = units

    def __repr__(self):
        return f"<UnitSystem at {id(self):0x}>"

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
