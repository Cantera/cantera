# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pint import UnitRegistry
from pint import set_application_registry as set_application_registry

from .solution import (
    Q_,
    CarbonDioxide,
    Heptane,
    Hfc134a,
    Hydrogen,
    Methane,
    Nitrogen,
    Oxygen,
    PureFluid,
    Solution,
    Water,
    units,
)

cantera_units_registry: UnitRegistry

__all__: list[str] = [
    "Q_",
    "CarbonDioxide",
    "Heptane",
    "Hfc134a",
    "Hydrogen",
    "Methane",
    "Nitrogen",
    "Oxygen",
    "PureFluid",
    "Solution",
    "Water",
    "cantera_units_registry",
    "set_application_registry",
    "units",
]
