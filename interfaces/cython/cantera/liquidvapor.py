# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from . import PureFluid, _cantera


def Water():
    """
    Create a `PureFluid` object using the equation of state for water and the
    `WaterTransport` class for viscosity and thermal conductivity.
    """
    class WaterWithTransport(PureFluid, _cantera.Transport):
        __slots__ = ()

    return WaterWithTransport('liquidvapor.yaml', 'water', transport_model='Water')


def Nitrogen():
    """Create a `PureFluid` object using the equation of state for nitrogen."""
    return PureFluid('liquidvapor.yaml', 'nitrogen')


def Methane():
    """Create a `PureFluid` object using the equation of state for methane."""
    return PureFluid('liquidvapor.yaml', 'methane')


def Hydrogen():
    """Create a `PureFluid` object using the equation of state for hydrogen."""
    return PureFluid('liquidvapor.yaml', 'hydrogen')


def Oxygen():
    """Create a `PureFluid` object using the equation of state for oxygen."""
    return PureFluid('liquidvapor.yaml', 'oxygen')


def Hfc134a():
    """Create a `PureFluid` object using the equation of state for HFC-134a."""
    return PureFluid('liquidvapor.yaml', 'HFC-134a')


def CarbonDioxide():
    """Create a `PureFluid` object using the equation of state for carbon dioxide."""
    return PureFluid('liquidvapor.yaml', 'carbon-dioxide')


def Heptane():
    """Create a `PureFluid` object using the equation of state for heptane."""
    return PureFluid('liquidvapor.yaml', 'heptane')
