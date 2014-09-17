"""Fluids with complete liquid/vapor equations of state..

These functions are defined for convenience only.  They simply call
function 'importPhase' to import the phase definition from file
'liquidvapor.cti' """

from importFromFile import importPhase

import os

from constants import *
from ThermoPhase import ThermoPhase
from set import setByName
import XML
import _cantera


class PureFluid(ThermoPhase):
    """
    A class for chemically-reacting solutions.

    Instances can be created to represent any type of solution -- a
    mixture of gases, a liquid solution, or a solid solution, for
    example.

    """

    def __init__(self, src="", id=""):

        self.ckin = 0
        self._owner = 0
        self.verbose = 1
        fname = os.path.basename(src)
        ff = os.path.splitext(fname)

        if src:
            root = XML.XML_Node(name = 'doc', src = src, preprocess = 1)

        if id:
            s = root.child(id = id)

        else:
            s = root.child(name = "phase")

        self._name = s['id']

        # initialize the equation of state
        ThermoPhase.__init__(self, xml_phase=s)


    def __del__(self):
        ThermoPhase.__del__(self)

    def __repr__(self):
        return _cantera.phase_report(self._phase_id, self.verbose)

    def name(self):
        return self._name

    def set(self, **options):
        """Set various properties.
        T       --- temperature [K]
        P       --- pressure [Pa]
        Rho     --- density [kg/m3]
        V       --- specific volume [m3/kg]
        H       --- specific enthalpy [J/kg]
        U       --- specific internal energy [J/kg]
        S       --- specific entropy [J/kg/K]
        X       --- mole fractions (string or array)
        Y       --- mass fractions (string or array)
        Vapor   --- saturated vapor fraction
        Liquid  --- saturated liquid fraction
        """
        setByName(self, options)

    def critTemperature(self):
        """Critical temperature [K]."""
        return _cantera.thermo_getfp(self._phase_id,50)

    def critPressure(self):
        """Critical pressure [Pa]."""
        return _cantera.thermo_getfp(self._phase_id,51)

    def critDensity(self):
        """Critical density [kg/m3]."""
        return _cantera.thermo_getfp(self._phase_id,52)

    def vaporFraction(self):
        """Vapor fraction."""
        return _cantera.thermo_getfp(self._phase_id,53)

    def setState_Psat(self, p, vaporFraction):
        """Set the state of a saturated liquid/vapor mixture by
        specifying the pressure and vapor fraction."""
        _cantera.thermo_setfp(self._phase_id,8, p, vaporFraction)

    def setState_Tsat(self, t, vaporFraction):
        """Set the state of a saturated liquid/vapor mixture by
        specifying the temperature and vapor fraction."""
        _cantera.thermo_setfp(self._phase_id,7, t, vaporFraction)



def Water():
    return PureFluid('liquidvapor.cti','water')

def Nitrogen():
    return PureFluid('liquidvapor.cti','nitrogen')

def Methane():
    return PureFluid('liquidvapor.cti','methane')

def Hydrogen():
    return PureFluid('liquidvapor.cti','hydrogen')

def Oxygen():
    return PureFluid('liquidvapor.cti','oxygen')

def HFC134a():
    return PureFluid('liquidvapor.cti','hfc134a')

def CarbonDioxide():
    return PureFluid('liquidvapor.cti','carbondioxide')

def Heptane():
    return PureFluid('liquidvapor.cti','heptane')
