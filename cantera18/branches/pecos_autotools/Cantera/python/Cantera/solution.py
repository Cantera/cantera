
import os

from constants import *
from ThermoPhase import ThermoPhase
from Kinetics import Kinetics
from Transport import Transport
from set import setByName
import XML
import _cantera

class Solution(ThermoPhase, Kinetics, Transport):    
    """
    A class for chemically-reacting solutions.

    Instances can be created to represent any type of solution -- a
    mixture of gases, a liquid solution, or a solid solution, for
    example.

    Class Solution derives from classes ThermoPhase, Kinetics, and
    Transport.  It defines very few methods of its own, and is
    provided largely for convenience, so that a single object can be
    used to compute thermodynamic, kinetic, and transport properties
    of a solution. Functions like IdealGasMix and others defined in
    module gases return objects of class Solution.

    """

    def __init__(self, src="", id="", loglevel = 0, debug = 0):

        self.ckin = 0
        self._owner = 0
        self.verbose = 1
        fname = os.path.basename(src)
        ff = os.path.splitext(fname)

        if src:
            root = XML.XML_Node(name = 'doc', src = src,
                                preprocess = 1, debug = debug)

        if id:
            s = root.child(id = id)
            
        else:
            s = root.child(name = "phase")

        self._name = s['id']
        
        # initialize the equation of state
        ThermoPhase.__init__(self, xml_phase=s)

        # initialize the kinetics model
        ph = [self]
        Kinetics.__init__(self, xml_phase=s, phases=ph)

        # initialize the transport model
        Transport.__init__(self, xml_phase=s, phase=self,
                           model = '', loglevel=loglevel)
        
    def __del__(self):
        Transport.__del__(self)
        Kinetics.__del__(self)
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
    
