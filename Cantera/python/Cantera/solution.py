"""
"""

#import string
import os

from constants import *
from ThermoPhase import ThermoPhase
from Kinetics import Kinetics
from Transport import Transport
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

        # initialize the equation of state
        ThermoPhase.__init__(self, xml_phase=s)

        # initialize the kinetics model
        ph = [self]
        Kinetics.__init__(self, xml_phase=s, phases=ph)

        # initialize the transport model
        Transport.__init__(self, xml_phase=s, phase=self,
                           model = '', loglevel=0)
        
    def __del__(self):
        Transport.__del__(self)
        Kinetics.__del__(self)
        ThermoPhase.__del__(self)

    def __repr__(self):
        return _cantera.phase_report(self._phase_id, self.verbose)


    
