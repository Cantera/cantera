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

    def __init__(self, src="", root=None,
                 transport = "", thermo_db = "",
                 transport_db = "", phases=[]):
        
        self.ckin = 0
        self._owner = 0
        self.verbose = 1
        fn = src.split('#')
        id = ""
        if len(fn) > 1:
            id = fn[1]
            fn = fn[0]
            
        fname = os.path.basename(fn)
        ff = os.path.splitext(fname)

        if src and not root:
            root = XML.XML_Node(name = 'doc', src = fn, preprocess = 1)
            
        if id:
            s = root.child(id = id)
        else:
            s = root.child(name = "phase")

        # get the equation of state model
        ThermoPhase.__init__(self, xml_phase=s)

        # get the kinetics model
        ph = [self]+list(phases)
        Kinetics.__init__(self, xml_phase=s, phases=ph)

        Transport.__init__(self, xml_phase=s, phase=self,
                           model = transport, loglevel=4)
        

    def __repr__(self):
        return _cantera.phase_report(self._phase_id, self.verbose)


    def __del__(self):
        Transport.__del__(self)
        Kinetics.__del__(self)
        ThermoPhase.__del__(self)
    
