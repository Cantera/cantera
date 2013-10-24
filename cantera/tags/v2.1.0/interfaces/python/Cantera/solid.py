"""
"""

import string
import os

from constants import *
from ThermoPhase import ThermoPhase
from Kinetics import Kinetics
from SolidTransport import SolidTransport
import XML
import _cantera

class Solid(ThermoPhase, Kinetics, SolidTransport):
    """
    """

    def __init__(self, src="", root=None):

        self.ckin = 0
        self._owner = 0
        self.verbose = 1

        # get the 'phase' element
        s = XML.find_XML(src=src, root=root, name="phase")

        # get the equation of state model
        ThermoPhase.__init__(self, xml_phase=s)

        # get the kinetics model
        Kinetics.__init__(self, xml_phase=s, phases=[self])

        SolidTransport.__init__(self, phase=self)

        #self.setState_TP(300.0, OneAtm)


    def __repr__(self):
        return _cantera.phase_report(self._phase_id, self.verbose)


    def __del__(self):
        SolidTransport.__del__(self)
        Kinetics.__del__(self)
        ThermoPhase.__del__(self)
