"""
"""

import string
import os

from constants import *
from SurfacePhase import SurfacePhase
from Kinetics import Kinetics
import XML

class Interface(SurfacePhase, Kinetics):    
    """
    ...
    """

    def __init__(self, src="", root=None, phases=[]):
        
        self.ckin = 0
        self._owner = 0
        self.verbose = 1

        fn = string.split(src,'#')
        fn = fn[0]
        fname = os.path.basename(fn)
        ff = os.path.splitext(fname)
        
        # get the 'phase' element
        s = XML.find_XML(src=src, root=root, name="phase")

        # get the equation of state model
        SurfacePhase.__init__(self, xml_phase=s)

        # get the kinetics model
        Kinetics.__init__(self, xml_phase=s, phases=[self]+phases)


    def __del__(self):
        Kinetics.__del__(self)
        SurfacePhase.__del__(self)
    
