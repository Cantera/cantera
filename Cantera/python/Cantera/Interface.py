
import string
import os

from constants import *
from SurfacePhase import SurfacePhase
from Kinetics import Kinetics
import XML

class Interface(SurfacePhase, Kinetics):    
    """
    Interface objects represent reacting 2D interfaces between bulk 3D phases. Use function
    importInterface to build an Interface object from a CTI file definition, rather than
    calling the Interface constructor directly.
    """
    def __init__(self, src="", root=None, phases=[]):
        
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
        
        # get the 'phase' element
        if src and not root:
            root = XML.XML_Node(name = 'doc', src = fn, preprocess = 1)
            
        if id:
            s = root.child(id = id)            
        else:
            s = root.child(name = "phase")            

        # get the equation of state model
        SurfacePhase.__init__(self, xml_phase=s)

        # get the kinetics model
        Kinetics.__init__(self, xml_phase=s, phases=phases+[self])


    def __del__(self):
        Kinetics.__del__(self)
        SurfacePhase.__del__(self)
    
