import string
import os

from constants import *
from SurfacePhase import SurfacePhase, EdgePhase
from Kinetics import Kinetics
import XML

class Interface(SurfacePhase, Kinetics):
    """
    Two-dimensional interfaces.

    Instances of class Interface represent reacting 2D interfaces
    between bulk 3D phases. Class Interface defines no methods of its
    own. All of its methods derive from either :class:`.SurfacePhase` or
    :class:`.Kinetics`.

    Function :func:`.importInterface` should usually be used to build an
    Interface object from a CTI file definition, rather than calling
    the Interface constructor directly.
    """
    def __init__(self, src="", root=None, phases=[], debug = 0):
        """
        :param src:
            CTML or CTI input file name. If more than one phase is
            defined in the file, src should be specified as ``filename#id``
            If the file is not CTML, it will be run through the CTI -> CTML
            preprocessor first.
        :param root:
            If a CTML tree has already been read in that contains the
            definition of this interface, the root of this tree can be
            specified instead of specifying *src*.
        :param phases:
            A list of all objects representing the neighboring phases which
            participate in the reaction mechanism.
        """
        self.ckin = 0
        self._owner = 0
        self.verbose = 1

        # src has the form '<filename>#<id>'
        fn = src.split('#')
        id = ""
        if len(fn) > 1:
            id = fn[1]
            fn = fn[0]

        # read in the root element of the tree if not building from
        # an already-built XML tree. Enable preprocessing if the film
        # is a .cti file instead of XML.
        if src and not root:
            root = XML.XML_Node(name = 'doc', src = fn, preprocess = 1, debug = debug)

        # If an 'id' tag was specified, find the node in the tree with
        # that tag
        if id:
            s = root.child(id = id)

        # otherwise, find the first element with tag name 'phase'
        # (both 2D and 3D phases use the CTML tag name 'phase'
        else:
            s = root.child(name = "phase")

        # build the surface phase
        SurfacePhase.__init__(self, xml_phase=s)

        # build the reaction mechanism. This object (representing the
        # surface phase) is added to the end of the list of phases
        Kinetics.__init__(self, xml_phase=s, phases=phases+[self])


    def __del__(self):
        """Delete the Interface instance."""
        Kinetics.__del__(self)
        SurfacePhase.__del__(self)
