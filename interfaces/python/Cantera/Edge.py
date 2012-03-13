import string
import os

from constants import *
from SurfacePhase import EdgePhase
from Kinetics import Kinetics
import XML

class Edge(EdgePhase, Kinetics):
    """
    One-dimensional edge between two surfaces.

    Instances of class Edge represent reacting 1D edges between
    between 2D surfaces. Class Edge defines no methods of its
    own. All of its methods derive from either :class:`.EdgePhase` or
    :class:`.Kinetics`.

    Function :func:`.importInterface` should usually be used to build an
    Edge object from a CTI file definition, rather than calling
    the :class:`.Edge` constructor directly.
    """
    def __init__(self, src="", root=None, surfaces=[]):
        """
        :param src:
            CTML or CTI input file name. If more than one phase is
            defined in the file, src should be specified as ``filename#id``
            If the file is not CTML, it will be run through the CTI -> CTML
            preprocessor first.
        :param root:
            If a CTML tree has already been read in that contains
            the definition of this interface, the root of this tree can be
            specified instead of specifying *src*.
        :param phases:
            A list of all objects representing the neighboring
            surface phases which participate in the reaction mechanism.
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
            root = XML.XML_Node(name = 'doc', src = fn, preprocess = 1)

        # If an 'id' tag was specified, find the node in the tree with
        # that tag
        if id:
            s = root.child(id = id)

        # otherwise, find the first element with tag name 'phase'
        # (1D, 2D and 3D phases use the CTML tag name 'phase'
        else:
            s = root.child(name = "phase")

        # build the surface phase
        EdgePhase.__init__(self, xml_phase=s)

        # build the reaction mechanism. This object (representing the
        # surface phase) is added to the end of the list of phases
        Kinetics.__init__(self, xml_phase=s, phases=surfaces+[self])


    def __del__(self):
        """Delete the Edge instance."""
        Kinetics.__del__(self)
        EdgePhase.__del__(self)
