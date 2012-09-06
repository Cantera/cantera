cdef class _SolutionBase:
    def __cinit__(self, infile, phaseid='', phases=()):
        rootNode = getCtmlTree(stringify(infile))

        # Get XML data
        cdef XML_Node* phaseNode
        if phaseid:
            phaseNode = rootNode.findID(stringify(phaseid))
        else:
            phaseNode = rootNode.findByName(stringify('phase'))
        if phaseNode is NULL:
            raise ValueError("Couldn't read phase node from XML file")

        # Thermo
        if isinstance(self, ThermoPhase):
            self.thermo = newPhase(deref(phaseNode))
        else:
            self.thermo = NULL

        # Kinetics
        cdef vector[CxxThermoPhase*] v
        cdef _SolutionBase phase

        if isinstance(self, Kinetics):
            v.push_back(self.thermo)
            for phase in phases:
                # adjacent bulk phases for a surface phase
                v.push_back(phase.thermo)
            self.kinetics = newKineticsMgr(deref(phaseNode), v)
        else:
            self.kinetics = NULL

        # Initialization of transport is deferred to Transport.__init__
        self.transport = NULL

    def __init__(self, *args, **kwargs):
        if isinstance(self, Transport):
            assert self.transport is not NULL

    def __dealloc__(self):
        del self.thermo
        del self.kinetics
        del self.transport
