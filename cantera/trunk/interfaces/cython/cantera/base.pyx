cdef class _SolutionBase:
    def __cinit__(self, infile, phaseid=''):
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

        if isinstance(self, Kinetics):
            v.push_back(self.thermo)
            self.kinetics = newKineticsMgr(deref(phaseNode), v)
        else:
            self.kinetics = NULL

        # Transport
        if isinstance(self, Transport):
            self.transport = newDefaultTransportMgr(self.thermo)
        else:
            self.transport = NULL

    def __dealloc__(self):
        del self.thermo
        del self.kinetics
        del self.transport
