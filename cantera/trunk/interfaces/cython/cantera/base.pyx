cdef class _SolutionBase:
    def __cinit__(self, infile='', phaseid='', phases=(), source=None):
        # Shallow copy of an existing Solution (for slicing support)
        cdef _SolutionBase other
        if source is not None:
            other = <_SolutionBase?>source

            # keep a reference to the parent to prevent the underlying
            # C++ objects from being deleted
            self.parent = other

            self.thermo = other.thermo
            self.kinetics = other.kinetics
            self.transport = other.transport

            self.thermoBasis = other.thermoBasis
            self._selectedSpecies = other._selectedSpecies.copy()
            return

        # Instantiate a set of new Cantera C++ objects
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

        self._selectedSpecies = np.ndarray(0, dtype=np.integer)

    def __init__(self, *args, **kwargs):
        if isinstance(self, Transport):
            assert self.transport is not NULL

    def __getitem__(self, selection):
        copy = self.__class__(source=self)
        if isinstance(selection, slice):
            selection = range(selection.start or 0,
                              selection.stop or self.nSpecies,
                              selection.step or 1)
        copy.selectedSpecies = selection
        return copy

    property selectedSpecies:
        def __get__(self):
            return list(self._selectedSpecies)
        def __set__(self, species):
            if isinstance(species, (str, int)):
                species = (species,)
            self._selectedSpecies.resize(len(species))
            for i,spec in enumerate(species):
                self._selectedSpecies[i] = self.speciesIndex(spec)

    def __dealloc__(self):
        # only delete the C++ objects if this is the parent object
        if self.parent is None:
            del self.thermo
            del self.kinetics
            del self.transport
