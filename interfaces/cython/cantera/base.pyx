cdef class _SolutionBase:
    def __cinit__(self, infile='', phaseid='', phases=(), origin=None,
                  source=None, thermo=None, species=(), **kwargs):
        # Shallow copy of an existing Solution (for slicing support)
        cdef _SolutionBase other
        if origin is not None:
            other = <_SolutionBase?>origin

            # keep a reference to the parent to prevent the underlying
            # C++ objects from being deleted
            self.parent = other

            self.thermo = other.thermo
            self.kinetics = other.kinetics
            self.transport = other.transport

            self.thermo_basis = other.thermo_basis
            self._selected_species = other._selected_species.copy()
            return

        if infile or source:
            self._init_cti_xml(infile, phaseid, phases, source)
        elif thermo and species:
            self._init_parts(thermo, species)
        else:
            raise ValueError("Arguments are insufficient to define a phase")

        # Initialization of transport is deferred to Transport.__init__
        self.transport = NULL

        self._selected_species = np.ndarray(0, dtype=np.integer)

    def __init__(self, *args, **kwargs):
        if isinstance(self, Transport):
            assert self.transport is not NULL

    def _init_cti_xml(self, infile, phaseid, phases, source):
        """
        Instantiate a set of new Cantera C++ objects from a CTI or XML
        phase definition
        """
        if infile:
            rootNode = CxxGetXmlFile(stringify(infile))
        elif source:
            rootNode = CxxGetXmlFromString(stringify(source))

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

    def _init_parts(self, thermo, species):
        """
        Instantiate a set of new Cantera C++ objects based on a string defining
        the model type and a list of Species objects.
        """
        self.thermo = newThermoPhase(stringify(thermo))
        self.thermo.addUndefinedElements()
        cdef Species S
        for S in species:
            self.thermo.addSpecies(S._species)
        self.thermo.initThermo()

        if isinstance(self, Kinetics):
            # Not yet implemented
            self.kinetics = CxxNewKinetics(stringify("none"))

    def __getitem__(self, selection):
        copy = self.__class__(origin=self)
        if isinstance(selection, slice):
            selection = range(selection.start or 0,
                              selection.stop or self.n_species,
                              selection.step or 1)
        copy.selected_species = selection
        return copy

    property selected_species:
        def __get__(self):
            return list(self._selected_species)
        def __set__(self, species):
            if isinstance(species, (str, unicode, int)):
                species = (species,)
            self._selected_species.resize(len(species))
            for i,spec in enumerate(species):
                self._selected_species[i] = self.species_index(spec)

    def __reduce__(self):
        raise NotImplementedError('Solution object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Solution object is not copyable')

    def __dealloc__(self):
        # only delete the C++ objects if this is the parent object
        if self.parent is None:
            del self.thermo
            del self.kinetics
            del self.transport
