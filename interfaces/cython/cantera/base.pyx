# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

cdef class _SolutionBase:
    def __cinit__(self, infile='', phaseid='', phases=(), origin=None,
                  source=None, yaml=None, thermo=None, species=(),
                  kinetics=None, reactions=(), electron=None,
                  electron_cross_sections=(), efile='', **kwargs):
        # Shallow copy of an existing Solution (for slicing support)
        cdef _SolutionBase other
        if origin is not None:
            other = <_SolutionBase?>origin

            self.thermo = other.thermo
            self.kinetics = other.kinetics
            self.transport = other.transport
            self.electron = other.electron
            self._thermo = other._thermo
            self._kinetics = other._kinetics
            self._transport = other._transport
            self._electron = other._electron

            self.thermo_basis = other.thermo_basis
            self._selected_species = other._selected_species.copy()
            return

        if infile.endswith('.yml') or infile.endswith('.yaml') or yaml:
            self._init_yaml(infile, phaseid, phases, yaml)
        elif infile or source:
            self._init_cti_xml(infile, phaseid, phases, source)
        elif thermo and species:
            self._init_parts(thermo, species, kinetics, phases, reactions)
        else:
            raise ValueError("Arguments are insufficient to define a phase")

        # Initialization of transport is deferred to Transport.__init__
        self.transport = NULL

        # Initiate electron
        if efile:
            self._init_efile(efile)

        if electron_cross_sections:
            self._init_electron(electron, electron_cross_sections)

        self._selected_species = np.ndarray(0, dtype=np.integer)

    def __init__(self, *args, **kwargs):
        if isinstance(self, Transport):
            assert self.transport is not NULL

    def _init_yaml(self, infile, phaseid, phases, source):
        """
        Instantiate a set of new Cantera C++ objects from a YAML
        phase definition
        """
        cdef CxxAnyMap root
        if infile:
            root = AnyMapFromYamlFile(stringify(infile))
        elif source:
            root = AnyMapFromYamlString(stringify(source))

        phaseNode = root["phases"].getMapWhere("name", stringify(phaseid))

        # Thermo
        if isinstance(self, ThermoPhase):
            self._thermo = newPhase(phaseNode, root)
            self.thermo = self._thermo.get()
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
            self._kinetics = newKinetics(v, phaseNode, root)
            self.kinetics = self._kinetics.get()
        else:
            self.kinetics = NULL

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
            self._thermo.reset(self.thermo)
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
            self._kinetics.reset(self.kinetics)
        else:
            self.kinetics = NULL

    def _init_parts(self, thermo, species, kinetics, phases, reactions):
        """
        Instantiate a set of new Cantera C++ objects based on a string defining
        the model type and a list of Species objects.
        """
        self.thermo = newThermoPhase(stringify(thermo))
        self._thermo.reset(self.thermo)
        self.thermo.addUndefinedElements()
        cdef Species S
        for S in species:
            self.thermo.addSpecies(S._species)
        self.thermo.initThermo()

        if not kinetics:
            kinetics = "none"

        cdef ThermoPhase phase
        cdef Reaction reaction
        if isinstance(self, Kinetics):
            self.kinetics = CxxNewKinetics(stringify(kinetics))
            self._kinetics.reset(self.kinetics)
            self.kinetics.addPhase(deref(self.thermo))
            for phase in phases:
                self.kinetics.addPhase(deref(phase.thermo))
            self.kinetics.init()
            self.kinetics.skipUndeclaredThirdBodies(True)
            for reaction in reactions:
                self.kinetics.addReaction(reaction._reaction)

    def _init_efile(self, efile):
        """
        Instantiate a new Electron object via a yaml file.
        """
        cdef CxxAnyMap root
        root = AnyMapFromYamlFile(stringify(efile))

        if isinstance(self, Electron):
            self._electron = newElectron(root, self.thermo)
            self.electron = self._electron.get()
        else:
            self.electron = NULL

    def _init_electron(self, electron, electron_cross_sections):
        """
        Instantiate a new Electron object.
        """
        cdef ElectronCrossSection ecs
        if isinstance(self, Electron):
            self.electron = newElectron(stringify(electron))
            self._electron.reset(self.electron)
            self.electron.init(self.thermo)
            for ecs in electron_cross_sections:
                self.electron.addElectronCrossSection(ecs._electron_cross_section)

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
            if isinstance(species, (str, int)):
                species = (species,)
            self._selected_species.resize(len(species))
            for i,spec in enumerate(species):
                self._selected_species[i] = self.species_index(spec)

    def __reduce__(self):
        raise NotImplementedError('Solution object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Solution object is not copyable')
