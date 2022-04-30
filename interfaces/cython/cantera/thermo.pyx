# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import warnings
import weakref
import numbers as _numbers

cdef enum ThermoBasisType:
    mass_basis = 0
    molar_basis = 1

ctypedef CxxPlasmaPhase* CxxPlasmaPhasePtr

class ThermoModelMethodError(Exception):
    """Exception raised for an invalid method used by a thermo model

    :param thermo_model:
        The thermo model used by class `ThermoPhase`

    """

    def __init__(self, thermo_model):
        self.thermo_model = thermo_model
        super().__init__(f"This method is invalid for {self.thermo_model}")


cdef class Species:
    """
    A class which stores data about a single chemical species that may be
    needed to add it to a `Solution` or `Interface` object (and to the
    underlying `ThermoPhase` and `Transport` objects).

    :param name:
        A string giving the name of the species, such as ``'CH4'``
    :param composition:
        The elemental composition of the species, given either as a dict or a
        composition string, such as ``{'C':1, 'H':4}`` or ``'C:1, H:4'``.
    :param charge:
        The electrical charge, in units of the elementary charge. Default 0.0.
    :param size:
        The effective size [m] of the species. Default 1.0.
    :param init:
        Used internally when wrapping :ct:`Species` objects returned from C++

    Example: creating an ideal gas phase with a single species::

        ch4 = ct.Species('CH4', 'C:1, H:4')
        ch4.thermo = ct.ConstantCp(300, 1000, 101325,
                                   (300, -7.453347e7, 1.865912e5, 3.576053e4))
        tran = ct.GasTransportData()
        tran.set_customary_units('nonlinear', 3.75, 141.40, 0.0, 2.60, 13.00)
        ch4.transport = tran
        gas = ct.Solution(thermo='IdealGas', species=[ch4])

    The static methods `fromYaml`, `fromCti`, `fromXml`, `listFromFile`,
    `listFromYaml`, `listFromCti`, and `listFromXml` can be used to create
    `Species` objects from existing definitions in the CTI or XML formats.
    Either of the following will produce a list of 53 `Species` objects
    containing the species defined in the GRI 3.0 mechanism::

        S = ct.Species.list_from_file("gri30.yaml")

        import pathlib
        S = ct.Species.listFromYaml(
            pathlib.Path('path/to/gri30.yaml').read_text(),
            section='species')

    """
    def __cinit__(self, *args, init=True, **kwargs):
        if init:
            self._species.reset(new CxxSpecies())
            self.species = self._species.get()

    def __init__(self, name=None, composition=None, charge=None, size=None,
                 *args, init=True, **kwargs):
        if not init:
            return

        if name is not None:
            self.species.name = stringify(name)

        if composition is not None:
            self.species.composition = comp_map(composition)

        if charge is not None:
            self.species.charge = charge

        if size is not None:
            self.species.size = size

    cdef _assign(self, shared_ptr[CxxSpecies] other):
        self._species = other
        self.species = self._species.get()

    @staticmethod
    def fromCti(text):
        """
        Create a Species object from its CTI string representation.

        .. deprecated:: 2.5

            The CTI input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_species = CxxGetSpecies(deref(CxxGetXmlFromString(stringify(text))))
        assert cxx_species.size() == 1, cxx_species.size()
        species = Species(init=False)
        species._assign(cxx_species[0])
        return species

    @staticmethod
    def fromXml(text):
        """
        Create a Species object from its XML string representation.

        .. deprecated:: 2.5

            The XML input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_species = CxxNewSpecies(deref(CxxGetXmlFromString(stringify(text))))
        species = Species(init=False)
        species._assign(cxx_species)
        return species

    @staticmethod
    def fromYaml(text):
        """
        Create a Species object from its YAML string representation.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by `Reaction.from_yaml`.
        """
        warnings.warn("Class method 'fromYaml' is renamed to 'from_yaml' "
            "and will be removed after Cantera 2.6.", DeprecationWarning)

        return Species.from_yaml(text)

    @staticmethod
    def from_yaml(text):
        """
        Create a `Species` object from its YAML string representation.
        """
        cxx_species = CxxNewSpecies(AnyMapFromYamlString(stringify(text)))
        species = Species(init=False)
        species._assign(cxx_species)
        return species

    @staticmethod
    def from_dict(data):
        """
        Create a `Species` object from a dictionary corresponding to its YAML
        representation.

        :param data:
            A dictionary corresponding to the YAML representation.
        """
        cdef CxxAnyMap any_map = dict_to_anymap(data)
        cxx_species = CxxNewSpecies(any_map)
        species = Species(init=False)
        species._assign(cxx_species)
        return species

    @staticmethod
    def listFromFile(filename, section='species'):
        """
        Create a list of Species objects from all of the species defined in a
        YAML, CTI or XML file. For YAML files, return species from the section
        ``section``.

        Directories on Cantera's input file path will be searched for the
        specified file.

        In the case of an XML file, the ``<species>`` nodes are assumed to be
        children of the ``<speciesData>`` node in a document with a ``<ctml>``
        root node, as in the XML files produced by conversion from CTI files.

        .. deprecated:: 2.5

            The CTI and XML input formats are deprecated and will be removed in
            Cantera 3.0.

        .. deprecated:: 2.6

            To be removed after Cantera 2.6. Replaced by 'Species.list_from_file'.
        """
        warnings.warn("Static method 'listFromFile' is renamed to 'list_from_file'."
            " The old name will be removed after Cantera 2.6.", DeprecationWarning)

        if filename.lower().split('.')[-1] in ('yml', 'yaml'):
            root = AnyMapFromYamlFile(stringify(filename))
            cxx_species = CxxGetSpecies(root[stringify(section)])
        else:
            cxx_species = CxxGetSpecies(deref(CxxGetXmlFile(stringify(filename))))

        species = []
        for a in cxx_species:
            b = Species(init=False)
            b._assign(a)
            species.append(b)
        return species

    @staticmethod
    def list_from_file(filename, section="species"):
        """
        Create a list of Species objects from all of the species defined in the section
        ``section`` of a YAML file. Directories on Cantera's input file path will be
        searched for the specified file.
        """
        root = AnyMapFromYamlFile(stringify(filename))
        species = []
        for a in CxxGetSpecies(root[stringify(section)]):
            b = Species(init=False)
            b._assign(a)
            species.append(b)
        return species

    @staticmethod
    def listFromXml(text):
        """
        Create a list of Species objects from all the species defined in an XML
        string. The ``<species>`` nodes are assumed to be children of the
        ``<speciesData>`` node in a document with a ``<ctml>`` root node, as in
        the XML files produced by conversion from CTI files.

        .. deprecated:: 2.5

            The XML input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_species = CxxGetSpecies(deref(CxxGetXmlFromString(stringify(text))))
        species = []
        for a in cxx_species:
            b = Species(init=False)
            b._assign(a)
            species.append(b)
        return species

    @staticmethod
    def listFromCti(text):
        """
        Create a list of Species objects from all the species defined in a CTI
        string.

        .. deprecated:: 2.5

            The CTI input format is deprecated and will be removed in Cantera 3.0.
        """
        # Currently identical to listFromXml since get_XML_from_string is able
        # to distinguish between CTI and XML.
        cxx_species = CxxGetSpecies(deref(CxxGetXmlFromString(stringify(text))))
        species = []
        for a in cxx_species:
            b = Species(init=False)
            b._assign(a)
            species.append(b)
        return species

    @staticmethod
    def listFromYaml(text, section=None):
        """
        Create a list of Species objects from all the species defined in a YAML string.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by `Reaction.list_from_yaml`.
        """
        warnings.warn("Class method 'listFromYaml' is renamed to 'list_from_yaml' "
            "and will be removed after Cantera 2.6.", DeprecationWarning)
        return Species.list_from_yaml(text, section)

    @staticmethod
    def list_from_yaml(text, section=None):
        """
        Create a list of Species objects from all the species defined in a YAML
        string. If ``text`` is a YAML mapping, the ``section`` name of the list
        to be read must be specified. If ``text`` is a YAML list, no ``section``
        name should be supplied.
        """
        root = AnyMapFromYamlString(stringify(text))

        # ``items`` is the pseudo-key used to access a list when it is at the
        # top level of a YAML document
        cxx_species = CxxGetSpecies(root[stringify(section or "items")])
        species = []
        for a in cxx_species:
            b = Species(init=False)
            b._assign(a)
            species.append(b)
        return species

    property name:
        """ The name of the species. """
        def __get__(self):
            return pystr(self.species.name)

    property composition:
        """
        A dict containing the elemental composition of the species. Keys are
        element names; values are the corresponding atomicities.
        """
        def __get__(self):
            return comp_map_to_dict(self.species.composition)

    property charge:
        """
        The electrical charge on the species, in units of the elementary charge.
        """
        def __get__(self):
            return self.species.charge

    property size:
        """ The effective size [m] of the species. """
        def __get__(self):
            return self.species.size

    property thermo:
        """
        Get/Set the species reference-state thermodynamic data, as an instance
        of class `SpeciesThermo`.
        """
        def __get__(self):
            if self.species.thermo.get() != NULL:
                return wrapSpeciesThermo(self.species.thermo)
            else:
                return None

        def __set__(self, SpeciesThermo spthermo):
            self.species.thermo = spthermo._spthermo

    property transport:
        """
        Get/Set the species transport parameters, as an instance of class
        `GasTransportData`.
        """
        def __get__(self):
            if self.species.transport.get() != NULL:
                data = GasTransportData(init=False)
                data._assign(self.species.transport)
                return data
            else:
                return None
        def __set__(self, GasTransportData tran):
            self.species.transport = tran._data

    property input_data:
        """
        Get input data defining this Species, along with any user-specified data
        provided with its input (YAML) definition.
        """
        def __get__(self):
            cdef CxxThermoPhase* phase = self._phase.thermo if self._phase else NULL
            return anymap_to_dict(self.species.parameters(phase))

    def update_user_data(self, data):
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.species.input.update(dict_to_anymap(data), False)

    def clear_user_data(self):
        """
        Clear all saved input data, so that the data given by `input_data` or
        `Solution.write_yaml` will only include values generated by Cantera based on
        the current object state.
        """
        self.species.input.clear()

    def __repr__(self):
        return '<Species {}>'.format(self.name)


cdef class ThermoPhase(_SolutionBase):
    """
    A phase with an equation of state.

    Class `ThermoPhase` may be used to represent the intensive thermodynamic
    state of a phase of matter, which might be a gas, liquid, or solid.

    Class `ThermoPhase` is not usually instantiated directly. It is used
    as a base class for classes `Solution` and `Interface`.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if 'source' not in kwargs:
            self.thermo_basis = mass_basis
        # In composite objects, the ThermoPhase constructor needs to be called first
        # to prevent instantiation of stand-alone 'Kinetics' or 'Transport' objects.
        # The following is used as a sentinel.
        self._references = weakref.WeakKeyDictionary()
        # validate plasma phase
        self._enable_plasma = False
        if dynamic_cast[CxxPlasmaPhasePtr](self.thermo):
            self._enable_plasma = True
            self.plasma = <CxxPlasmaPhase*>self.thermo

    property thermo_model:
        """
        Return thermodynamic model describing phase.
        """
        def __get__(self):
            return pystr(self.thermo.type())

    property phase_of_matter:
        """
        Get the thermodynamic phase (gas, liquid, etc.) at the current conditions.
        """
        def __get__(self):
            return pystr(self.thermo.phaseOfMatter())

    def report(self, show_thermo=True, float threshold=1e-14):
        """
        Generate a report describing the thermodynamic state of this phase. To
        print the report to the terminal, simply call the phase object. The
        following two statements are equivalent::

        >>> phase()
        >>> print(phase.report())
        """
        return pystr(self.thermo.report(bool(show_thermo), threshold))

    def __call__(self, *args, **kwargs):
        print(self.report(*args, **kwargs))

    property is_pure:
        """
        Returns true if the phase represents a pure (fixed composition) substance
        """
        def __get__(self):
            return self.thermo.isPure()

    property has_phase_transition:
        """
        Returns true if the phase represents a substance with phase transitions
        """
        def __get__(self):
            return self.thermo.hasPhaseTransition()

    property is_compressible:
        """
        Returns true if the density of the phase is an independent variable defining
        the thermodynamic state of a substance
        """
        def __get__(self):
            return self.thermo.isCompressible()

    property _native_state:
        """
        Default properties defining a state
        """
        def __get__(self):
            cdef pair[string, size_t] item
            native = {pystr(item.first): item.second for item in self.thermo.nativeState()}
            return tuple([i for i, j in sorted(native.items(), key=lambda kv: kv[1])])

    property _full_states:
        """
        Sets of parameters which set the full thermodynamic state
        """
        def __get__(self):
            states = self.thermo.fullStates()
            states = [pystr(s) for s in states]
            return {frozenset(k): k for k in states}

    property _partial_states:
        """
        Sets of parameters which set a valid partial thermodynamic state
        """
        def __get__(self):
            states = self.thermo.partialStates()
            states = [pystr(s) for s in states]
            return {frozenset(k): k for k in states}

    property basis:
        """
        Determines whether intensive thermodynamic properties are treated on a
        ``mass`` (per kg) or ``molar`` (per kmol) basis. This affects the values
        returned by the properties `h`, `u`, `s`, `g`, `v`, `density`, `cv`,
        and `cp`, as well as the values used with the state-setting properties
        such as `HPX` and `UV`.
        """
        def __get__(self):
            if self.thermo_basis == mass_basis:
                return 'mass'
            else:
                return 'molar'

        def __set__(self, value):
            if value == 'mass':
                self.thermo_basis = mass_basis
            elif value == 'molar':
                self.thermo_basis = molar_basis
            else:
                raise ValueError("Valid choices are 'mass' or 'molar'."
                                 " Got {!r}.".format(value))

    cdef double _mass_factor(self):
        """ Conversion factor from current basis to kg """
        if self.thermo_basis == molar_basis:
            return self.thermo.meanMolecularWeight()
        else:
            return 1.0

    cdef double _mole_factor(self):
        """ Conversion factor from current basis to moles """
        if self.thermo_basis == mass_basis:
            return 1.0/self.thermo.meanMolecularWeight()
        else:
            return 1.0

    def equilibrate(self, XY, solver='auto', double rtol=1e-9,
                    int max_steps=1000, int max_iter=100, int estimate_equil=0,
                    int log_level=0):
        """
        Set to a state of chemical equilibrium holding property pair
        ``XY`` constant.

        :param XY:
            A two-letter string, which must be one of the set::

                ['TP','TV','HP','SP','SV','UV']

        :param solver:
            Specifies the equilibrium solver to use. May be one of the following:

            * ``'element_potential'`` - a fast solver using the element potential
              method
            * ``'gibbs'`` - a slower but more robust Gibbs minimization solver
            * ``'vcs'`` - the VCS non-ideal equilibrium solver
            * ``'auto'`` - The element potential solver will be tried first, then
              if it fails the Gibbs solver will be tried.
        :param rtol:
            The relative error tolerance.
        :param max_steps:
            The maximum number of steps in composition to take to find a converged
            solution.
        :param max_iter:
            For the Gibbs minimization solver, this specifies the number of
            outer iterations on T or P when some property pair other
            than TP is specified.
        :param estimate_equil:
            Integer indicating whether the solver should estimate its own
            initial condition. If 0, the initial mole fraction vector in the
            `ThermoPhase` object is used as the initial condition. If 1, the
            initial mole fraction vector is used if the element abundances are
            satisfied. If -1, the initial mole fraction vector is thrown out,
            and an estimate is formulated.
        :param log_level:
            Set to a value greater than 0 to write diagnostic output.
        """
        self.thermo.equilibrate(stringify(XY.upper()), stringify(solver), rtol,
                                max_steps, max_iter, estimate_equil, log_level)

    ####### Composition, species, and elements ########

    property n_elements:
        """Number of elements."""
        def __get__(self):
            return self.thermo.nElements()

    cpdef int element_index(self, element) except *:
        """
        The index of element ``element``, which may be specified as a string or
        an integer. In the latter case, the index is checked for validity and
        returned. If no such element is present, an exception is thrown.
        """
        if isinstance(element, (str, bytes)):
            index = self.thermo.elementIndex(stringify(element))
        elif isinstance(element, (int, float)):
            index = <int>element
        else:
            raise TypeError("'element' must be a string or a number."
                            " Got {!r}.".format(element))

        if not 0 <= index < self.n_elements:
            raise ValueError('No such element {!r}.'.format(element))

        return index

    def element_name(self, m):
        """Name of the element with index ``m``."""
        return pystr(self.thermo.elementName(m))

    property element_names:
        """A list of all the element names."""
        def __get__(self):
            return [self.element_name(m) for m in range(self.n_elements)]

    def atomic_weight(self, m):
        """Atomic weight [kg/kmol] of element ``m``"""
        return self.thermo.atomicWeight(self.element_index(m))

    property atomic_weights:
        """Array of atomic weight [kg/kmol] for each element in the mixture."""
        def __get__(self):
            return np.array([self.thermo.atomicWeight(m) for m in range(self.n_elements)])

    property n_species:
        """Number of species."""
        def __get__(self):
            return self.thermo.nSpecies()

    property n_selected_species:
        """
        Number of species selected for output (by slicing of Solution object)
        """
        def __get__(self):
            return self._selected_species.size or self.n_species

    def species_name(self, k):
        """Name of the species with index ``k``."""
        return pystr(self.thermo.speciesName(k))

    property species_names:
        """A list of all the species names."""
        def __get__(self):
            if self._selected_species.size:
                indices = self._selected_species
            else:
                indices = range(self.n_species)
            return [self.species_name(k) for k in indices]

    cpdef int species_index(self, species) except *:
        """
        The index of species ``species``, which may be specified as a string or
        an integer. In the latter case, the index is checked for validity and
        returned. If no such species is present, an exception is thrown.
        """
        if isinstance(species, (str, bytes)):
            index = self.thermo.speciesIndex(stringify(species))
        elif isinstance(species, (int, float)):
            index = <int>species
        else:
            raise TypeError("'species' must be a string or a number."
                            " Got {!r}.".format(species))

        if not 0 <= index < self.n_species:
            raise ValueError('No such species {!r}.'.format(species))

        return index

    property case_sensitive_species_names:
        """Enforce case-sensitivity for look up of species names"""
        def __get__(self):
            return self.thermo.caseSensitiveSpecies()
        def __set__(self, val):
            self.thermo.setCaseSensitiveSpecies(bool(val))

    def species(self, k=None):
        """
        Return the `Species` object for species ``k``, where ``k`` is either the
        species index or the species name. If ``k`` is not specified, a list of
        all species objects is returned. Changes to this object do not affect
        the `ThermoPhase` or `Solution` object until the `modify_species`
        function is called.
        """
        if k is None:
            return [self.species(i) for i in range(self.n_species)]

        s = Species(init=False)
        s._phase = self
        if isinstance(k, (str, bytes)):
            s._assign(self.thermo.species(stringify(k)))
        elif isinstance(k, (int, float)):
            s._assign(self.thermo.species(<int>k))
        else:
            raise TypeError("Argument must be a string or a number."
                            " Got {!r}.".format(k))
        return s

    def modify_species(self, k, Species species):
        """
        Modify the thermodynamic data associated with a species. The species name,
        elemental composition, and type of thermo parameterization must be unchanged.
        """
        self.thermo.modifySpecies(k, species._species)
        if self.kinetics:
            self.kinetics.invalidateCache()

    def add_species(self, Species species):
        """
        Add a new species to this phase. Missing elements will be added
        automatically.
        """
        if self._references:
            raise CanteraError('Cannot add species to ThermoPhase object if it'
                ' is linked to a Reactor, Domain1D (flame), or Mixture object.')
        self.thermo.addUndefinedElements()
        self.thermo.addSpecies(species._species)
        species._phase = self
        self.thermo.initThermo()
        if self.kinetics:
            self.kinetics.invalidateCache()

    def add_species_alias(self, name, alias):
        """
        Add the alternate species name ``alias`` for an original species ``name``.
        """
        self.thermo.addSpeciesAlias(stringify(name), stringify(alias))

    def find_isomers(self, comp):
        """
        Find species/isomers matching a composition specified by ``comp``.
        """

        if isinstance(comp, dict):
            iso = self.thermo.findIsomers(comp_map(comp))
        elif isinstance(comp, (str, bytes)):
            iso = self.thermo.findIsomers(stringify(comp))
        else:
            raise CanteraError('Invalid composition')

        return [pystr(b) for b in iso]

    def n_atoms(self, species, element):
        """
        Number of atoms of element ``element`` in species ``species``. The element
        and species may be specified by name or by index.

        >>> phase.n_atoms('CH4','H')
        4
        """
        return self.thermo.nAtoms(self.species_index(species),
                                  self.element_index(element))

    cdef np.ndarray _getArray1(self, thermoMethod1d method):
        cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.n_species)
        method(self.thermo, &data[0])
        if self._selected_species.size:
            return data[self._selected_species]
        else:
            return data

    cdef void _setArray1(self, thermoMethod1d method, values) except *:
        cdef np.ndarray[np.double_t, ndim=1] data

        values = np.squeeze(values)
        if values.ndim == 0:
            values = values[np.newaxis] # corner case for single-species phases

        if len(values) == self.n_species:
            data = np.ascontiguousarray(values, dtype=np.double)
        elif len(values) == len(self._selected_species) != 0:
            data = np.zeros(self.n_species, dtype=np.double)
            for i,k in enumerate(self._selected_species):
                data[k] = values[i]
        else:
            msg = "Got {}. Expected {}".format(len(values), self.n_species)
            if len(self._selected_species):
                msg += ' or {}'.format(len(self._selected_species))
            raise ValueError('Array has incorrect length. ' + msg + '.')
        method(self.thermo, &data[0])

    property molecular_weights:
        """Array of species molecular weights (molar masses) [kg/kmol]."""
        def __get__(self):
            return self._getArray1(thermo_getMolecularWeights)

    property charges:
        """Array of species charges [elem. charge]."""
        def __get__(self):
            return self._getArray1(thermo_getCharges)

    property mean_molecular_weight:
        """The mean molecular weight (molar mass) [kg/kmol]."""
        def __get__(self):
            return self.thermo.meanMolecularWeight()

    property Y:
        """
        Get/Set the species mass fractions. Can be set as an array, as a dictionary,
        or as a string. Always returns an array::

            >>> phase.Y = [0.1, 0, 0, 0.4, 0, 0, 0, 0, 0.5]
            >>> phase.Y = {'H2':0.1, 'O2':0.4, 'AR':0.5}
            >>> phase.Y = 'H2:0.1, O2:0.4, AR:0.5'
            >>> phase.Y
            array([0.1, 0, 0, 0.4, 0, 0, 0, 0, 0.5])
        """
        def __get__(self):
            return self._getArray1(thermo_getMassFractions)
        def __set__(self, Y):
            if isinstance(Y, (str, bytes)):
                self.thermo.setMassFractionsByName(stringify(Y))
            elif isinstance(Y, dict):
                self.thermo.setMassFractionsByName(comp_map(Y))
            else:
                self._setArray1(thermo_setMassFractions, Y)

    property X:
        """
        Get/Set the species mole fractions. Can be set as an array, as a dictionary,
        or as a string. Always returns an array::

            >>> phase.X = [0.1, 0, 0, 0.4, 0, 0, 0, 0, 0.5]
            >>> phase.X = {'H2':0.1, 'O2':0.4, 'AR':0.5}
            >>> phase.X = 'H2:0.1, O2:0.4, AR:0.5'
            >>> phase.X
            array([0.1, 0, 0, 0.4, 0, 0, 0, 0, 0.5])
        """
        def __get__(self):
            return self._getArray1(thermo_getMoleFractions)
        def __set__(self, X):
            if isinstance(X, (str, bytes)):
                self.thermo.setMoleFractionsByName(stringify(X))
            elif isinstance(X, dict):
                self.thermo.setMoleFractionsByName(comp_map(X))
            else:
                self._setArray1(thermo_setMoleFractions, X)

    property concentrations:
        """Get/Set the species concentrations [kmol/m^3]."""
        def __get__(self):
            return self._getArray1(thermo_getConcentrations)
        def __set__(self, C):
            self._setArray1(thermo_setConcentrations, C)

    def __composition_to_array(self, comp, basis):
        """take a mixture composition in mole or mass fraction as string,
        dict or array and return array (for internal use)"""
        if (isinstance(comp, str) and ':' not in comp
            and comp in self.species_names):
            comp += ':1.0'

        original_state = self.state

        if basis == 'mole':
            self.TPX = None, None, comp
            arr = np.copy(self.X)
        elif basis == 'mass':
            self.TPY = None, None, comp
            arr = np.copy(self.Y)
        else:
            raise ValueError("basis must either be 'mass' or mole'.")

        self.state = original_state
        return arr

    def set_equivalence_ratio(self, phi, fuel, oxidizer, basis="mole", *, diluent=None,
                              fraction=None):
        """
        Set the composition to a mixture of ``fuel`` and ``oxidizer`` at the
        specified equivalence ratio ``phi``, holding temperature and pressure
        constant. Considers the oxidation of C to CO2, H to H2O and S to SO2.
        Other elements are assumed not to participate in oxidation (that is,
        N ends up as N2). The ``basis`` determines the fuel and oxidizer
        compositions: ``basis='mole'`` means mole fractions (default),
        ``basis='mass'`` means mass fractions. The fuel/oxidizer mixture can be
        be diluted by a ``diluent`` based on a mixing ``fraction``. The amount of
        diluent is quantified as a fraction of fuel, oxidizer or the fuel/oxidizer
        mixture. For more information, see `Python example
        <https://cantera.org/examples/python/thermo/equivalenceRatio.py.html>`_ ::

            >>> gas.set_equivalence_ratio(0.5, 'CH4', 'O2:1.0, N2:3.76', basis='mole')
            >>> gas.mass_fraction_dict()
            {'CH4': 0.02837633052851681, 'N2': 0.7452356312613029, 'O2': 0.22638803821018036}
            >>> gas.set_equivalence_ratio(1.2, {'NH3':0.8, 'CO':0.2}, 'O2:1.0', basis='mole')
            >>> gas.mass_fraction_dict()
            {'CO': 0.14784006249290754, 'NH3': 0.35956645545401045, 'O2': 0.49259348205308207}

        :param phi:
            Equivalence ratio
        :param fuel:
            Fuel species name or mole/mass fractions as string, array, or dict.
        :param oxidizer:
            Oxidizer species name or mole/mass fractions as a string, array, or dict.
        :param basis:
            Determines if ``fuel`` and ``oxidizer`` are given in mole
            fractions (``basis='mole'``) or mass fractions (``basis='mass'``)
        :param: diluent:
            Optional parameter. Required if dilution is used. Specifies the composition
            of the diluent in mole/mass fractions as a string, array or dict
        :param: fraction:
            Optional parameter. Dilutes the fuel/oxidizer mixture with the diluent
            according to ``fraction``. Fraction can refer to the fraction of diluent in
            the  mixture (for example ``fraction="diluent:0.7`` will create a mixture
            with 30 % fuel/oxidizer and 70 % diluent), the fraction of fuel in the
            mixture (for example ``fraction="fuel:0.1" means that the mixture contains
            10 % fuel. The amount of oxidizer is determined from the equivalence ratio
            and the remaining mixture is the diluent) or fraction of oxidizer in the
            mixture (for example ``fraction="oxidizer:0.1")``. The fraction itself is
            interpreted as mole or mass fraction based on ``basis``. The diluent is not
            considered in the computation of the equivalence ratio. Default is no
            dilution or ``fraction=None``. May be given as string or dictionary (for
            example ``fraction={"fuel":0.7})``
        """
        cdef np.ndarray[np.double_t, ndim=1] fuel_comp = np.ascontiguousarray(
                self.__composition_to_array(fuel, basis), dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] ox_comp = np.ascontiguousarray(
                self.__composition_to_array(oxidizer, basis), dtype=np.double)

        self.thermo.setEquivalenceRatio(phi, &fuel_comp[0], &ox_comp[0],
                                        ThermoBasis.mass if basis == "mass"
                                                         else ThermoBasis.molar)

        if (fraction is None) != (diluent is None):
            raise ValueError("If dilution is used, both 'fraction' and 'diluent' "
                             "parameters are required.")

        if fraction is None:
            return

        if isinstance(fraction, str):
            fraction_dict = comp_map_to_dict(parseCompString(stringify(fraction)))
        elif isinstance(fraction, dict):
            fraction_dict = fraction
        else:
            raise ValueError("The fraction argument must be given as string or "
                             "dictionary.")

        if len(fraction_dict) != 1:
            raise ValueError("Invalid format for the fraction. Must be provided for "
                             "example as fraction='fuel:0.1'")

        fraction_type  = list(fraction_dict.keys())[0]
        fraction_value = float(list(fraction_dict.values())[0])

        if fraction_value < 0 or fraction_value > 1:
            raise ValueError("The fraction must be between 0 and 1")

        if fraction_type not in ["fuel", "oxidizer", "diluent"]:
            raise ValueError("The fraction must specify 'fuel', 'oxidizer' or "
                             "'diluent'")

        cdef np.ndarray[np.double_t, ndim=1] diluent_comp = np.ascontiguousarray(
                self.__composition_to_array(diluent, basis), dtype=np.double)

        # this function changes the composition and fixes temperature and pressure
        T, P = self.T, self.P
        # if 'fraction' is specified for diluent, just scale the mass or mole fractions
        # of the fuel/oxidizer mixture accordingly
        if fraction_type == "diluent":
            if basis == "mole":
                X_fuelox = self.X
                self.X = diluent_comp
                self.X = (1.0 - fraction_value) * X_fuelox + fraction_value * self.X
            else:
                Y_fuelox = self.Y
                self.Y = diluent_comp
                self.Y = (1.0 - fraction_value) * Y_fuelox + fraction_value * self.Y
            self.TP = T, P
            return

        # get the mixture fraction before scaling / diluent addition
        Z_fuel = self.mixture_fraction(fuel, oxidizer, basis)

        if Z_fuel == 0.0 and fraction_type == "fuel":
            raise ValueError("No fuel in the fuel/oxidizer mixture")

        if Z_fuel == 1.0 and fraction_type == "oxidzer":
            raise ValueError("No oxidizer in the fuel/oxidizer mixture")

        if basis == "mass": # for mass basis, it is straight forward
            if fraction_type == "fuel":
                Z = Z_fuel
            else:  # oxidizer
                Z = 1.0 - Z_fuel
            if fraction_value > Z:
                raise ValueError(f"The {fraction_type} fraction after dilution cannot "
                                 "be higher than {fraction_type} fraction in the "
                                 "original mixture.")
            Y_mix = self.Y
            self.Y = diluent_comp
            factor = fraction_value / Z
            self.Y = factor*Y_mix + (1.0 - factor) * self.Y
        else:
            # convert mass based mixture fraction to molar one, Z = kg fuel / kg mixture
            X_mix = self.X
            M_mix = self.mean_molecular_weight
            self.X = fuel_comp
            M_fuel = self.mean_molecular_weight
            Z_fuel_mole = Z_fuel * M_mix / M_fuel # mol fuel / mol mixture
            if fraction_type == "fuel":
                Z = Z_fuel_mole
            else: # oxidizer
                Z = 1.0 - Z_fuel_mole
            if fraction_value > Z:
                raise ValueError(f"The {fraction_type} fuel or oxidizer fraction after "
                                 "dilution cannot be higher than {fraction_type} "
                                 "fraction in the original mixture.")
            self.X = diluent_comp
            factor = fraction_value / Z
            self.X = factor * X_mix + (1.0 - factor) * self.X
        self.TP = T, P

    def set_mixture_fraction(self, mixture_fraction, fuel, oxidizer, basis='mole'):
        """
        Set the composition to a mixture of ``fuel`` and ``oxidizer`` at the
        specified mixture fraction ``mixture_fraction`` (kg fuel / kg mixture), holding
        temperature and pressure constant. Considers the oxidation of C to CO2,
        H to H2O and S to SO2. Other elements are assumed not to participate in
        oxidation (that is, N ends up as N2). The ``basis`` determines the composition
        of fuel and oxidizer: ``basis='mole'`` (default) means mole fractions,
        ``basis='mass'`` means mass fractions::

            >>> gas.set_mixture_fraction(0.5, 'CH4', 'O2:1.0, N2:3.76')
            >>> gas.mass_fraction_dict()
            {'CH4': 0.5, 'N2': 0.38350014242997776, 'O2': 0.11649985757002226}
            >>> gas.set_mixture_fraction(0.5, {'NH3':0.8, 'CO':0.2}, 'O2:1.0')
            >>> gas.mass_fraction_dict()
            {'CO': 0.145682068778996, 'NH3': 0.354317931221004, 'O2': 0.5}

        :param mixture_fraction:
            Mixture fraction (kg fuel / kg mixture)
        :param fuel:
            Fuel species name or mole/mass fractions as string, array, or dict.
        :param oxidizer:
            Oxidizer species name or mole/mass fractions as a string, array, or
            dict.
        :param basis: determines if ``fuel`` and ``oxidizer`` are given in mole
            fractions (``basis='mole'``) or mass fractions (``basis='mass'``)
        """
        cdef np.ndarray[np.double_t, ndim=1] f = \
                np.ascontiguousarray(self.__composition_to_array(fuel, basis), dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] o = \
                np.ascontiguousarray(self.__composition_to_array(oxidizer, basis), dtype=np.double)

        self.thermo.setMixtureFraction(mixture_fraction, &f[0], &o[0], ThermoBasis.mass if basis == 'mass' else ThermoBasis.molar)

    def equivalence_ratio(self, fuel=None, oxidizer=None, basis="mole",
                          include_species=None):
        """
        Get the equivalence ratio of the current mixture, which is a
        conserved quantity. Considers the oxidation of C to CO2, H to H2O
        and S to SO2. Other elements are assumed not to participate in oxidation
        (that is, N ends up as N2). If fuel and oxidizer are not specified, the
        equivalence ratio is computed from the available oxygen and the
        required oxygen for complete oxidation. The ``basis`` determines the
        composition of fuel and oxidizer: ``basis='mole'`` (default) means mole
        fractions, ``basis='mass'`` means mass fractions. Additionally, a
        list of species can be provided with ``include_species``. This means that
        only these species are considered for the computation of the equivalence
        ratio. For more information, see `Python example
        <https://cantera.org/examples/python/thermo/equivalenceRatio.py.html>`_ ::

            >>> gas.set_equivalence_ratio(0.5, fuel='CH3:0.5, CH3OH:.5, N2:0.125', oxidizer='O2:0.21, N2:0.79, NO:0.01')
            >>> gas.equivalence_ratio(fuel='CH3:0.5, CH3OH:.5, N2:0.125', oxidizer='O2:0.21, N2:0.79, NO:0.01')
            0.5

        :param fuel:
            Fuel species name or mole/mass fractions as string, array, or dict.
        :param oxidizer:
            Oxidizer species name or mole/mass fractions as a string, array, or dict.
        :param basis:
            Determines if ``fuel`` and ``oxidizer`` are given in mole fractions
            (``basis="mole"``) or mass fractions (``basis="mass"``)
        :param include_species:
            List of species names (optional). Only these species are considered for the
            computation of the equivalence ratio. By default, all species are considered
        """
        if include_species is not None:
            # remove unwanted species temporarily
            Y = np.zeros(self.n_species)
            indices = [self.species_index(s) for s in include_species]
            for k in indices:
                Y[k] = self.Y[k]
            T_orig, P_orig, Y_orig = self.T, self.P, self.Y
            self.Y = Y

        if fuel is None and oxidizer is None:
            phi = self.thermo.equivalenceRatio()
            if include_species is not None:
                self.TPY = T_orig, P_orig, Y_orig
            return phi

        cdef np.ndarray[np.double_t, ndim=1] f = np.ascontiguousarray(
                self.__composition_to_array(fuel, basis), dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] o = np.ascontiguousarray(
                self.__composition_to_array(oxidizer, basis), dtype=np.double)

        phi = self.thermo.equivalenceRatio(&f[0], &o[0],
                                           ThermoBasis.mass if basis=="mass"
                                                            else ThermoBasis.molar)
        if include_species is not None:
            self.TPY = T_orig, P_orig, Y_orig
        return phi

    def mixture_fraction(self, fuel, oxidizer, basis='mole', element="Bilger"):
        """
        Get the mixture fraction of the current mixture (kg fuel / (kg oxidizer + kg fuel))
        This is a quantity that is conserved after oxidation. Considers the
        oxidation of C to CO2, H to H2O and S to SO2. Other elements are assumed
        not to participate in oxidation (that is, N ends up as N2).
        The ``basis`` determines the composition of fuel and oxidizer:
        ``basis="mole"`` (default) means mole fractions, ``basis="mass"`` means mass fractions.
        The mixture fraction can be computed from a single element (for example, carbon
        with ``element="C"``) or from all elements, which is the Bilger mixture
        fraction (``element="Bilger"``). The Bilger mixture fraction is computed by default::

            >>> gas.set_mixture_fraction(0.5, 'CH3:0.5, CH3OH:.5, N2:0.125', 'O2:0.21, N2:0.79, NO:0.01')
            >>> gas.mixture_fraction('CH3:0.5, CH3OH:.5, N2:0.125', 'O2:0.21, N2:0.79, NO:0.01')
            0.5

        :param fuel:
            Fuel species name or mole/mass fractions as string, array, or dict.
        :param oxidizer:
            Oxidizer species name or mole/mass fractions as a string, array, or
            dict.
        :param basis:
            Determines if ``fuel`` and ``oxidizer`` are given in mole
            fractions (``basis='mole'``) or mass fractions (``basis='mass'``)
        :param element:
            Computes the mixture fraction from the specified elemental
            mass fraction (given by element name or element index) or as
            the Bilger mixture fraction (default)
        """
        cdef np.ndarray[np.double_t, ndim=1] f = \
                np.ascontiguousarray(self.__composition_to_array(fuel, basis), dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] o = \
                np.ascontiguousarray(self.__composition_to_array(oxidizer, basis), dtype=np.double)

        if isinstance(element, (str, bytes)):
            e_name = element
        else:
            e_name = self.element_name(self.element_index(element))

        return self.thermo.mixtureFraction(&f[0], &o[0], ThermoBasis.mass if basis=='mass' else ThermoBasis.molar, stringify(e_name))

    def stoich_air_fuel_ratio(self, fuel, oxidizer, basis='mole'):
        """
        Get the stoichiometric air to fuel ratio (kg oxidizer / kg fuel). Considers the
        oxidation of C to CO2, H to H2O and S to SO2. Other elements are assumed
        not to participate in oxidation (that is, N ends up as N2).
        The ``basis`` determines the composition of fuel and oxidizer: ``basis='mole'`` (default)
        means mole fractions, ``basis='mass'`` means mass fractions::

            >>> gas.set_mixture_fraction(0.5, 'CH3:0.5, CH3OH:.5, N2:0.125', 'O2:0.21, N2:0.79, NO:0.01')
            >>> gas.stoich_air_fuel_ratio('CH3:0.5, CH3OH:.5, N2:0.125', 'O2:0.21, N2:0.79, NO:0.01')
            8.148040722239438

        :param fuel:
            Fuel species name or mole/mass fractions as string, array, or dict.
        :param oxidizer:
            Oxidizer species name or mole/mass fractions as a string, array, or
            dict.
        :param basis:
            Determines if ``fuel`` and ``oxidizer`` are given in mole
            fractions (``basis='mole'``) or mass fractions (``basis='mass'``)

        """
        cdef np.ndarray[np.double_t, ndim=1] f = \
                np.ascontiguousarray(self.__composition_to_array(fuel, basis), dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] o = \
                np.ascontiguousarray(self.__composition_to_array(oxidizer, basis), dtype=np.double)

        return self.thermo.stoichAirFuelRatio(&f[0], &o[0], ThermoBasis.mass if basis=='mass' else ThermoBasis.molar)

    def elemental_mass_fraction(self, m):
        r"""
        Get the elemental mass fraction :math:`Z_{\mathrm{mass},m}` of element
        :math:`m` as defined by:

        .. math:: Z_{\mathrm{mass},m} = \sum_k \frac{a_{m,k} M_m}{M_k} Y_k

        with :math:`a_{m,k}` being the number of atoms of element :math:`m` in
        species :math:`k`, :math:`M_m` the atomic weight of element :math:`m`,
        :math:`M_k` the molecular weight of species :math:`k`, and :math:`Y_k`
        the mass fraction of species :math:`k`::

            >>> phase.elemental_mass_fraction('H')
            1.0

        :param m:
            Base element, may be specified by name or by index.
        """
        return self.thermo.elementalMassFraction(self.element_index(m))

    def elemental_mole_fraction(self, m):
        r"""
        Get the elemental mole fraction :math:`Z_{\mathrm{mole},m}` of element
        :math:`m` (the number of atoms of element m divided by the total number
        of atoms) as defined by:

        .. math:: Z_{\mathrm{mole},m} = \frac{\sum_k a_{m,k} X_k}
                                             {\sum_k \sum_j a_{j,k} X_k}

        with :math:`a_{m,k}` being the number of atoms of element :math:`m` in
        species :math:`k`, :math:`\sum_j` being a sum over all elements, and
        :math:`X_k` being the mole fraction of species :math:`k`::

            >>> phase.elemental_mole_fraction('H')
            1.0

        :param m:
            Base element, may be specified by name or by index.
        """
        return self.thermo.elementalMoleFraction(self.element_index(m))

    def set_unnormalized_mass_fractions(self, Y):
        """
        Set the mass fractions without normalizing to force ``sum(Y) == 1.0``.
        Useful primarily when calculating derivatives with respect to ``Y[k]`` by
        finite difference.
        """
        cdef np.ndarray[np.double_t, ndim=1] data
        if len(Y) == self.n_species:
            data = np.ascontiguousarray(Y, dtype=np.double)
        else:
            raise ValueError("Array has incorrect length."
                 " Got {}, expected {}.".format(len(Y), self.n_species))
        self.thermo.setMassFractions_NoNorm(&data[0])

    def set_unnormalized_mole_fractions(self, X):
        """
        Set the mole fractions without normalizing to force ``sum(X) == 1.0``.
        Useful primarily when calculating derivatives with respect to ``X[k]``
        by finite difference.
        """
        cdef np.ndarray[np.double_t, ndim=1] data
        if len(X) == self.n_species:
            data = np.ascontiguousarray(X, dtype=np.double)
        else:
            raise ValueError("Array has incorrect length."
                " Got {}, expected {}.".format(len(X), self.n_species))
        self.thermo.setMoleFractions_NoNorm(&data[0])

    def mass_fraction_dict(self, double threshold=0.0):
        """
        Return a dictionary giving the mass fraction for each species by name where the
        mass fraction is greater than ``threshold``.
        """
        cdef pair[string,double] item
        Y = self.thermo.getMassFractionsByName(threshold)
        return {pystr(item.first):item.second for item in Y}

    def mole_fraction_dict(self, double threshold=0.0):
        """
        Return a dictionary giving the mole fraction for each species by name where the
        mole fraction is greater than ``threshold``.
        """
        cdef pair[string,double] item
        X = self.thermo.getMoleFractionsByName(threshold)
        return {pystr(item.first):item.second for item in X}

    ######## Read-only thermodynamic properties ########

    property P:
        """Pressure [Pa]."""
        def __get__(self):
            return self.thermo.pressure()

    property T:
        """Temperature [K]."""
        def __get__(self):
            return self.thermo.temperature()

    property density:
        """Density [kg/m^3 or kmol/m^3] depending on `basis`."""
        def __get__(self):
            return self.thermo.density() / self._mass_factor()

    property density_mass:
        """(Mass) density [kg/m^3]."""
        def __get__(self):
            return self.thermo.density()

    property density_mole:
        """Molar density [kmol/m^3]."""
        def __get__(self):
            return self.thermo.molarDensity()

    property v:
        """Specific volume [m^3/kg or m^3/kmol] depending on `basis`."""
        def __get__(self):
            return self._mass_factor() / self.thermo.density()

    property volume_mass:
        """Specific volume [m^3/kg]."""
        def __get__(self):
            return 1.0 / self.thermo.density()

    property volume_mole:
        """Molar volume [m^3/kmol]."""
        def __get__(self):
            return self.thermo.molarVolume()

    property u:
        """Internal energy in [J/kg or J/kmol]."""
        def __get__(self):
            return self.thermo.intEnergy_mole() * self._mole_factor()

    property int_energy_mole:
        """Molar internal energy [J/kmol]."""
        def __get__(self):
            return self.thermo.intEnergy_mole()

    property int_energy_mass:
        """Specific internal energy [J/kg]."""
        def __get__(self):
            return self.thermo.intEnergy_mass()

    property h:
        """Enthalpy [J/kg or J/kmol] depending on `basis`."""
        def __get__(self):
            return self.thermo.enthalpy_mole() * self._mole_factor()

    property enthalpy_mole:
        """Molar enthalpy [J/kmol]."""
        def __get__(self):
            return self.thermo.enthalpy_mole()

    property enthalpy_mass:
        """Specific enthalpy [J/kg]."""
        def __get__(self):
            return self.thermo.enthalpy_mass()

    property s:
        """Entropy [J/kg/K or J/kmol/K] depending on `basis`."""
        def __get__(self):
            return self.thermo.entropy_mole() * self._mole_factor()

    property entropy_mole:
        """Molar entropy [J/kmol/K]."""
        def __get__(self):
            return self.thermo.entropy_mole()

    property entropy_mass:
        """Specific entropy [J/kg/K]."""
        def __get__(self):
            return self.thermo.entropy_mass()

    property g:
        """Gibbs free energy [J/kg or J/kmol] depending on `basis`."""
        def __get__(self):
            return self.thermo.gibbs_mole() * self._mole_factor()

    property gibbs_mole:
        """Molar Gibbs free energy [J/kmol]."""
        def __get__(self):
            return self.thermo.gibbs_mole()

    property gibbs_mass:
        """Specific Gibbs free energy [J/kg]."""
        def __get__(self):
            return self.thermo.gibbs_mass()

    property cv:
        """
        Heat capacity at constant volume [J/kg/K or J/kmol/K] depending on
        `basis`.
        """
        def __get__(self):
            return self.thermo.cv_mole() * self._mole_factor()

    property cv_mole:
        """Molar heat capacity at constant volume [J/kmol/K]."""
        def __get__(self):
            return self.thermo.cv_mole()

    property cv_mass:
        """Specific heat capacity at constant volume [J/kg/K]."""
        def __get__(self):
            return self.thermo.cv_mass()

    property cp:
        """
        Heat capacity at constant pressure [J/kg/K or J/kmol/K] depending
        on `basis`.
        """
        def __get__(self):
            return self.thermo.cp_mole() * self._mole_factor()

    property cp_mole:
        """Molar heat capacity at constant pressure [J/kmol/K]."""
        def __get__(self):
            return self.thermo.cp_mole()

    property cp_mass:
        """Specific heat capacity at constant pressure [J/kg/K]."""
        def __get__(self):
            return self.thermo.cp_mass()

    property critical_temperature:
        """Critical temperature [K]."""
        def __get__(self):
            return self.thermo.critTemperature()

    property critical_pressure:
        """Critical pressure [Pa]."""
        def __get__(self):
            return self.thermo.critPressure()

    property critical_density:
        """Critical density [kg/m^3 or kmol/m^3] depending on `basis`."""
        def __get__(self):
            return self.thermo.critDensity() / self._mass_factor()

    property P_sat:
        """Saturation pressure [Pa] at the current temperature."""
        def __get__(self):
            return self.thermo.satPressure(self.T)

    property T_sat:
        """Saturation temperature [K] at the current pressure."""
        def __get__(self):
            return self.thermo.satTemperature(self.P)

    ######## Methods to get/set the complete thermodynamic state ########

    property state_size:
        """
        Return size of vector defining internal state of the phase.
        """
        def __get__(self):
            return self.thermo.stateSize()

    property state:
        """
        Get/Set the full thermodynamic state as a single array, arranged as
        [temperature, density, mass fractions] for most phases. Useful mainly
        in cases where it is desired to store many states in a multidimensional
        array.
        """
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] state = np.empty(self.state_size)
            self.thermo.saveState(len(state), &state[0])
            return state

        def __set__(self, state):
            cdef np.ndarray[np.double_t, ndim=1] cstate = np.asarray(state)
            self.thermo.restoreState(len(state), &cstate[0])

    property TD:
        """Get/Set temperature [K] and density [kg/m^3 or kmol/m^3]."""
        def __get__(self):
            return self.T, self.density
        def __set__(self, values):
            assert len(values) == 2, 'incorrect number of values'
            T = values[0] if values[0] is not None else self.T
            D = values[1] if values[1] is not None else self.density
            self.thermo.setState_TR(T, D * self._mass_factor())

    property TDX:
        """
        Get/Set temperature [K], density [kg/m^3 or kmol/m^3], and mole
        fractions.
        """
        def __get__(self):
            return self.T, self.density, self.X
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            T = values[0] if values[0] is not None else self.T
            D = values[1] if values[1] is not None else self.density
            self.X = values[2]
            self.thermo.setState_TR(T, D * self._mass_factor())

    property TDY:
        """
        Get/Set temperature [K] and density [kg/m^3 or kmol/m^3], and mass
        fractions.
        """
        def __get__(self):
            return self.T, self.density, self.Y
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            T = values[0] if values[0] is not None else self.T
            D = values[1] if values[1] is not None else self.density
            self.Y = values[2]
            self.thermo.setState_TR(T, D * self._mass_factor())

    property TP:
        """Get/Set temperature [K] and pressure [Pa]."""
        def __get__(self):
            return self.T, self.P
        def __set__(self, values):
            assert len(values) == 2, 'incorrect number of values'
            T = values[0] if values[0] is not None else self.T
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_TP(T, P)

    property TPX:
        """Get/Set temperature [K], pressure [Pa], and mole fractions."""
        def __get__(self):
            return self.T, self.P, self.X
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            T = values[0] if values[0] is not None else self.T
            P = values[1] if values[1] is not None else self.P
            self.X = values[2]
            self.thermo.setState_TP(T, P)

    property TPY:
        """Get/Set temperature [K], pressure [Pa], and mass fractions."""
        def __get__(self):
            return self.T, self.P, self.Y
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            T = values[0] if values[0] is not None else self.T
            P = values[1] if values[1] is not None else self.P
            self.Y = values[2]
            self.thermo.setState_TP(T, P)

    property UV:
        """
        Get/Set internal energy [J/kg or J/kmol] and specific volume
        [m^3/kg or m^3/kmol].
        """
        def __get__(self):
            return self.u, self.v
        def __set__(self, values):
            assert len(values) == 2, 'incorrect number of values'
            U = values[0] if values[0] is not None else self.u
            V = values[1] if values[1] is not None else self.v
            self.thermo.setState_UV(U / self._mass_factor(),
                                    V / self._mass_factor())

    property UVX:
        """
        Get/Set internal energy [J/kg or J/kmol], specific volume
        [m^3/kg or m^3/kmol], and mole fractions.
        """
        def __get__(self):
            return self.u, self.v, self.X
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            U = values[0] if values[0] is not None else self.u
            V = values[1] if values[1] is not None else self.v
            self.X = values[2]
            self.thermo.setState_UV(U / self._mass_factor(),
                                    V / self._mass_factor())

    property UVY:
        """
        Get/Set internal energy [J/kg or J/kmol], specific volume
        [m^3/kg or m^3/kmol], and mass fractions.
        """
        def __get__(self):
            return self.u, self.v, self.Y
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            U = values[0] if values[0] is not None else self.u
            V = values[1] if values[1] is not None else self.v
            self.Y = values[2]
            self.thermo.setState_UV(U / self._mass_factor(),
                                    V / self._mass_factor())

    property DP:
        """Get/Set density [kg/m^3] and pressure [Pa]."""
        def __get__(self):
            return self.density, self.P
        def __set__(self, values):
            assert len(values) == 2, 'incorrect number of values'
            D = values[0] if values[0] is not None else self.density
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_RP(D*self._mass_factor(), P)

    property DPX:
        """Get/Set density [kg/m^3], pressure [Pa], and mole fractions."""
        def __get__(self):
            return self.density, self.P, self.X
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            D = values[0] if values[0] is not None else self.density
            P = values[1] if values[1] is not None else self.P
            self.X = values[2]
            self.thermo.setState_RP(D*self._mass_factor(), P)

    property DPY:
        """Get/Set density [kg/m^3], pressure [Pa], and mass fractions."""
        def __get__(self):
            return self.density, self.P, self.Y
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            D = values[0] if values[0] is not None else self.density
            P = values[1] if values[1] is not None else self.P
            self.Y = values[2]
            self.thermo.setState_RP(D*self._mass_factor(), P)

    property HP:
        """Get/Set enthalpy [J/kg or J/kmol] and pressure [Pa]."""
        def __get__(self):
            return self.h, self.P
        def __set__(self, values):
            assert len(values) == 2, 'incorrect number of values'
            H = values[0] if values[0] is not None else self.h
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_HP(H / self._mass_factor(), P)

    property HPX:
        """Get/Set enthalpy [J/kg or J/kmol], pressure [Pa] and mole fractions."""
        def __get__(self):
            return self.h, self.P, self.X
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            H = values[0] if values[0] is not None else self.h
            P = values[1] if values[1] is not None else self.P
            self.X = values[2]
            self.thermo.setState_HP(H / self._mass_factor(), P)

    property HPY:
        """Get/Set enthalpy [J/kg or J/kmol], pressure [Pa] and mass fractions."""
        def __get__(self):
            return self.h, self.P, self.Y
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            H = values[0] if values[0] is not None else self.h
            P = values[1] if values[1] is not None else self.P
            self.Y = values[2]
            self.thermo.setState_HP(H / self._mass_factor(), P)

    property SP:
        """Get/Set entropy [J/kg/K or J/kmol/K] and pressure [Pa]."""
        def __get__(self):
            return self.s, self.P
        def __set__(self, values):
            assert len(values) == 2, 'incorrect number of values'
            S = values[0] if values[0] is not None else self.s
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_SP(S / self._mass_factor(), P)

    property SPX:
        """Get/Set entropy [J/kg/K or J/kmol/K], pressure [Pa], and mole fractions."""
        def __get__(self):
            return self.s, self.P, self.X
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            S = values[0] if values[0] is not None else self.s
            P = values[1] if values[1] is not None else self.P
            self.X = values[2]
            self.thermo.setState_SP(S / self._mass_factor(), P)

    property SPY:
        """Get/Set entropy [J/kg/K or J/kmol/K], pressure [Pa], and mass fractions."""
        def __get__(self):
            return self.s, self.P, self.Y
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            S = values[0] if values[0] is not None else self.s
            P = values[1] if values[1] is not None else self.P
            self.Y = values[2]
            self.thermo.setState_SP(S / self._mass_factor(), P)

    property SV:
        """
        Get/Set entropy [J/kg/K or J/kmol/K] and specific volume [m^3/kg or
        m^3/kmol].
        """
        def __get__(self):
            return self.s, self.v
        def __set__(self, values):
            assert len(values) == 2, 'incorrect number of values'
            S = values[0] if values[0] is not None else self.s
            V = values[1] if values[1] is not None else self.v
            self.thermo.setState_SV(S / self._mass_factor(),
                                    V / self._mass_factor())

    property SVX:
        """
        Get/Set entropy [J/kg/K or J/kmol/K], specific volume [m^3/kg or
        m^3/kmol], and mole fractions.
        """
        def __get__(self):
            return self.s, self.v, self.X
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            S = values[0] if values[0] is not None else self.s
            V = values[1] if values[1] is not None else self.v
            self.X = values[2]
            self.thermo.setState_SV(S / self._mass_factor(),
                                    V / self._mass_factor())

    property SVY:
        """
        Get/Set entropy [J/kg/K or J/kmol/K], specific volume [m^3/kg or
        m^3/kmol], and mass fractions.
        """
        def __get__(self):
            return self.s, self.v, self.Y
        def __set__(self, values):
            assert len(values) == 3, 'incorrect number of values'
            S = values[0] if values[0] is not None else self.s
            V = values[1] if values[1] is not None else self.v
            self.Y = values[2]
            self.thermo.setState_SV(S / self._mass_factor(),
                                    V / self._mass_factor())

    # partial molar / non-dimensional properties
    property partial_molar_enthalpies:
        """Array of species partial molar enthalpies [J/kmol]."""
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarEnthalpies)

    property partial_molar_entropies:
        """Array of species partial molar entropies [J/kmol/K]."""
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarEntropies)

    property partial_molar_int_energies:
        """Array of species partial molar internal energies [J/kmol]."""
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarIntEnergies)

    property chemical_potentials:
        """Array of species chemical potentials [J/kmol]."""
        def __get__(self):
            return self._getArray1(thermo_getChemPotentials)

    property electrochemical_potentials:
        """Array of species electrochemical potentials [J/kmol]."""
        def __get__(self):
            return self._getArray1(thermo_getElectrochemPotentials)

    property partial_molar_cp:
        """
        Array of species partial molar specific heat capacities at constant
        pressure [J/kmol/K].
        """
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarCp)

    property partial_molar_volumes:
        """Array of species partial molar volumes [m^3/kmol]."""
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarVolumes)

    property standard_enthalpies_RT:
        """
        Array of nondimensional species standard-state enthalpies at the
        current temperature and pressure.
        """
        def __get__(self):
            return self._getArray1(thermo_getEnthalpy_RT)

    property standard_entropies_R:
        """
        Array of nondimensional species standard-state entropies at the
        current temperature and pressure.
        """
        def __get__(self):
            return self._getArray1(thermo_getEntropy_R)

    property standard_int_energies_RT:
        """
        Array of nondimensional species standard-state internal energies at the
        current temperature and pressure.
        """
        def __get__(self):
            return self._getArray1(thermo_getIntEnergy_RT)

    property standard_gibbs_RT:
        """
        Array of nondimensional species standard-state Gibbs free energies at
        the current temperature and pressure.
        """
        def __get__(self):
            return self._getArray1(thermo_getGibbs_RT)

    property standard_cp_R:
        """
        Array of nondimensional species standard-state specific heat capacities
        at constant pressure at the current temperature and pressure.
        """
        def __get__(self):
            return self._getArray1(thermo_getCp_R)

    property activities:
        """
        Array of nondimensional activities. Returns either molar or molal
        activities depending on the convention of the thermodynamic model.
        """
        def __get__(self):
            return self._getArray1(thermo_getActivities)

    property activity_coefficients:
        """
        Array of nondimensional, molar activity coefficients.
        """
        def __get__(self):
            return self._getArray1(thermo_getActivityCoefficients)

    ######## Miscellaneous properties ########
    property isothermal_compressibility:
        """Isothermal compressibility [1/Pa]."""
        def __get__(self):
            return self.thermo.isothermalCompressibility()

    property thermal_expansion_coeff:
        """Thermal expansion coefficient [1/K]."""
        def __get__(self):
            return self.thermo.thermalExpansionCoeff()

    property min_temp:
        """
        Minimum temperature for which the thermodynamic data for the phase are
        valid.
        """
        def __get__(self):
            return self.thermo.minTemp()

    property max_temp:
        """
        Maximum temperature for which the thermodynamic data for the phase are
        valid.
        """
        def __get__(self):
            return self.thermo.maxTemp()

    property reference_pressure:
        """Reference state pressure [Pa]."""
        def __get__(self):
            return self.thermo.refPressure()

    property electric_potential:
        """Get/Set the electric potential [V] for this phase."""
        def __get__(self):
            return self.thermo.electricPotential()
        def __set__(self, double value):
            self.thermo.setElectricPotential(value)

    property standard_concentration_units:
        """Get standard concentration units for this phase."""
        def __get__(self):
            cdef CxxUnits units = self.thermo.standardConcentrationUnits()
            return Units.copy(units)

    # methods for plasma
    property Te:
        """Get/Set electron Temperature [K]."""
        def __get__(self):
                return self.thermo.electronTemperature()

        def __set__(self, value):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            self.plasma.setElectronTemperature(value)

    def set_discretized_electron_energy_distribution(self, levels, distribution):
        """
        Set electron energy distribution. When this method is used, electron
        temperature is calculated from the distribution.

        :param levels:
            vector of electron energy levels [eV]
        :param distribution:
            vector of distribution
        """
        if not self._enable_plasma:
            raise TypeError('This method is invalid for '
                            f'thermo model: {self.thermo_model}.')
        # check length
        if (len(levels) != len(distribution)):
            raise ValueError('Length of levels and distribution are different')

        cdef np.ndarray[np.double_t, ndim=1] data_levels = \
            np.ascontiguousarray(levels, dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] data_dist = \
            np.ascontiguousarray(distribution, dtype=np.double)

        self.plasma.setDiscretizedElectronEnergyDist(&data_levels[0],
                                                     &data_dist[0],
                                                     len(levels))

    property n_electron_energy_levels:
        """ Number of electron energy levels """
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            return self.plasma.nElectronEnergyLevels()

    property electron_energy_levels:
        """ Electron energy levels [eV]"""
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(
                self.n_electron_energy_levels)
            self.plasma.getElectronEnergyLevels(&data[0])
            return data
        def __set__(self, levels):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(levels, dtype=np.double)
            self.plasma.setElectronEnergyLevels(&data[0], len(levels))

    property electron_energy_distribution:
        """ Electron energy distribution """
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(
                self.n_electron_energy_levels)
            self.plasma.getElectronEnergyDistribution(&data[0])
            return data

    property isotropic_shape_factor:
        """ Shape factor of isotropic-velocity distribution for electron energy """
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            return self.plasma.isotropicShapeFactor()
        def __set__(self, x):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            self.plasma.setIsotropicShapeFactor(x)

    property electron_energy_distribution_type:
        """ Electron energy distribution type """
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            return pystr(self.plasma.electronEnergyDistributionType())
        def __set__(self, distribution_type):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            self.plasma.setElectronEnergyDistributionType(distribution_type)

    property mean_electron_energy:
        """ Mean electron energy [eV] """
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            return self.plasma.meanElectronEnergy()
        def __set__(self, double energy):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            self.plasma.setMeanElectronEnergy(energy)

    property quadrature_method:
        """ Quadrature method """
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            return pystr(self.plasma.quadratureMethod())
        def __set__(self, method):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            self.plasma.setQuadratureMethod(stringify(method))

    property normalize_electron_energy_distribution_enabled:
        """ Automatically normalize electron energy distribuion """
        def __get__(self):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            return self.plasma.normalizeElectronEnergyDistEnabled()
        def __set__(self, enable):
            if not self._enable_plasma:
                raise ThermoModelMethodError(self.thermo_model)
            self.plasma.enableNormalizeElectronEnergyDist(enable)


cdef class InterfacePhase(ThermoPhase):
    """ A class representing a surface or edge phase"""
    def __cinit__(self, *args, **kwargs):
        if not kwargs.get("init", True):
            return
        if pystr(self.thermo.type()) not in ("Surf", "Edge"):
            raise TypeError('Underlying ThermoPhase object is of the wrong type.')
        self.surf = <CxxSurfPhase*>(self.thermo)

    property adjacent:
        """
        A dictionary containing higher-dimensional phases adjacent to this interface,
        for example bulk phases adjacent to a surface.
        """
        def __get__(self):
            return self._adjacent

    property site_density:
        """
        Get/Set the site density. [kmol/m^2] for surface phases; [kmol/m] for
        edge phases.
        """
        def __get__(self):
            return self.surf.siteDensity()
        def __set__(self, double value):
            self.surf.setSiteDensity(value)

    property coverages:
        """Get/Set the fraction of sites covered by each species."""
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.n_species)
            self.surf.getCoverages(&data[0])
            if self._selected_species.size:
                return data[self._selected_species]
            else:
                return data

        def __set__(self, theta):
            if isinstance(theta, (dict, str, bytes)):
                self.surf.setCoveragesByName(comp_map(theta))
                return

            if len(theta) != self.n_species:
                raise ValueError("Array has incorrect length."
                    " Got {}, expected {}".format(len(theta), self.n_species))
            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(theta, dtype=np.double)
            self.surf.setCoverages(&data[0])

    def set_unnormalized_coverages(self, cov):
        """
        Set the surface coverages without normalizing to force sum(cov) == 1.0.
        Useful primarily when calculating derivatives with respect to cov[k] by
        finite difference.
        """
        cdef np.ndarray[np.double_t, ndim=1] data
        if len(cov) == self.n_species:
            data = np.ascontiguousarray(cov, dtype=np.double)
        else:
            raise ValueError("Array has incorrect length."
                 " Got {}, expected {}.".format(len(cov), self.n_species))
        self.surf.setCoveragesNoNorm(&data[0])


cdef class PureFluid(ThermoPhase):
    """
    A pure substance that can  be a gas, a liquid, a mixed gas-liquid fluid,
    or a fluid beyond its critical point.
    """

    property Q:
        """
        Get/Set vapor fraction (quality). Can be set only when in the two-phase
        region.
        """
        def __get__(self):
            return self.thermo.vaporFraction()
        def __set__(self, Q):
            if (self.P >= self.critical_pressure or
                abs(self.P-self.P_sat)/self.P > 1e-4):
                raise ValueError('Cannot set vapor quality outside the'
                                 'two-phase region')
            self.thermo.setState_Psat(self.P, Q)

    property TQ:
        """Get/Set the temperature [K] and vapor fraction of a two-phase state."""
        def __get__(self):
            return self.T, self.Q
        def __set__(self, values):
            T = values[0] if values[0] is not None else self.T
            Q = values[1] if values[1] is not None else self.Q
            self.thermo.setState_Tsat(T, Q)

    property PQ:
        """Get/Set the pressure [Pa] and vapor fraction of a two-phase state."""
        def __get__(self):
            return self.P, self.Q
        def __set__(self, values):
            P = values[0] if values[0] is not None else self.P
            Q = values[1] if values[1] is not None else self.Q
            self.thermo.setState_Psat(P, Q)

    property ST:
        """Get/Set the entropy [J/kg/K] and temperature [K] of a PureFluid."""
        def __get__(self):
            return self.s, self.T
        def __set__(self, values):
            S = values[0] if values[0] is not None else self.s
            T = values[1] if values[1] is not None else self.T
            self.thermo.setState_ST(S / self._mass_factor(), T)

    property TV:
        """
        Get/Set the temperature [K] and specific volume [m^3/kg] of
        a PureFluid.
        """
        def __get__(self):
            return self.T, self.v
        def __set__(self, values):
            T = values[0] if values[0] is not None else self.T
            V = values[1] if values[1] is not None else self.v
            self.thermo.setState_TV(T, V / self._mass_factor())

    property PV:
        """
        Get/Set the pressure [Pa] and specific volume [m^3/kg] of
        a PureFluid.
        """
        def __get__(self):
            return self.P, self.v
        def __set__(self, values):
            P = values[0] if values[0] is not None else self.P
            V = values[1] if values[1] is not None else self.v
            self.thermo.setState_PV(P, V / self._mass_factor())

    property UP:
        """
        Get/Set the specific internal energy [J/kg] and the
        pressure [Pa] of a PureFluid.
        """
        def __get__(self):
            return self.u, self.P
        def __set__(self, values):
            U = values[0] if values[0] is not None else self.u
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_UP(U / self._mass_factor(), P)

    property VH:
        """
        Get/Set the specific volume [m^3/kg] and the specific
        enthalpy [J/kg] of a PureFluid.
        """
        def __get__(self):
            return self.v, self.h
        def __set__(self, values):
            V = values[0] if values[0] is not None else self.v
            H = values[1] if values[1] is not None else self.h
            self.thermo.setState_VH(V/self._mass_factor(), H/self._mass_factor())

    property TH:
        """
        Get/Set the temperature [K] and the specific enthalpy [J/kg]
        of a PureFluid.
        """
        def __get__(self):
            return self.T, self.h
        def __set__(self, values):
            T = values[0] if values[0] is not None else self.T
            H = values[1] if values[1] is not None else self.h
            self.thermo.setState_TH(T, H / self._mass_factor())

    property SH:
        """
        Get/Set the specific entropy [J/kg/K] and the specific
        enthalpy [J/kg] of a PureFluid.
        """
        def __get__(self):
            return self.s, self.h
        def __set__(self, values):
            S = values[0] if values[0] is not None else self.s
            H = values[1] if values[1] is not None else self.h
            self.thermo.setState_SH(S/self._mass_factor(), H/self._mass_factor())

    property TDQ:
        """
        Get the temperature [K], density [kg/m^3 or kmol/m^3], and vapor
        fraction.
        """
        def __get__(self):
            return self.T, self.density, self.Q

    property TPQ:
        """
        Get/Set the temperature [K], pressure [Pa], and vapor fraction of a
        PureFluid.

        An Exception is raised if the thermodynamic state is not consistent.
        """
        def __get__(self):
            return self.T, self.P, self.Q
        def __set__(self, values):
            T = values[0] if values[0] is not None else self.T
            P = values[1] if values[1] is not None else self.P
            Q = values[2] if values[2] is not None else self.Q
            self.thermo.setState_TPQ(T, P, Q)

    property UVQ:
        """
        Get the internal energy [J/kg or J/kmol], specific volume
        [m^3/kg or m^3/kmol], and vapor fraction.
        """
        def __get__(self):
            return self.u, self.v, self.Q

    property DPQ:
        """Get the density [kg/m^3], pressure [Pa], and vapor fraction."""
        def __get__(self):
            return self.density, self.P, self.Q

    property HPQ:
        """
        Get the enthalpy [J/kg or J/kmol], pressure [Pa] and vapor fraction.
        """
        def __get__(self):
            return self.h, self.P, self.Q

    property SPQ:
        """
        Get the entropy [J/kg/K or J/kmol/K], pressure [Pa], and vapor fraction.
        """
        def __get__(self):
            return self.s, self.P, self.Q

    property SVQ:
        """
        Get the entropy [J/kg/K or J/kmol/K], specific volume [m^3/kg or
        m^3/kmol], and vapor fraction.
        """
        def __get__(self):
            return self.s, self.v, self.Q


class Element:
    """
    An element or a named isotope defined in Cantera.

    Class `Element` gets data for the elements and isotopes defined in
    :ct:`Elements.cpp`. This class can be used in two ways. The
    first way is to get information about all of the elements stored in
    Cantera. The three attributes `num_elements_defined`,
    `element_symbols`, and `element_names` can be accessed by::

        >>> ct.Element.num_elements_defined
        >>> ct.Element.element_symbols
        >>> ct.Element.element_names

    Otherwise, if the class `Element` is called with an argument, it
    stores the data about that particular element. For example::

        >>> ar_sym = ct.Element('Ar')
        >>> ar_name = ct.Element('argon')
        >>> ar_num = ct.Element(18)

    would all create instances with the information for argon. The
    available argument options to create an instance of the `Element`
    class with the element information are the `name`, `symbol`, and
    `atomic_number`. Once an instance of the class is made, the `name`,
    `atomic_number`, `symbol`, and atomic `weight` can be accessed as
    attributes of the instance of the `Element` class.

        >>> ar_sym.name
        'argon'
        >>> ar_sym.weight
        39.948
        >>> ar_sym.atomic_number
        18
        >>> ar_sym.symbol
        'Ar'

    The elements available are listed below, in the `element_symbols`
    and `element_names` attribute documentation.
    """

    #: The number of named elements (not isotopes) defined in Cantera
    num_elements_defined = numElementsDefined()

    #: A list of the symbols of all the elements (not isotopes) defined
    #: in Cantera
    element_symbols = [pystr(getElementSymbol(<int>(m+1)))
                       for m in range(num_elements_defined)]

    #: A list of the names of all the elements (not isotopes) defined
    #: in Cantera
    element_names = [pystr(getElementName(<int>m+1))
                     for m in range(num_elements_defined)]

    def __init__(self, arg):
        if isinstance(arg, (str, bytes)):
            try:
                # Assume the argument is the element symbol and try to get the name
                self._name = pystr(getElementName(stringify(arg)))
            except CanteraError:
                # If getting the name failed, the argument must be the name
                self._symbol = pystr(getElementSymbol(stringify(arg)))
                self._name = arg.lower()
            else:
                self._symbol = arg

            self._atomic_number = getAtomicNumber(stringify(arg))
            self._weight = getElementWeight(stringify(arg))
        elif isinstance(arg, int):
            self._atomic_number = arg
            self._name = pystr(getElementName(<int>arg))
            self._symbol = pystr(getElementSymbol(<int>arg))
            self._weight = getElementWeight(<int>arg)
        else:
            raise TypeError('The input argument to Element must be a string '
                            'or an integer')

    @property
    def name(self):
        """The name of the element or isotope."""
        return self._name

    @property
    def atomic_number(self):
        """The atomic number of the element or isotope."""
        return self._atomic_number

    @property
    def symbol(self):
        """The symbol of the element or isotope."""
        return self._symbol

    @property
    def weight(self):
        """The atomic weight of the element or isotope."""
        return self._weight
