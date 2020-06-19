# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import warnings
import weakref
import numbers as _numbers

cdef enum ThermoBasis:
    mass_basis = 0
    molar_basis = 1


cdef class Species:
    """
    A class which stores data about a single chemical species that may be
    needed to add it to a `Solution` or `Interface` object (and to the
    underlying `ThermoPhase` and `Transport` objects).

    :param name:
        A string giving the name of the species, e.g. ``'CH4'``
    :param composition:
        The elemental composition of the species, given either as a dict or a
        composition string, e.g. ``{'C':1, 'H':4}`` or ``'C:1, H:4'``.
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

        S = ct.Species.listFromFile('gri30.yaml')

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
        """
        cxx_species = CxxNewSpecies(AnyMapFromYamlString(stringify(text)))
        species = Species(init=False)
        species._assign(cxx_species)
        return species

    @staticmethod
    def listFromFile(filename, section='species'):
        """
        Create a list of Species objects from all of the species defined in a
        YAML, CTI or XML file. For YAML files, return species from the section
        *section*.

        Directories on Cantera's input file path will be searched for the
        specified file.

        In the case of an XML file, the ``<species>`` nodes are assumed to be
        children of the ``<speciesData>`` node in a document with a ``<ctml>``
        root node, as in the XML files produced by conversion from CTI files.

        .. deprecated:: 2.5

            The CTI and XML input formats are deprecated and will be removed in
            Cantera 3.0.
        """
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

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if 'source' not in kwargs:
            self.thermo_basis = mass_basis
        self._references = weakref.WeakKeyDictionary()

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

    property ID:
        """
        The identifier of the object. The default value corresponds to the
        CTI/XML/YAML input file phase entry.

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Usage merged with `name`.
        """
        def __get__(self):
            warnings.warn("To be removed after Cantera 2.5. "
                          "Use 'name' attribute instead", DeprecationWarning)
            return pystr(self.base.name())
        def __set__(self, id_):
            warnings.warn("To be removed after Cantera 2.5. "
                          "Use 'name' attribute instead", DeprecationWarning)
            self.base.setName(stringify(id_))

    property basis:
        """
        Determines whether intensive thermodynamic properties are treated on a
        `mass` (per kg) or `molar` (per kmol) basis. This affects the values
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
                    int log_level=0, **kwargs):
        """
        Set to a state of chemical equilibrium holding property pair
        *XY* constant.

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
        if 'maxsteps' in kwargs:
            max_steps = kwargs['maxsteps']
            warnings.warn(
                "Keyword argument 'maxsteps' is deprecated and will be removed after "
                "Cantera 2.5. Use argument 'max_steps' instead.", DeprecationWarning,
            )

        if 'maxiter' in kwargs:
            max_iter = kwargs['maxiter']
            warnings.warn(
                "Keyword argument 'maxiter' is deprecated and will be removed after "
                "Cantera 2.5. Use argument 'max_iter' instead.", DeprecationWarning,
            )

        if 'loglevel' in kwargs:
            log_level = kwargs['loglevel']
            warnings.warn(
                "Keyword argument 'loglevel' is deprecated and will be removed after "
                "Cantera 2.5. Use argument 'log_level' instead.", DeprecationWarning,
            )

        self.thermo.equilibrate(stringify(XY.upper()), stringify(solver), rtol,
                                max_steps, max_iter, estimate_equil, log_level)

    ####### Composition, species, and elements ########

    property n_elements:
        """Number of elements."""
        def __get__(self):
            return self.thermo.nElements()

    cpdef int element_index(self, element) except *:
        """
        The index of element *element*, which may be specified as a string or
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
        """Name of the element with index *m*."""
        return pystr(self.thermo.elementName(m))

    property element_names:
        """A list of all the element names."""
        def __get__(self):
            return [self.element_name(m) for m in range(self.n_elements)]

    def atomic_weight(self, m):
        """Atomic weight [kg/kmol] of element *m*"""
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
        """Name of the species with index *k*."""
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
        The index of species *species*, which may be specified as a string or
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
        Return the `Species` object for species *k*, where *k* is either the
        species index or the species name. If *k* is not specified, a list of
        all species objects is returned. Changes to this object do not affect
        the `ThermoPhase` or `Solution` object until the `modify_species`
        function is called.
        """
        if k is None:
            return [self.species(i) for i in range(self.n_species)]

        s = Species(init=False)
        if isinstance(k, (str, bytes)):
            s._assign(self.thermo.species(stringify(k)))
        elif isinstance(k, (int, float)):
            s._assign(self.thermo.species(<int>k))
        else:
            raise TypeError("Argument must be a string or a number."
                            " Got {!r}.".format(k))
        return s

    def modify_species(self, k, Species species):
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
        self.thermo.initThermo()
        if self.kinetics:
            self.kinetics.invalidateCache()

    def add_species_alias(self, name, alias):
        """
        Add the alternate species name *alias* for an original species *name*.
        """
        self.thermo.addSpeciesAlias(stringify(name), stringify(alias))

    def find_isomers(self, comp):
        """
        Find species/isomers matching a composition specified by *comp*.
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
        Number of atoms of element *element* in species *species*. The element
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

    def set_equivalence_ratio(self, phi, fuel, oxidizer):
        """
        Set the composition to a mixture of *fuel* and *oxidizer* at the
        specified equivalence ratio *phi*, holding temperature and pressure
        constant. Considers the oxidation of C to CO2, H to H2O and S to SO2.
        Other elements are assumed not to participate in oxidation (i.e. N ends up as
        N2)::

            >>> gas.set_equivalence_ratio(0.5, 'CH4', 'O2:1.0, N2:3.76')
            >>> gas.mole_fraction_dict()
            {'CH4': 0.049900199, 'N2': 0.750499001, 'O2': 0.199600798}

            >>> gas.set_equivalence_ratio(1.2, {'NH3':0.8, 'CO':0.2}, 'O2:1.0')
            >>> gas.mole_fraction_dict()
            {'CO': 0.1263157894, 'NH3': 0.505263157, 'O2': 0.36842105}

        :param phi: Equivalence ratio
        :param fuel:
            Fuel species name or molar composition as string, array, or dict.
        :param oxidizer:
            Oxidizer species name or molar composition as a string, array, or
            dict.
        """
        if (isinstance(fuel, str) and ':' not in fuel
            and fuel in self.species_names):
            fuel += ':1.0'

        if (isinstance(oxidizer, str) and ':' not in oxidizer
            and oxidizer in self.species_names):
            oxidizer += ':1.0'

        self.TPX = None, None, fuel
        Xf = self.X
        self.TPX = None, None, oxidizer
        Xo = self.X

        nO = np.array([self.n_atoms(k, 'O') for k in range(self.n_species)])

        if 'C' in self.element_names:
            nC = np.array([self.n_atoms(k, 'C') for k in range(self.n_species)])
        else:
            nC = np.zeros(self.n_species)

        if 'H' in self.element_names:
            nH = np.array([self.n_atoms(k, 'H') for k in range(self.n_species)])
        else:
            nH = np.zeros(self.n_species)

        if 'S' in self.element_names:
            nS = np.array([self.n_atoms(k, 'S') for k in range(self.n_species)])
        else:
            nS = np.zeros(self.n_species)

        Cf = nC.dot(Xf)
        Co = nC.dot(Xo)
        Of = nO.dot(Xf)
        Oo = nO.dot(Xo)
        Hf = nH.dot(Xf)
        Ho = nH.dot(Xo)
        Sf = nS.dot(Xf)
        So = nS.dot(Xo)

        stoichAirFuelRatio = - (Of - 2*Cf - 2*Sf - Hf/2.0) / (Oo - 2*Co - 2*So - Ho/2.0)
        Xr = phi * Xf + stoichAirFuelRatio * Xo
        self.TPX = None, None, Xr

    def get_equivalence_ratio(self, oxidizers=[], ignore=[]):
        """
        Get the composition of a fuel/oxidizer mixture. This gives the
        equivalence ratio of an unburned mixture. This is not a quantity that is
        conserved after oxidation. Considers the oxidation of C to CO2, H to H2O
        and S to SO2. Other elements are assumed not to participate in oxidation
        (i.e. N ends up as N2).

        :param oxidizers:
            List of oxidizer species names as strings. Default: with
            ``oxidizers=[]``, every species that contains O but does not contain
            H, C, or S is considered to be an oxidizer.
        :param ignore:
            List of species names as strings to ignore.

        >>> gas.set_equivalence_ratio(0.5, 'CH3:0.5, CH3OH:.5, N2:0.125', 'O2:0.21, N2:0.79, NO:0.01')
        >>> gas.get_equivalence_ratio()
        0.50000000000000011
        >>> gas.get_equivalence_ratio(['O2'])  # Only consider O2 as the oxidizer instead of O2 and NO
        0.48809523809523814
        >>> gas.X = 'CH4:1, O2:2, NO:0.1'
        >>> gas.get_equivalence_ratio(ignore=['NO'])
        1.0
        """
        if not oxidizers:
            # Default behavior, find all possible oxidizers
            oxidizers = []
            for s in self.species():
                if all(y not in s.composition for y in ['C', 'H', 'S']):
                    oxidizers.append(s.name)

        alpha = 0
        mol_O = 0
        for k, s in enumerate(self.species()):
            if s.name in ignore:
                continue
            elif s.name in oxidizers:
                mol_O += s.composition.get('O', 0) * self.X[k]
            else:
                nC = s.composition.get('C', 0)
                nH = s.composition.get('H', 0)
                nO = s.composition.get('O', 0)
                nS = s.composition.get('S', 0)

                alpha += (2 * nC + nH / 2 + 2 * nS - nO) * self.X[k]

        if mol_O == 0:
            return float('inf')
        else:
            return alpha / mol_O

    def elemental_mass_fraction(self, m):
        r"""
        Get the elemental mass fraction :math:`Z_{\mathrm{mass},m}` of element
        :math:`m` as defined by:

        .. math:: Z_{\mathrm{mass},m} = \sum_k \frac{a_{m,k} M_m}{M_k} Y_k

        with :math:`a_{m,k}` being the number of atoms of element :math:`m` in
        species :math:`k`, :math:`M_m` the atomic weight of element :math:`m`,
        :math:`M_k` the molecular weight of species :math:`k`, and :math:`Y_k`
        the mass fraction of species :math:`k`.

        :param m:
            Base element, may be specified by name or by index.

        >>> phase.elemental_mass_fraction('H')
        1.0
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
        :math:`X_k` being the mole fraction of species :math:`k`.

        :param m:
            Base element, may be specified by name or by index.

        >>> phase.elemental_mole_fraction('H')
        1.0
        """
        return self.thermo.elementalMoleFraction(self.element_index(m))

    def set_unnormalized_mass_fractions(self, Y):
        """
        Set the mass fractions without normalizing to force sum(Y) == 1.0.
        Useful primarily when calculating derivatives with respect to Y[k] by
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
        Set the mole fractions without normalizing to force sum(X) == 1.0.
        Useful primarily when calculating derivatives with respect to X[k]
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
        Y = self.thermo.getMassFractionsByName(threshold)
        return {pystr(item.first):item.second for item in Y}

    def mole_fraction_dict(self, double threshold=0.0):
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
        """Specific entropy [J/kg]."""
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


cdef class InterfacePhase(ThermoPhase):
    """ A class representing a surface or edge phase"""
    def __cinit__(self, *args, **kwargs):
        if pystr(self.thermo.type()) not in ("Surf", "Edge"):
            raise TypeError('Underlying ThermoPhase object is of the wrong type.')
        self.surf = <CxxSurfPhase*>(self.thermo)

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

    property X:
        """
        Get/Set vapor fraction (quality). Can be set only when in the two-phase
        region.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `Q`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'Q'", DeprecationWarning)
            return self.Q
        def __set__(self, X):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'Q'", DeprecationWarning)
            self.Q = X

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

    property TX:
        """Get/Set the temperature [K] and vapor fraction of a two-phase state.

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Renamed to `TQ`.
        """
        def __get__(self):
            warnings.warn("To be removed after Cantera 2.5. "
                          "Attribute renamed to 'TQ'", DeprecationWarning)
            return self.TQ
        def __set__(self, values):
            warnings.warn("To be removed after Cantera 2.5. "
                          "Attribute renamed to 'TQ'", DeprecationWarning)
            self.TQ = values

    property TQ:
        """Get/Set the temperature [K] and vapor fraction of a two-phase state."""
        def __get__(self):
            return self.T, self.Q
        def __set__(self, values):
            T = values[0] if values[0] is not None else self.T
            Q = values[1] if values[1] is not None else self.Q
            self.thermo.setState_Tsat(T, Q)

    property PX:
        """Get/Set the pressure [Pa] and vapor fraction of a two-phase state.

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Renamed to `PQ`.
        """
        def __get__(self):
            warnings.warn("To be removed after Cantera 2.5. "
                          "Attribute renamed to 'PQ'", DeprecationWarning)
            return self.PQ
        def __set__(self, values):
            warnings.warn("To be removed after Cantera 2.5. "
                          "Attribute renamed to 'PQ'", DeprecationWarning)
            self.PQ = values

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

    property TDX:
        """
        Get the temperature [K], density [kg/m^3 or kmol/m^3], and vapor
        fraction.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `TDQ`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'TDQ'", DeprecationWarning)
            return self.TDQ

    property TDQ:
        """
        Get the temperature [K], density [kg/m^3 or kmol/m^3], and vapor
        fraction.
        """
        def __get__(self):
            return self.T, self.density, self.Q

    property TPX:
        """
        Get/Set the temperature [K], pressure [Pa], and vapor fraction of a
        PureFluid.

        An Exception is raised if the thermodynamic state is not consistent.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `TPQ`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'TPQ'", DeprecationWarning)
            return self.TPQ
        def __set__(self, values):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'TPQ'", DeprecationWarning)
            self.TPQ = values

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

    property UVX:
        """
        Get the internal energy [J/kg or J/kmol], specific volume
        [m^3/kg or m^3/kmol], and vapor fraction.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `UVQ`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'UVQ'", DeprecationWarning)
            return self.UVQ

    property UVQ:
        """
        Get the internal energy [J/kg or J/kmol], specific volume
        [m^3/kg or m^3/kmol], and vapor fraction.
        """
        def __get__(self):
            return self.u, self.v, self.Q

    property DPX:
        """Get the density [kg/m^3], pressure [Pa], and vapor fraction.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `DPQ`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'DPQ'", DeprecationWarning)
            return self.DPQ

    property DPQ:
        """Get the density [kg/m^3], pressure [Pa], and vapor fraction."""
        def __get__(self):
            return self.density, self.P, self.Q

    property HPX:
        """
        Get the enthalpy [J/kg or J/kmol], pressure [Pa] and vapor fraction.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `HPQ`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'HPQ'", DeprecationWarning)
            return self.HPQ

    property HPQ:
        """
        Get the enthalpy [J/kg or J/kmol], pressure [Pa] and vapor fraction.
        """
        def __get__(self):
            return self.h, self.P, self.Q

    property SPX:
        """
        Get the entropy [J/kg/K or J/kmol/K], pressure [Pa], and vapor fraction.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `SPQ`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'SPQ'", DeprecationWarning)
            return self.SPQ

    property SPQ:
        """
        Get the entropy [J/kg/K or J/kmol/K], pressure [Pa], and vapor fraction.
        """
        def __get__(self):
            return self.s, self.P, self.Q

    property SVX:
        """
        Get the entropy [J/kg/K or J/kmol/K], specific volume [m^3/kg or
        m^3/kmol], and vapor fraction.

        .. deprecated:: 2.5

             Behavior changes after version 2.5, when `X` will refer to mole
             fraction. Renamed to `SVQ`.
        """
        def __get__(self):
            warnings.warn("Behavior changes after Cantera 2.5, "
                          "when 'X' will refer to mole fraction. "
                          "Attribute renamed to 'SVQ'", DeprecationWarning)
            return self.SVQ

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
    `src/thermo/Elements.cpp`. This class can be used in two ways. The
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
