import warnings

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

    The static methods `fromCti`, `fromXml`, `listFromFile`, `listFromCti`, and
    `listFromXml` can be used to create `Species` objects from existing
    definitions in the CTI or XML formats. All of the following will produce a
    list of 53 `Species` objects containing the species defined in the GRI 3.0
    mechanism::

        S = ct.Species.listFromFile('gri30.cti')
        S = ct.Species.listFromCti(open('path/to/gri30.cti').read())
        S = ct.Species.listFromXml(open('path/to/gri30.xml').read())

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
        """
        cxx_species = CxxNewSpecies(deref(CxxGetXmlFromString(stringify(text))))
        species = Species(init=False)
        species._assign(cxx_species)
        return species

    @staticmethod
    def listFromFile(filename):
        """
        Create a list of Species objects from all of the species defined in a
        CTI or XML file.

        Directories on Cantera's input file path will be searched for the
        specified file.

        In the case of an XML file, the ``<species>`` nodes are assumed to be
        children of the ``<speciesData>`` node in a document with a ``<ctml>``
        root node, as in the XML files produced by conversion from CTI files.
        """
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

    property name:
        """
        The name assigned to this phase. The default is taken from the CTI/XML
        input file.
        """
        def __get__(self):
            return pystr(self.thermo.name())
        def __set__(self, name):
            self.thermo.setName(stringify(name))

    property ID:
        """
        The ID of the phase. The default is taken from the CTI/XML input file.
        """
        def __get__(self):
            return pystr(self.thermo.id())
        def __set__(self, id_):
            self.thermo.setID(stringify(id_))


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
                raise ValueError("Valid choices are 'mass' or 'molar'.")

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
                    int maxsteps=1000, int maxiter=100, int estimate_equil=0,
                    int loglevel=0):
        """
        Set to a state of chemical equilibrium holding property pair
        *XY* constant.

        :param XY:
            A two-letter string, which must be one of the set::

                ['TP','TV','HP','SP','SV','UV']

        :param solver:
            Specifies the equilibrium solver to use. May be one of the following:

            * ''element_potential'' - a fast solver using the element potential
              method
            * 'gibbs' - a slower but more robust Gibbs minimization solver
            * 'vcs' - the VCS non-ideal equilibrium solver
            * "auto" - The element potential solver will be tried first, then
              if it fails the Gibbs solver will be tried.
        :param rtol:
            the relative error tolerance.
        :param maxsteps:
            maximum number of steps in composition to take to find a converged
            solution.
        :param maxiter:
            For the Gibbs minimization solver, this specifies the number of
            'outer' iterations on T or P when some property pair other
            than TP is specified.
        :param estimate_equil:
            Integer indicating whether the solver should estimate its own
            initial condition. If 0, the initial mole fraction vector in the
            ThermoPhase object is used as the initial condition. If 1, the
            initial mole fraction vector is used if the element abundances are
            satisfied. If -1, the initial mole fraction vector is thrown out,
            and an estimate is formulated.
        :param loglevel:
            Set to a value > 0 to write diagnostic output.
            """
        if isinstance(solver, int):
            warnings.warn('ThermoPhase.equilibrate: Using integer solver '
                'flags is deprecated, and will be disabled after Cantera 2.2.')
            if solver == -1:
                solver = 'auto'
            elif solver == 0:
                solver = 'element_potential'
            elif solver == 1:
                solver = 'gibbs'
            elif solver == 2:
                solver = 'vcs'
            else:
                raise ValueError('Invalid equilibrium solver specified: '
                    '"{0}"'.format(solver))

        self.thermo.equilibrate(stringify(XY.upper()), stringify(solver), rtol,
                                maxsteps, maxiter, estimate_equil, loglevel)

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
        if isinstance(element, (str, unicode, bytes)):
            index = self.thermo.elementIndex(stringify(element))
        elif isinstance(element, (int, float)):
            index = <int>element
        else:
            raise TypeError("'element' must be a string or a number")

        if not 0 <= index < self.n_elements:
            raise ValueError('No such element.')

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
        if isinstance(species, (str, unicode, bytes)):
            index = self.thermo.speciesIndex(stringify(species))
        elif isinstance(species, (int, float)):
            index = <int>species
        else:
            raise TypeError("'species' must be a string or a number")

        if not 0 <= index < self.n_species:
            raise ValueError('No such species.')

        return index

    def species(self, k=None):
        """
        Return the `Species` object for species *k*, where *k* is either the
        species index or the species name. If *k* is not specified, a list of
        all species objects is returned.
        """
        if k is None:
            return [self.species(i) for i in range(self.n_species)]

        s = Species(init=False)
        if isinstance(k, (str, unicode, bytes)):
            s._assign(self.thermo.species(stringify(k)))
        elif isinstance(k, (int, float)):
            s._assign(self.thermo.species(<int>k))
        else:
            raise TypeError("Argument must be a string or a number")
        return s

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

        if len(values) == self.n_species:
            data = np.ascontiguousarray(values, dtype=np.double)
        elif len(values) == len(self._selected_species):
            data = np.zeros(self.n_species, dtype=np.double)
            for i,k in enumerate(self._selected_species):
                data[k] = values[i]
        else:
            raise ValueError("Array has incorrect length")
        method(self.thermo, &data[0])

    property molecular_weights:
        """Array of species molecular weights (molar masses) [kg/kmol]."""
        def __get__(self):
            return self._getArray1(thermo_getMolecularWeights)

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
            if isinstance(Y, (str, unicode, bytes)):
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
            if isinstance(X, (str, unicode, bytes)):
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
        :math:`m` as defined by:

        .. math:: Z_{\mathrm{mole},m} = \sum_k \frac{a_{m,k}}{\sum_j a_{j,k}} X_k

        with :math:`a_{m,k}` being the number of atoms of element :math:`m` in
        species :math:`k` and :math:`X_k` the mole fraction of species
        :math:`k`.

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
            raise ValueError("Array has incorrect length")
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
            raise ValueError("Array has incorrect length")
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

    property TD:
        """Get/Set temperature [K] and density [kg/m^3 or kmol/m^3]."""
        def __get__(self):
            return self.T, self.density
        def __set__(self, values):
            assert len(values) == 2
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
            assert len(values) == 3
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
            assert len(values) == 3
            T = values[0] if values[0] is not None else self.T
            D = values[1] if values[1] is not None else self.density
            self.Y = values[2]
            self.thermo.setState_TR(T, D * self._mass_factor())

    property TP:
        """Get/Set temperature [K] and pressure [Pa]."""
        def __get__(self):
            return self.T, self.P
        def __set__(self, values):
            assert len(values) == 2
            T = values[0] if values[0] is not None else self.T
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_TP(T, P)

    property TPX:
        """Get/Set temperature [K], pressure [Pa], and mole fractions."""
        def __get__(self):
            return self.T, self.P, self.X
        def __set__(self, values):
            assert len(values) == 3
            T = values[0] if values[0] is not None else self.T
            P = values[1] if values[1] is not None else self.P
            self.X = values[2]
            self.thermo.setState_TP(T, P)

    property TPY:
        """Get/Set temperature [K], pressure [Pa], and mass fractions."""
        def __get__(self):
            return self.T, self.P, self.Y
        def __set__(self, values):
            assert len(values) == 3
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
            assert len(values) == 2
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
            assert len(values) == 3
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
            assert len(values) == 3
            U = values[0] if values[0] is not None else self.u
            V = values[1] if values[1] is not None else self.v
            self.Y = values[2]
            self.thermo.setState_UV(U / self._mass_factor(),
                                    V / self._mass_factor())

    property HP:
        """Get/Set enthalpy [J/kg or J/kmol] and pressure [Pa]."""
        def __get__(self):
            return self.h, self.P
        def __set__(self, values):
            assert len(values) == 2
            H = values[0] if values[0] is not None else self.h
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_HP(H / self._mass_factor(), P)

    property HPX:
        """Get/Set enthalpy [J/kg or J/kmol], pressure [Pa] and mole fractions."""
        def __get__(self):
            return self.h, self.P, self.X
        def __set__(self, values):
            assert len(values) == 3
            H = values[0] if values[0] is not None else self.h
            P = values[1] if values[1] is not None else self.P
            self.X = values[2]
            self.thermo.setState_HP(H / self._mass_factor(), P)

    property HPY:
        """Get/Set enthalpy [J/kg or J/kmol], pressure [Pa] and mass fractions."""
        def __get__(self):
            return self.h, self.P, self.Y
        def __set__(self, values):
            assert len(values) == 3
            H = values[0] if values[0] is not None else self.h
            P = values[1] if values[1] is not None else self.P
            self.Y = values[2]
            self.thermo.setState_HP(H / self._mass_factor(), P)

    property SP:
        """Get/Set entropy [J/kg/K or J/kmol/K] and pressure [Pa]."""
        def __get__(self):
            return self.s, self.P
        def __set__(self, values):
            assert len(values) == 2
            S = values[0] if values[0] is not None else self.s
            P = values[1] if values[1] is not None else self.P
            self.thermo.setState_SP(S / self._mass_factor(), P)

    property SPX:
        """Get/Set entropy [J/kg/K or J/kmol/K], pressure [Pa], and mole fractions."""
        def __get__(self):
            return self.s, self.P, self.X
        def __set__(self, values):
            assert len(values) == 3
            S = values[0] if values[0] is not None else self.s
            P = values[1] if values[1] is not None else self.P
            self.X = values[2]
            self.thermo.setState_SP(S / self._mass_factor(), P)

    property SPY:
        """Get/Set entropy [J/kg/K or J/kmol/K], pressure [Pa], and mass fractions."""
        def __get__(self):
            return self.s, self.P, self.Y
        def __set__(self, values):
            assert len(values) == 3
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
            assert len(values) == 2
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
            assert len(values) == 3
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
            assert len(values) == 3
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

    def element_potentials(self):
        """
        Get the array of element potentials. The element potentials are only
        defined for equilibrium states. This method first sets the composition
        to a state of equilibrium at constant T and P, then computes the
        element potentials for this equilibrium state.
        """
        self.equilibrate('TP')
        cdef np.ndarray[np.double_t, ndim=1] data = np.zeros(self.n_elements)
        self.thermo.getElementPotentials(&data[0])
        return data


cdef class InterfacePhase(ThermoPhase):
    """ A class representing a surface or edge phase"""
    def __cinit__(self, *args, **kwargs):
        if self.thermo.eosType() not in (thermo_type_surf, thermo_type_edge):
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
            if isinstance(theta, (dict, str, unicode, bytes)):
                self.surf.setCoveragesByName(comp_map(theta))
                return

            if len(theta) != self.n_species:
                raise ValueError("Array has incorrect length")
            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(theta, dtype=np.double)
            self.surf.setCoverages(&data[0])


cdef class PureFluid(ThermoPhase):
    """
    A pure substance that can  be a gas, a liquid, a mixed gas-liquid fluid,
    or a fluid beyond its critical point.
    """
    property X:
        """
        Get/Set vapor fraction (quality). Can be set only when in the two-phase
        region.
        """
        def __get__(self):
            return self.thermo.vaporFraction()
        def __set__(self, X):
            if (self.P >= self.critical_pressure or
                abs(self.P-self.P_sat)/self.P > 1e-4):
                raise ValueError('Cannot set vapor quality outside the'
                                 'two-phase region')
            self.thermo.setState_Psat(self.P, X)

    property TX:
        """Get/Set the temperature and vapor fraction of a two-phase state."""
        def __get__(self):
            return self.T, self.X
        def __set__(self, values):
            T = values[0] if values[0] is not None else self.T
            X = values[1] if values[1] is not None else self.X
            self.thermo.setState_Tsat(T, X)

    property PX:
        """Get/Set the pressure and vapor fraction of a two-phase state."""
        def __get__(self):
            return self.P, self.X
        def __set__(self, values):
            P = values[0] if values[0] is not None else self.P
            X = values[1] if values[1] is not None else self.X
            self.thermo.setState_Psat(P, X)
