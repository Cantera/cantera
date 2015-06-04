# NOTE: These cdef functions cannot be members of Kinetics because they would
# cause "layout conflicts" when creating derived classes with multiple bases,
# e.g. class Solution. [Cython 0.16]
cdef np.ndarray get_species_array(Kinetics kin, kineticsMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(kin.n_total_species)
    method(kin.kinetics, &data[0])
    # @TODO: Fix _selected_species to work with interface kinetics
    if kin._selected_species.size:
        return data[kin._selected_species]
    else:
        return data

cdef np.ndarray get_reaction_array(Kinetics kin, kineticsMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(kin.n_reactions)
    method(kin.kinetics, &data[0])
    return data


cdef class Kinetics(_SolutionBase):
    """
    Instances of class `Kinetics` are responsible for evaluating reaction rates
    of progress, species production rates, and other quantities pertaining to
    a reaction mechanism.
    """

    property n_total_species:
        """
        Total number of species in all phases participating in the kinetics
        mechanism.
        """
        def __get__(self):
            return self.kinetics.nTotalSpecies()

    property n_reactions:
        """Number of reactions in the reaction mechanism."""
        def __get__(self):
            return self.kinetics.nReactions()

    property n_phases:
        """Number of phases in the reaction mechanism."""
        def __get__(self):
            return self.kinetics.nPhases()

    property reaction_phase_index:
        """The index of the phase where the reactions occur."""
        def __get__(self):
            return self.kinetics.reactionPhaseIndex()

    def _check_phase_index(self, n):
        if not 0 <= n < self.n_phases:
            raise ValueError("Phase index ({0}) out of range".format(n))

    def _check_reaction_index(self, n):
        if not 0 <= n < self.n_reactions:
            raise ValueError("Reaction index ({0}) out of range".format(n))

    def _check_kinetics_species_index(self, n):
        if not 0 <= n < self.n_total_species:
            raise ValueError("Kinetics Species index ({0}) out of range".format(n))

    def kinetics_species_index(self, species, int phase=0):
        """
        The index of species *species* of phase *phase* within arrays returned
        by methods of class `Kinetics`. If *species* is a string, the *phase*
        argument is unused.
        """
        cdef int k
        if isinstance(species, (str, unicode, bytes)):
            return self.kinetics.kineticsSpeciesIndex(stringify(species))
        else:
            k = species
            self._check_kinetics_species_index(k)
            self._check_phase_index(k)
            return self.kinetics.kineticsSpeciesIndex(k, phase)

    def reaction(self, int i_reaction):
        """
        Return a `Reaction` object representing the reaction with index
        ``i_reaction``.
        """
        return wrapReaction(self.kinetics.reaction(i_reaction))

    def reactions(self):
        """
        Return a list of all `Reaction` objects
        """
        return [self.reaction(i) for i in range(self.n_reactions)]

    def modify_reaction(self, int irxn, Reaction rxn):
        """
        Modify the `Reaction` with index ``irxn`` to have the same rate
        parameters as ``rxn``. ``rxn`` must have the same reactants and products
        and be of the same type (i.e. `ElementaryReaction`, `FalloffReaction`,
        `PlogReaction`, etc.) as the existing reaction. This method does not
        modify the third-body efficiencies, reaction orders, or reversibility of
        the reaction.
        """
        self.kinetics.modifyReaction(irxn, rxn._reaction)

    def is_reversible(self, int i_reaction):
        """True if reaction `i_reaction` is reversible."""
        self._check_reaction_index(i_reaction)
        return self.kinetics.isReversible(i_reaction)

    def multiplier(self, int i_reaction):
        """
        A scaling factor applied to the rate coefficient for reaction
        *i_reaction*. Can be used to carry out sensitivity analysis or to
        selectively disable a particular reaction. See `set_multiplier`.
        """
        self._check_reaction_index(i_reaction)
        return self.kinetics.multiplier(i_reaction)

    def set_multiplier(self, double value, int i_reaction=-1):
        """
        Set the multiplier for for reaction *i_reaction* to *value*.
        If *i_reaction* is not specified, then the multiplier for all reactions
        is set to *value*. See `multiplier`.
        """
        if i_reaction == -1:
            for i_reaction in range(self.n_reactions):
                self.kinetics.setMultiplier(i_reaction, value)
        else:
            self._check_reaction_index(i_reaction)
            self.kinetics.setMultiplier(i_reaction, value)

    def reaction_type(self, int i_reaction):
        """Type of reaction *i_reaction*."""
        self._check_reaction_index(i_reaction)
        return self.kinetics.reactionType(i_reaction)

    def reaction_equation(self, int i_reaction):
        """The equation for the specified reaction. See also `reaction_equations`."""
        self._check_reaction_index(i_reaction)
        return pystr(self.kinetics.reactionString(i_reaction))

    def reactants(self, int i_reaction):
        """The reactants portion of the reaction equation"""
        self._check_reaction_index(i_reaction)
        return pystr(self.kinetics.reactantString(i_reaction))

    def products(self, int i_reaction):
        """The products portion of the reaction equation"""
        self._check_reaction_index(i_reaction)
        return pystr(self.kinetics.productString(i_reaction))

    def reaction_equations(self, indices=None):
        """
        Returns a list containing the reaction equation for all reactions in the
        mechanism (if *indices* is unspecified) or the equations for each
        reaction in the sequence *indices*. For example::

            >>> gas.reaction_equations()
            ['2 O + M <=> O2 + M', 'O + H + M <=> OH + M', 'O + H2 <=> H + OH', ...]
            >>> gas.reaction_equations([2,3])
            ['O + H + M <=> OH + M', 'O + H2 <=> H + OH']

        See also `reaction_equation`.
        """
        if indices is None:
            return [self.reaction_equation(i) for i in range(self.n_reactions)]
        else:
            return [self.reaction_equation(i) for i in indices]

    def reactant_stoich_coeff(self, k_spec, int i_reaction):
        """
        The stoichiometric coefficient of species *k_spec* as a reactant in
        reaction *i_reaction*.
        """
        cdef int k
        if isinstance(k_spec, (str, unicode, bytes)):
            k = self.kinetics_species_index(k_spec)
        else:
            k = k_spec
            self._check_kinetics_species_index(k_spec)

        self._check_reaction_index(i_reaction)
        return self.kinetics.reactantStoichCoeff(k, i_reaction)

    def product_stoich_coeff(self, k_spec, int i_reaction):
        """
        The stoichiometric coefficient of species *k_spec* as a product in
        reaction *i_reaction*.
        """
        cdef int k
        if isinstance(k_spec, (str, unicode, bytes)):
            k = self.kinetics_species_index(k_spec)
        else:
            k = k_spec
            self._check_kinetics_species_index(k_spec)

        self._check_reaction_index(i_reaction)
        return self.kinetics.productStoichCoeff(k, i_reaction)

    def reactant_stoich_coeffs(self):
        """
        The array of reactant stoichiometric coefficients. Element *[k,i]* of
        this array is the reactant stoichiometric coefficient of species *k* in
        reaction *i*.
        """
        cdef np.ndarray[np.double_t, ndim=2] data = np.empty((self.n_total_species,
                                                              self.n_reactions))
        cdef int i,k
        for i in range(self.n_reactions):
            for k in range(self.n_total_species):
                data[k,i] = self.kinetics.reactantStoichCoeff(k,i)
        return data

    def product_stoich_coeffs(self):
        """
        The array of product stoichiometric coefficients. Element *[k,i]* of
        this array is the product stoichiometric coefficient of species *k* in
        reaction *i*.
        """
        cdef np.ndarray[np.double_t, ndim=2] data = np.empty((self.n_total_species,
                                                              self.n_reactions))
        cdef int i,k
        for i in range(self.n_reactions):
            for k in range(self.n_total_species):
                data[k,i] = self.kinetics.productStoichCoeff(k,i)
        return data

    property forward_rates_of_progress:
        """
        Forward rates of progress for the reactions. [kmol/m^3/s] for bulk
        phases or [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRatesOfProgress)

    property reverse_rates_of_progress:
        """
        Reverse rates of progress for the reactions. [kmol/m^3/s] for bulk
        phases or [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getRevRatesOfProgress)

    property net_rates_of_progress:
        """
        Net rates of progress for the reactions. [kmol/m^3/s] for bulk phases
        or [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getNetRatesOfProgress)

    property equilibrium_constants:
        """Equilibrium constants in concentration units for all reactions."""
        def __get__(self):
            return get_reaction_array(self, kin_getEquilibriumConstants)

    property forward_rate_constants:
        """
        Forward rate constants for all reactions. Units are a combination of
        kmol, m^3 and s, that depend on the rate expression for the reaction.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRateConstants)

    property reverse_rate_constants:
        """
        Reverse rate constants for all reactions. Units are a combination of
        kmol, m^3 and s, that depend on the rate expression for the reaction.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getRevRateConstants)

    property creation_rates:
        """
        Creation rates for each species. [kmol/m^3/s] for bulk phases or
        [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_species_array(self, kin_getCreationRates)

    property destruction_rates:
        """
        Destruction rates for each species. [kmol/m^3/s] for bulk phases or
        [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_species_array(self, kin_getDestructionRates)

    property net_production_rates:
        """
        Net production rates for each species. [kmol/m^3/s] for bulk phases or
        [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_species_array(self, kin_getNetProductionRates)

    property delta_enthalpy:
        """Change in enthalpy for each reaction [J/kmol]."""
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaEnthalpy)

    property delta_gibbs:
        """Change in Gibbs free energy for each reaction [J/kmol]."""
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaGibbs)

    property delta_entropy:
        """Change in entropy for each reaction [J/kmol/K]."""
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaEntropy)

    property delta_standard_enthalpy:
        """
        Change in standard-state enthalpy (independent of composition) for
        each reaction [J/kmol].
        """
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSEnthalpy)

    property delta_standard_gibbs:
        """
        Change in standard-state Gibbs free energy (independent of composition)
        for each reaction [J/kmol].
        """
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSGibbs)

    property delta_standard_entropy:
        """
        Change in standard-state entropy (independent of composition) for
        each reaction [J/kmol/K].
        """
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSEntropy)


cdef class InterfaceKinetics(Kinetics):
    """
    A kinetics manager for heterogeneous reaction mechanisms. The
    reactions are assumed to occur at an interface between bulk phases.
    """
    def __init__(self, infile='', phaseid='', phases=(), *args, **kwargs):
        super().__init__(infile, phaseid, phases, *args, **kwargs)
        if self.kinetics.type() not in (kinetics_type_interface,
                                        kinetics_type_edge):
            raise TypeError("Underlying Kinetics class is not of the correct type.")

        self._phase_indices = {}
        for phase in [self] + list(phases):
            i = self.kinetics.phaseIndex(stringify(phase.name))
            self._phase_indices[phase] = i
            self._phase_indices[phase.name] = i
            self._phase_indices[i] = i

    def advance_coverages(self, double dt):
        """
        This method carries out a time-accurate advancement of the surface
        coverages for a specified amount of time.
        """
        (<CxxInterfaceKinetics*>self.kinetics).advanceCoverages(dt)

    def phase_index(self, phase):
        """
        Get the index of the phase *phase*, where *phase* may specified using
        the phase object, the name, or the index itself.
        """
        return self._phase_indices[phase]

    def _phase_slice(self, phase):
        p = self.phase_index(phase)
        k1 = self.kinetics_species_index(0, p)

        if p == self.n_phases - 1:
            k2 = self.n_total_species
        else:
            k2 = self.kinetics_species_index(0, p+1)

        return slice(k1,k2)

    def get_creation_rates(self, phase):
        """
        Creation rates for each species in phase *phase*. Use the
        `creation_rates` property to get the creation rates for species in all
        phases.
        """
        return self.creation_rates[self._phase_slice(phase)]

    def get_destruction_rates(self, phase):
        """
        Destruction rates for each species in phase *phase*. Use the
        `destruction_rates` property to get the destruction rates for species
        in all phases.
        """
        return self.destruction_rates[self._phase_slice(phase)]

    def get_net_production_rates(self, phase):
        """
        Net production rates for each species in phase *phase*. Use the
        `net_production_rates` property to get the net_production rates for
        species in all phases.
        """
        return self.net_production_rates[self._phase_slice(phase)]
