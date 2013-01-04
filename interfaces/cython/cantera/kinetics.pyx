ctypedef void (*kineticsMethod1d)(CxxKinetics*, double*) except +

# NOTE: These cdef functions cannot be members of Kinetics because they would
# cause "layout conflicts" when creating derived classes with multiple bases,
# e.g. class Solution. [Cython 0.16]
cdef np.ndarray get_species_array(Kinetics kin, kineticsMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(kin.nTotalSpecies)
    method(kin.kinetics, &data[0])
    # @TODO: Fix _selectedSpecies to work with interface kinetics
    if kin._selectedSpecies.size:
        return data[kin._selectedSpecies]
    else:
        return data

cdef np.ndarray get_reaction_array(Kinetics kin, kineticsMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(kin.nReactions)
    method(kin.kinetics, &data[0])
    return data


cdef class Kinetics(_SolutionBase):
    """
    Instances of class `Kinetics` are responsible for evaluating reaction rates
    of progress, species production rates, and other quantities pertaining to
    a reaction mechanism.
    """

    property nTotalSpecies:
        """
        Total number of species in all phases participating in the kinetics
        mechanism.
        """
        def __get__(self):
            return self.kinetics.nTotalSpecies()

    property nReactions:
        """Number of reactions in the reaction mechanism."""
        def __get__(self):
            return self.kinetics.nReactions()

    property nPhases:
        """Number of phases in the reaction mechanism."""
        def __get__(self):
            return self.kinetics.nPhases()

    property reactionPhaseIndex:
        """The index of the phase where the reactions occur."""
        def __get__(self):
            return self.kinetics.reactionPhaseIndex()

    def _checkPhaseIndex(self, n):
        if not 0 <= n < self.nPhases:
            raise ValueError("Phase index ({}) out of range".format(n))

    def _checkReactionIndex(self, n):
        if not 0 <= n < self.nReactions:
            raise ValueError("Reaction index ({}) out of range".format(n))

    def _checkKineticsSpeciesIndex(self, n):
        if not 0 <= n < self.nTotalSpecies:
            raise ValueError("Kinetics Species index ({}) out of range".format(n))

    def kineticsSpeciesIndex(self, int species, int phase):
        """
        The index of species *species* of phase *phase* within arrays returned
        by methods of class `Kinetics`.
        """
        self._checkKineticsSpeciesIndex(species)
        self._checkPhaseIndex(phase)
        return self.kinetics.kineticsSpeciesIndex(species, phase)

    def isReversible(self, int iReaction):
        """True if reaction `iReaction` is reversible."""
        self._checkReactionIndex(iReaction)
        return self.kinetics.isReversible(iReaction)

    def multiplier(self, int iReaction):
        """
        A scaling factor applied to the rate coefficient for reaction
        *iReaction*. Can be used to carry out sensitivity analysis or to
        selectively disable a particular reaction. See `setMultiplier`.
        """
        self._checkReactionIndex(iReaction)
        return self.kinetics.multiplier(iReaction)

    def setMultiplier(self, double value, int iReaction=-1):
        """
        Set the multiplier for for reaction *iReaction* to *value*.
        If *iReaction* is not specified, then the multiplier for all reactions
        is set to *value*. See `multiplier`.
        """
        if iReaction == -1:
            for iReaction in range(self.nReactions):
                self.kinetics.setMultiplier(iReaction, value)
        else:
            self._checkReactionIndex(iReaction)
            self.kinetics.setMultiplier(iReaction, value)

    def reactionType(self, int iReaction):
        """Type of reaction *iReaction*."""
        self._checkReactionIndex(iReaction)
        return self.kinetics.reactionType(iReaction)

    def reactionEquation(self, int iReaction):
        """The equation for the specified reaction. See also `reactionEquations`."""
        self._checkReactionIndex(iReaction)
        return pystr(self.kinetics.reactionString(iReaction))

    def reactionEquations(self, indices=None):
        """
        Returns a list containing the reaction equation for all reactions in the
        mechanism (if *indices* is unspecified) or the equations for each
        reaction in the sequence *indices*. For example::

            >>> gas.reactionEquations()
            ['2 O + M <=> O2 + M', 'O + H + M <=> OH + M', 'O + H2 <=> H + OH', ...]
            >>> gas.reactionEquations([2,3])
            ['O + H + M <=> OH + M', 'O + H2 <=> H + OH']

        See also `reactionEquation`.
        """
        if indices is None:
            return [self.reactionEquation(i) for i in range(self.nReactions)]
        else:
            return [self.reactionEquation(i) for i in indices]

    def reactantStoichCoeff(self, int kSpec, int iReaction):
        """
        The stoichiometric coefficient of species *kSpec* as a reactant in
        reaction *iReaction*.
        """
        self._checkKineticsSpeciesIndex(kSpec)
        self._checkReactionIndex(iReaction)
        return self.kinetics.reactantStoichCoeff(kSpec, iReaction)

    def productStoichCoeff(self, int kSpec, int iReaction):
        """
        The stoichiometric coefficient of species *kSpec* as a product in
        reaction *iReaction*.
        """
        self._checkKineticsSpeciesIndex(kSpec)
        self._checkReactionIndex(iReaction)
        return self.kinetics.productStoichCoeff(kSpec, iReaction)

    def reactantStoichCoeffs(self):
        """
        The array of reactant stoichiometric coefficients. Element *[k,i]* of
        this array is the reactant stoichiometric coefficient of species *k* in
        reaction *i*.
        """
        cdef np.ndarray[np.double_t, ndim=2] data = np.empty((self.nTotalSpecies,
                                                              self.nReactions))
        cdef int i,k
        for i in range(self.nReactions):
            for k in range(self.nTotalSpecies):
                data[k,i] = self.kinetics.reactantStoichCoeff(k,i)
        return data

    def productStoichCoeffs(self):
        """
        The array of product stoichiometric coefficients. Element *[k,i]* of
        this array is the product stoichiometric coefficient of species *k* in
        reaction *i*.
        """
        cdef np.ndarray[np.double_t, ndim=2] data = np.empty((self.nTotalSpecies,
                                                              self.nReactions))
        cdef int i,k
        for i in range(self.nReactions):
            for k in range(self.nTotalSpecies):
                data[k,i] = self.kinetics.productStoichCoeff(k,i)
        return data

    property fwdRatesOfProgress:
        """
        Forward rates of progress for the reactions. [kmol/m^3/s] for bulk
        phases or [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRatesOfProgress)

    property revRatesOfProgress:
        """
        Reverse rates of progress for the reactions. [kmol/m^3/s] for bulk
        phases or [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getRevRatesOfProgress)

    property netRatesOfProgress:
        """
        Net rates of progress for the reactions. [kmol/m^3/s] for bulk phases
        or [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getNetRatesOfProgress)

    property equilibriumConstants:
        """Equilibrium constants in concentration units for all reactions."""
        def __get__(self):
            return get_reaction_array(self, kin_getEquilibriumConstants)

    property activationEnergies:
        """Activation energies for all reactions [K]."""
        def __get__(self):
            return get_reaction_array(self, kin_getActivationEnergies)

    property fwdRateConstants:
        """
        Forward rate constants for all reactions. Units are a combination of
        kmol, m^3 and s, that depend on the rate expression for the reaction.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRateConstants)

    property revRateConstants:
        """
        Reverse rate constants for all reactions. Units are a combination of
        kmol, m^3 and s, that depend on the rate expression for the reaction.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getRevRateConstants)

    property creationRates:
        """
        Creation rates for each species. [kmol/m^3/s] for bulk phases or
        [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_species_array(self, kin_getCreationRates)

    property destructionRates:
        """
        Destruction rates for each species. [kmol/m^3/s] for bulk phases or
        [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_species_array(self, kin_getDestructionRates)

    property netProductionRates:
        """
        Net production rates for each species. [kmol/m^3/s] for bulk phases or
        [kmol/m^2/s] for surface phases.
        """
        def __get__(self):
            return get_species_array(self, kin_getNetProductionRates)

    property deltaEnthalpy:
        """Change in enthalpy for each reaction [J/kmol]."""
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaEnthalpy)

    property deltaGibbs:
        """Change in Gibbs free energy for each reaction [J/kmol]."""
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaGibbs)

    property deltaEntropy:
        """Change in entropy for each reaction [J/kmol/K]."""
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaEntropy)

    property deltaStandardEnthalpy:
        """
        Change in standard-state enthalpy (independent of composition) for
        each reaction [J/kmol].
        """
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSEnthalpy)

    property deltaStandardGibbs:
        """
        Change in standard-state Gibbs free energy (independent of composition)
        for each reaction [J/kmol].
        """
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSGibbs)

    property deltaStandardEntropy:
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
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.kinetics.type() not in (kinetics_type_interface,
                                        kinetics_type_edge):
            raise TypeError("Underlying Kinetics class is not of the correct type.")

    def advanceCoverages(self, double dt):
        """
        This method carries out a time-accurate advancement of the surface
        coverages for a specified amount of time.
        """
        (<CxxInterfaceKinetics*>self.kinetics).advanceCoverages(dt)
