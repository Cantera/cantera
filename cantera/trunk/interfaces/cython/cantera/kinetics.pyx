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
    property nTotalSpecies:
        def __get__(self):
            return self.kinetics.nTotalSpecies()

    property nReactions:
        def __get__(self):
            return self.kinetics.nReactions()

    property nPhases:
        def __get__(self):
            return self.kinetics.nPhases()

    property reactionPhaseIndex:
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
        self._checkKineticsSpeciesIndex(species)
        self._checkPhaseIndex(phase)
        return self.kinetics.kineticsSpeciesIndex(species, phase)

    def isReversible(self, int iReaction):
        self._checkReactionIndex(iReaction)
        return self.kinetics.isReversible(iReaction)

    def multiplier(self, int iReaction):
        self._checkReactionIndex(iReaction)
        return self.kinetics.multiplier(iReaction)

    def setMultiplier(self, double value, int iReaction=-1):
        if iReaction == -1:
            for iReaction in range(self.nReactions):
                self.kinetics.setMultiplier(iReaction, value)
        else:
            self._checkReactionIndex(iReaction)
            self.kinetics.setMultiplier(iReaction, value)

    def reactionType(self, int iReaction):
        self._checkReactionIndex(iReaction)
        return self.kinetics.reactionType(iReaction)

    def reactionEquation(self, int iReaction):
        self._checkReactionIndex(iReaction)
        return pystr(self.kinetics.reactionString(iReaction))

    def reactionEquations(self, indices=None):
        if indices is None:
            return [self.reactionEquation(i) for i in range(self.nReactions)]
        else:
            return [self.reactionEquation(i) for i in indices]

    def reactantStoichCoeff(self, int kSpec, int iReaction):
        self._checkKineticsSpeciesIndex(kSpec)
        self._checkReactionIndex(iReaction)
        return self.kinetics.reactantStoichCoeff(kSpec, iReaction)

    def productStoichCoeff(self, int kSpec, int iReaction):
        self._checkKineticsSpeciesIndex(kSpec)
        self._checkReactionIndex(iReaction)
        return self.kinetics.productStoichCoeff(kSpec, iReaction)

    def reactantStoichCoeffs(self):
        cdef np.ndarray[np.double_t, ndim=2] data = np.empty((self.nTotalSpecies,
                                                              self.nReactions))
        cdef int i,k
        for i in range(self.nReactions):
            for k in range(self.nTotalSpecies):
                data[k,i] = self.kinetics.reactantStoichCoeff(k,i)
        return data

    def productStoichCoeffs(self):
        cdef np.ndarray[np.double_t, ndim=2] data = np.empty((self.nTotalSpecies,
                                                              self.nReactions))
        cdef int i,k
        for i in range(self.nReactions):
            for k in range(self.nTotalSpecies):
                data[k,i] = self.kinetics.productStoichCoeff(k,i)
        return data

    property fwdRatesOfProgress:
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRatesOfProgress)

    property revRatesOfProgress:
        def __get__(self):
            return get_reaction_array(self, kin_getRevRatesOfProgress)

    property netRatesOfProgress:
        def __get__(self):
            return get_reaction_array(self, kin_getNetRatesOfProgress)

    property equilibriumConstants:
        def __get__(self):
            return get_reaction_array(self, kin_getEquilibriumConstants)

    property activationEnergies:
        def __get__(self):
            return get_reaction_array(self, kin_getActivationEnergies)

    property fwdRateConstants:
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRateConstants)

    property revRateConstants:
        def __get__(self):
            return get_reaction_array(self, kin_getRevRateConstants)

    property creationRates:
        def __get__(self):
            return get_species_array(self, kin_getCreationRates)

    property destructionRates:
        def __get__(self):
            return get_species_array(self, kin_getDestructionRates)

    property netProductionRates:
        def __get__(self):
            return get_species_array(self, kin_getNetProductionRates)

    property deltaEnthalpy:
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaEnthalpy)

    property deltaGibbs:
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaGibbs)

    property deltaEntropy:
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaEntropy)

    property deltaStandardEnthalpy:
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSEnthalpy)

    property deltaStandardGibbs:
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSGibbs)

    property deltaStandardEntropy:
        def __get__(self):
            return get_reaction_array(self, kin_getDeltaSSEntropy)


cdef class InterfaceKinetics(Kinetics):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.kinetics.type() not in (kinetics_type_interface,
                                        kinetics_type_edge):
            raise TypeError("Underlying Kinetics class is not of the correct type.")

    def advanceCoverages(self, double dt):
        (<CxxInterfaceKinetics*>self.kinetics).advanceCoverages(dt)
