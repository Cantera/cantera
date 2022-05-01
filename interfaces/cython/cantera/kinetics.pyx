# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from ctypes import c_int

# NOTE: These cdef functions cannot be members of Kinetics because they would
# cause "layout conflicts" when creating derived classes with multiple bases,
# such as class Solution. [Cython 0.16]
cdef np.ndarray get_species_array(Kinetics kin, kineticsMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(kin.n_total_species)
    if kin.n_total_species == 0:
        return data
    method(kin.kinetics, &data[0])
    # @todo: Fix _selected_species to work with interface kinetics
    if kin._selected_species.size:
        return data[kin._selected_species]
    else:
        return data

cdef np.ndarray get_reaction_array(Kinetics kin, kineticsMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(kin.n_reactions)
    if kin.n_reactions == 0:
        return data
    method(kin.kinetics, &data[0])
    return data

cdef np.ndarray get_dense(Kinetics kin, kineticsMethodSparse method):
    cdef CxxSparseMatrix smat = method(kin.kinetics)
    cdef size_t length = smat.nonZeros()
    if length == 0:
        return np.zeros((kin.n_reactions, 0))

    # index/value triplets
    cdef np.ndarray[int, ndim=1, mode="c"] rows = np.empty(length, dtype=c_int)
    cdef np.ndarray[int, ndim=1, mode="c"] cols = np.empty(length, dtype=c_int)
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(length)

    size = CxxSparseTriplets(smat, &rows[0], &cols[0], &data[0], length)
    out = np.zeros((smat.rows(), smat.cols()))
    for i in xrange(length):
        out[rows[i], cols[i]] = data[i]
    return out

cdef tuple get_sparse(Kinetics kin, kineticsMethodSparse method):
    # retrieve sparse matrix
    cdef CxxSparseMatrix smat = method(kin.kinetics)

    # pointers to values and inner indices of CSC storage
    cdef size_t length = smat.nonZeros()
    cdef np.ndarray[np.double_t, ndim=1] value = np.empty(length)
    cdef np.ndarray[int, ndim=1, mode="c"] inner = np.empty(length, dtype=c_int)

    # pointers outer indices of CSC storage
    cdef size_t ncols = smat.outerSize()
    cdef np.ndarray[int, ndim=1, mode="c"] outer = np.empty(ncols + 1, dtype=c_int)

    CxxSparseCscData(smat, &value[0], &inner[0], &outer[0])
    return value, inner, outer

cdef class Kinetics(_SolutionBase):
    """
    Instances of class `Kinetics` are responsible for evaluating reaction rates
    of progress, species production rates, and other quantities pertaining to
    a reaction mechanism.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self._references is None:
            raise ValueError(
                "Cannot instantiate stand-alone 'Kinetics' object as it requires an "
                "associated thermo phase.\nAll 'Kinetics' methods should be accessed "
                "from a 'Solution' object.")

    property kinetics_model:
        """
        Return type of kinetics.
        """
        def __get__(self):
            return pystr(self.kinetics.kineticsType())

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
        The index of species ``species`` of phase ``phase`` within arrays returned
        by methods of class `Kinetics`. If ``species`` is a string, the ``phase``
        argument is unused.
        """
        cdef int k
        if isinstance(species, (str, bytes)):
            return self.kinetics.kineticsSpeciesIndex(stringify(species))
        else:
            k = species
            self._check_kinetics_species_index(k)
            self._check_phase_index(k)
            return self.kinetics.kineticsSpeciesIndex(k, phase)

    def kinetics_species_name(self, k):
        """
        Name of the species with index ``k`` in the arrays returned by methods
        of class `Kinetics`.
        """
        return pystr(self.kinetics.kineticsSpeciesName(k))

    property kinetics_species_names:
        """
        A list of all species names, corresponding to the arrays returned by
        methods of class `Kinetics`.
        """
        def __get__(self):
            return [self.kinetics_species_name(k)
                    for k in range(self.n_total_species)]

    def reaction(self, int i_reaction):
        """
        Return a `Reaction` object representing the reaction with index
        ``i_reaction``. Changes to this object do not affect the `Kinetics` or
        `Solution` object until the `modify_reaction` function is called.
        """
        return Reaction.wrap(self.kinetics.reaction(i_reaction))

    def reactions(self):
        """
        Return a list of all `Reaction` objects. Changes to these objects do not
        affect the `Kinetics` or `Solution` object until the `modify_reaction`
        function is called.
        """
        return [self.reaction(i) for i in range(self.n_reactions)]

    def modify_reaction(self, int irxn, Reaction rxn):
        """
        Modify the `Reaction` with index ``irxn`` to have the same rate
        parameters as ``rxn``. ``rxn`` must have the same reactants and products
        and be of the same type (for example, `ElementaryReaction`, `FalloffReaction`,
        `PlogReaction`, etc.) as the existing reaction. This method does not
        modify the third-body efficiencies, reaction orders, or reversibility of
        the reaction.
        """
        self.kinetics.modifyReaction(irxn, rxn._reaction)

    def add_reaction(self, Reaction rxn):
        """ Add a new reaction to this phase. """
        self.kinetics.addReaction(rxn._reaction)

    def is_reversible(self, int i_reaction):
        """
        True if reaction ``i_reaction`` is reversible.

        .. deprecated:: 2.6

            Replaced by property `Reaction.reversible`.
            Example: ``gas.is_reversible(0)`` is replaced by
            ``gas.reaction(0).reversible``
        """
        rxn = self.reaction(i_reaction)
        warnings.warn(
            "'is_reversible' is deprecated and will be removed after Cantera 2.6.\n"
            "Replaceable by property 'reversible' of the corresponding "
            "reaction object.", DeprecationWarning)
        return rxn.reversible

    def multiplier(self, int i_reaction):
        """
        A scaling factor applied to the rate coefficient for reaction
        ``i_reaction``. Can be used to carry out sensitivity analysis or to
        selectively disable a particular reaction. See `set_multiplier`.
        """
        self._check_reaction_index(i_reaction)
        return self.kinetics.multiplier(i_reaction)

    def set_multiplier(self, double value, int i_reaction=-1):
        """
        Set the multiplier for for reaction ``i_reaction`` to ``value``.
        If ``i_reaction`` is not specified, then the multiplier for all reactions
        is set to ``value``. See `multiplier`.
        """
        if i_reaction == -1:
            for i_reaction in range(self.n_reactions):
                self.kinetics.setMultiplier(i_reaction, value)
        else:
            self._check_reaction_index(i_reaction)
            self.kinetics.setMultiplier(i_reaction, value)

    def reaction_type(self, int i_reaction):
        """
        Type code of reaction ``i_reaction``.

        .. deprecated:: 2.6

            Replaced by properties `Reaction.type` and `Reaction.rate.type`.
            Example: ``gas.reaction_type(0)`` is replaced by
            ``gas.reaction(0).reaction_type`` and ``gas.reaction(0).rate.type``
        """
        rxn = self.reaction(i_reaction)
        if not rxn.uses_legacy:
            warnings.warn(
                "'reaction_type' is deprecated and will be removed after "
                "Cantera 2.6.\nReplaceable by property 'reaction_type' of the "
                "corresponding reaction object (or property 'type' of the\n"
                "associated 'rate').", DeprecationWarning)
        return rxn.type

    def reaction_equation(self, int i_reaction):
        """
        The equation for the specified reaction. See also `reaction_equations`.

        .. deprecated:: 2.6

            Replaced by property `Reaction.equation`.
            Example: ``gas.reaction_equation(0)`` is replaced by
            ``gas.reaction(0).equation``
        """
        rxn = self.reaction(i_reaction)
        warnings.warn(
            "'reaction_equation' is deprecated and will be removed after "
            "Cantera 2.6.\nReplaceable by property 'equation' of the corresponding "
            "reaction object.", DeprecationWarning)
        return rxn.equation

    def reactants(self, int i_reaction):
        """
        The reactants portion of the reaction equation

        .. deprecated:: 2.6

            Replaced by property `Reaction.reactants`.
            Example: ``gas.reactants(0)`` is replaced by ``gas.reaction(0).reactants``
        """
        rxn = self.reaction(i_reaction)
        warnings.warn(
            "'reactants' is deprecated and will be removed after Cantera 2.6.\n"
            "Replaceable by property 'reactant_string' of the corresponding "
            "reaction object.", DeprecationWarning)
        return rxn.reactant_string

    def products(self, int i_reaction):
        """
        The products portion of the reaction equation

        .. deprecated:: 2.6

            Replaced by property `Reaction.products`.
            Example: ``gas.products(0)`` is replaced by ``gas.reaction(0).products``
        """
        rxn = self.reaction(i_reaction)
        warnings.warn(
            "'products' is deprecated and will be removed after Cantera 2.6.\n"
            "Replaceable by property 'product_string' of the corresponding "
            "reaction object.", DeprecationWarning)
        return rxn.product_string

    def reaction_equations(self, indices=None):
        """
        Returns a list containing the reaction equation for all reactions in the
        mechanism if ``indices`` is unspecified, or the equations for each
        reaction in the sequence ``indices``. For example::

            >>> gas.reaction_equations()
            ['2 O + M <=> O2 + M', 'O + H + M <=> OH + M', 'O + H2 <=> H + OH', ...]
            >>> gas.reaction_equations([2,3])
            ['O + H + M <=> OH + M', 'O + H2 <=> H + OH']

        See also `reaction_equation`.
        """
        if indices is None:
            return list(r.equation for r in self.reactions())
        else:
            return [self.reaction(i).equation for i in indices]

    def reactant_stoich_coeff(self, k_spec, int i_reaction):
        """
        The stoichiometric coefficient of species ``k_spec`` as a reactant in
        reaction ``i_reaction``.
        """
        cdef int k
        if isinstance(k_spec, (str, bytes)):
            k = self.kinetics_species_index(k_spec)
        else:
            k = k_spec
            self._check_kinetics_species_index(k_spec)

        self._check_reaction_index(i_reaction)
        return self.kinetics.reactantStoichCoeff(k, i_reaction)

    def product_stoich_coeff(self, k_spec, int i_reaction):
        """
        The stoichiometric coefficient of species ``k_spec`` as a product in
        reaction ``i_reaction``.
        """
        cdef int k
        if isinstance(k_spec, (str, bytes)):
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

        .. deprecated:: 2.6

            Behavior to change after Cantera 2.6; for new behavior, see property
            `Kinetics.reactant_stoich_coeffs3`.
        """
        warnings.warn("Behavior to change after Cantera 2.6; for new behavior, see "
                      "property 'reactant_stoich_coeffs3'.", DeprecationWarning)
        return self.reactant_stoich_coeffs3

    property reactant_stoich_coeffs3:
        """
        The array of reactant stoichiometric coefficients. Element ``[k,i]`` of
        this array is the reactant stoichiometric coefficient of species ``k`` in
        reaction ``i``.

        For sparse output, set ``ct.use_sparse(True)``.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_reactantStoichCoeffs)
                shape = self.n_total_species, self.n_reactions
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_reactantStoichCoeffs)

    def product_stoich_coeffs(self):
        """
        The array of product stoichiometric coefficients. Element *[k,i]* of
        this array is the product stoichiometric coefficient of species *k* in
        reaction *i*.

        .. deprecated:: 2.6

            Behavior to change after Cantera 2.6; for new behavior, see property
            `Kinetics.reactant_stoich_coeffs3`.
        """
        warnings.warn("Behavior to change after Cantera 2.6; for new behavior, see "
                      "property 'product_stoich_coeffs3'.", DeprecationWarning)
        return self.product_stoich_coeffs3

    property product_stoich_coeffs3:
        """
        The array of product stoichiometric coefficients. Element ``[k,i]`` of
        this array is the product stoichiometric coefficient of species ``k`` in
        reaction ``i``.

        For sparse output, set ``ct.use_sparse(True)``.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_productStoichCoeffs)
                shape = self.n_total_species, self.n_reactions
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_productStoichCoeffs)

    property product_stoich_coeffs_reversible:
        """
        The array of product stoichiometric coefficients of reversible reactions.
        Element ``[k,i]`` of this array is the product stoichiometric coefficient
        of species ``k`` in reaction ``i``.

        For sparse output, set ``ct.use_sparse(True)``.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_revProductStoichCoeffs)
                shape = self.n_total_species, self.n_reactions
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_revProductStoichCoeffs)

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
        Forward rate constants for all reactions. The computed values include
        all temperature-dependent, pressure-dependent, and third body
        contributions. Units are a combination of kmol, m^3 and s, that depend
        on the rate expression for the reaction.

        .. deprecated:: 2.6

            Behavior to change after Cantera 2.6; for Cantera 2.6, rate constants of
            three-body reactions are multiplied with third-body concentrations
            (no change to legacy behavior). After Cantera 2.6, results will no longer
            include third-body concentrations and be consistent with conventional
            definitions (see Eq. 9.75 in Kee, Coltrin, and Glarborg, *Chemically
            Reacting Flow*, Wiley Interscience, 2003).
            To switch to new behavior, run ``ct.use_legacy_rate_constants(False)``.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRateConstants)

    property reverse_rate_constants:
        """
        Reverse rate constants for all reactions. The computed values include
        all temperature-dependent, pressure-dependent, and third body
        contributions. Units are a combination of kmol, m^3 and s, that depend
        on the rate expression for the reaction.

        .. deprecated:: 2.6

            Behavior to change after Cantera 2.6; for Cantera 2.6, rate constants of
            three-body reactions are multiplied with third-body concentrations
            (no change to legacy behavior). After Cantera 2.6, results will no longer
            include third-body concentrations and be consistent with conventional
            definitions (see Eq. 9.75 in Kee, Coltrin and Glarborg, *Chemically
            Reacting Flow*, Wiley Interscience, 2003).
            To switch to new behavior, run ``ct.use_legacy_rate_constants(False)``.
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

    property derivative_settings:
        """
        Property setting behavior of derivative evaluation.

        For ``GasKinetics``, the following keyword/value pairs are supported:

        -  ``skip-third-bodies`` (boolean) ... if `False` (default), third body
           concentrations are considered for the evaluation of derivatives

        -  ``skip-falloff`` (boolean) ... if `True` (default), third-body effects
           on reaction rates are not considered.

        -  ``rtol-delta`` (double) ... relative tolerance used to perturb properties
           when calculating numerical derivatives. The default value is 1e-8.

        Derivative settings are updated using a dictionary::

            >>> gas.derivative_settings = {"skip-falloff": True}

        Passing an empty dictionary will reset all values to their defaults.
        """
        def __get__(self):
            cdef CxxAnyMap settings
            self.kinetics.getDerivativeSettings(settings)
            return anymap_to_dict(settings)
        def __set__(self, settings):
            self.kinetics.setDerivativeSettings(dict_to_anymap(settings))

    property forward_rate_constants_ddT:
        """
        Calculate derivatives for forward rate constants with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRateConstants_ddT)

    property forward_rate_constants_ddP:
        """
        Calculate derivatives for forward rate constants with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRateConstants_ddP)

    property forward_rate_constants_ddC:
        """
        Calculate derivatives for forward rate constants with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRateConstants_ddC)

    property forward_rates_of_progress_ddT:
        """
        Calculate derivatives for forward rates-of-progress with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRatesOfProgress_ddT)

    property forward_rates_of_progress_ddP:
        """
        Calculate derivatives for forward rates-of-progress with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRatesOfProgress_ddP)

    property forward_rates_of_progress_ddC:
        """
        Calculate derivatives for forward rates-of-progress with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getFwdRatesOfProgress_ddC)

    property forward_rates_of_progress_ddX:
        """
        Calculate derivatives for forward rates-of-progress with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_fwdRatesOfProgress_ddX)
                shape = self.n_reactions, self.n_total_species
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_fwdRatesOfProgress_ddX)

    property reverse_rates_of_progress_ddT:
        """
        Calculate derivatives for reverse rates-of-progress with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getRevRatesOfProgress_ddT)

    property reverse_rates_of_progress_ddP:
        """
        Calculate derivatives for reverse rates-of-progress with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getRevRatesOfProgress_ddP)

    property reverse_rates_of_progress_ddC:
        """
        Calculate derivatives for reverse rates-of-progress with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getRevRatesOfProgress_ddC)

    property reverse_rates_of_progress_ddX:
        """
        Calculate derivatives for reverse rates-of-progress with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_revRatesOfProgress_ddX)
                shape = self.n_reactions, self.n_total_species
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_revRatesOfProgress_ddX)

    property net_rates_of_progress_ddT:
        """
        Calculate derivatives for net rates-of-progress with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getNetRatesOfProgress_ddT)

    property net_rates_of_progress_ddP:
        """
        Calculate derivatives for net rates-of-progress with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getNetRatesOfProgress_ddP)

    property net_rates_of_progress_ddC:
        """
        Calculate derivatives for net rates-of-progress with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getNetRatesOfProgress_ddC)

    property net_rates_of_progress_ddX:
        """
        Calculate derivatives for net rates-of-progress with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_netRatesOfProgress_ddX)
                shape = self.n_reactions, self.n_total_species
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_netRatesOfProgress_ddX)

    property creation_rates_ddT:
        """
        Calculate derivatives of species creation rates with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_species_array(self, kin_getCreationRates_ddT)

    property creation_rates_ddP:
        """
        Calculate derivatives of species creation rates with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_species_array(self, kin_getCreationRates_ddP)

    property creation_rates_ddC:
        """
        Calculate derivatives of species creation rates with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            return get_species_array(self, kin_getCreationRates_ddC)

    property creation_rates_ddX:
        """
        Calculate derivatives for species creation rates with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_creationRates_ddX)
                shape = self.n_total_species, self.n_total_species
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_creationRates_ddX)

    property destruction_rates_ddT:
        """
        Calculate derivatives of species destruction rates with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_species_array(self, kin_getDestructionRates_ddT)

    property destruction_rates_ddP:
        """
        Calculate derivatives of species destruction rates with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_species_array(self, kin_getDestructionRates_ddP)

    property destruction_rates_ddC:
        """
        Calculate derivatives of species destruction rates with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            return get_species_array(self, kin_getDestructionRates_ddC)

    property destruction_rates_ddX:
        """
        Calculate derivatives for species destruction rates with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_destructionRates_ddX)
                shape = self.n_total_species, self.n_total_species
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_destructionRates_ddX)

    property net_production_rates_ddT:
        """
        Calculate derivatives of species net production rates with respect to
        temperature at constant pressure, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_species_array(self, kin_getNetProductionRates_ddT)

    property net_production_rates_ddP:
        """
        Calculate derivatives of species net production rates with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        def __get__(self):
            return get_species_array(self, kin_getNetProductionRates_ddP)

    property net_production_rates_ddC:
        """
        Calculate derivatives of species net production rates with respect to molar
        density at constant temperature, pressure and mole fractions.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            return get_species_array(self, kin_getNetProductionRates_ddC)

    property net_production_rates_ddX:
        """
        Calculate derivatives for species net production rates with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        **Warning:** this property is an experimental part of the Cantera API and
        may be changed or removed without notice.
        """
        def __get__(self):
            if _USE_SPARSE:
                tup = get_sparse(self, kin_netProductionRates_ddX)
                shape = self.n_total_species, self.n_total_species
                return _scipy_sparse.csc_matrix(tup, shape=shape)
            return get_dense(self, kin_netProductionRates_ddX)

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

    property third_body_concentrations:
        """
        Effective third-body concentrations used by individual reactions; values
        are only defined for reactions involving third-bodies and are set to
        not-a-number otherwise.
        """
        def __get__(self):
            return get_reaction_array(self, kin_getThirdBodyConcentrations)

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

    property heat_release_rate:
        """
        Get the total volumetric heat release rate [W/m^3].
        """
        def __get__(self):
            return - np.sum(self.partial_molar_enthalpies *
                            self.net_production_rates, 0)

    property heat_production_rates:
        """
        Get the volumetric heat production rates [W/m^3] on a per-reaction
        basis. The sum over all reactions results in the total volumetric heat
        release rate.
        Example: C. K. Law: Combustion Physics (2006), Fig. 7.8.6

        >>> gas.heat_production_rates[1]  # heat production rate of the 2nd reaction
        """
        def __get__(self):
            return - self.net_rates_of_progress * self.delta_enthalpy


cdef class InterfaceKinetics(Kinetics):
    """
    A kinetics manager for heterogeneous reaction mechanisms. The
    reactions are assumed to occur at an interface between bulk phases.
    """
    def __init__(self, infile='', name='', adjacent=(), *args, **kwargs):
        super().__init__(infile, name, *args, **kwargs)
        if not kwargs.get("init", True):
            return
        if pystr(self.kinetics.kineticsType()) not in ("Surf", "Edge"):
            raise TypeError("Underlying Kinetics class is not of the correct type.")
        self._setup_phase_indices()

    def _setup_phase_indices(self):
        self._phase_indices = {}
        for name, phase in list(self.adjacent.items()) + [(self.name, self)]:
            i = self.kinetics.phaseIndex(stringify(name))
            self._phase_indices[phase] = i
            self._phase_indices[name] = i
            self._phase_indices[i] = i

    def advance_coverages(self, double dt, double rtol=1e-7, double atol=1e-14,
                          double max_step_size=0, int max_steps=20000,
                          int max_error_test_failures=7):
        """
        This method carries out a time-accurate advancement of the surface
        coverages for a specified amount of time.
        """
        (<CxxInterfaceKinetics*>self.kinetics).advanceCoverages(
            dt, rtol, atol, max_step_size, max_steps, max_error_test_failures)

    def advance_coverages_to_steady_state(self):
        """
        This method advances the surface coverages to steady state.
        """
        (<CxxInterfaceKinetics*>self.kinetics).solvePseudoSteadyStateProblem()

    def phase_index(self, phase):
        """
        Get the index of the phase ``phase``, where ``phase`` may specified using
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
        Creation rates for each species in phase ``phase``. Use the
        `creation_rates` property to get the creation rates for species in all
        phases.
        """
        return self.creation_rates[self._phase_slice(phase)]

    def get_destruction_rates(self, phase):
        """
        Destruction rates for each species in phase ``phase``. Use the
        `destruction_rates` property to get the destruction rates for species
        in all phases.
        """
        return self.destruction_rates[self._phase_slice(phase)]

    def get_net_production_rates(self, phase):
        """
        Net production rates for each species in phase ``phase``. Use the
        `net_production_rates` property to get the net production rates for
        species in all phases.
        """
        return self.net_production_rates[self._phase_slice(phase)]

    def write_yaml(self, filename, phases=None, units=None, precision=None,
                   skip_user_defined=None):
        """
        See `_SolutionBase.write_yaml <cantera._cantera._SolutionBase.write_yaml>`.
        """
        if phases is not None:
            phases = list(phases)
        else:
            phases = []

        for phase in self._phase_indices:
            if isinstance(phase, _SolutionBase) and phase is not self:
                phases.append(phase)

        super().write_yaml(filename, phases, units, precision, skip_user_defined)
