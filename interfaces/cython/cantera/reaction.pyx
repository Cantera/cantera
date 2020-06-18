# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

cdef extern from "cantera/kinetics/reaction_defs.h" namespace "Cantera":
    cdef int ELEMENTARY_RXN
    cdef int THREE_BODY_RXN
    cdef int FALLOFF_RXN
    cdef int PLOG_RXN
    cdef int CHEBYSHEV_RXN
    cdef int CHEMACT_RXN
    cdef int INTERFACE_RXN


cdef class Reaction:
    """
    A class which stores data about a reaction and its rate parameterization so
    that it can be added to a `Kinetics` object.

    :param reactants:
        Value used to set `reactants`
    :param products:
        Value used to set `products`

    The static methods `listFromFile`, `listFromYaml`, `listFromCti`, and
    `listFromXml` can be used to create lists of `Reaction` objects from
    existing definitions in the YAML, CTI, or XML formats. All of the following
    will produce a list of the 325 reactions which make up the GRI 3.0
    mechanism::

        R = ct.Reaction.listFromFile('gri30.yaml', gas)
        R = ct.Reaction.listFromCti(open('path/to/gri30.cti').read())
        R = ct.Reaction.listFromXml(open('path/to/gri30.xml').read())

    where `gas` is a `Solution` object with the appropriate thermodynamic model,
    which is the `ideal-gas` model in this case.

    The static method `listFromYaml` can be used to create lists of `Reaction`
    objects from a YAML list::

        rxns = '''
          - equation: O + H2 <=> H + OH
            rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
          - equation: O + HO2 <=> OH + O2
            rate-constant: {A: 2.0e+13, b: 0.0, Ea: 0.0}
        '''
        R = ct.Reaction.listFromYaml(rxns, gas)

    The methods `fromYaml`, `fromCti`, and `fromXml` can be used to create
    individual `Reaction` objects from definitions in these formats. In the case
    of using YAML or CTI definitions, it is important to verify that either the
    pre-exponential factor and activation energy are supplied in SI units, or
    that they have their units specified::

        R = ct.Reaction.fromYaml('''{equation: O + H2 <=> H + OH,
                rate-constant: {A: 3.87e+04 cm^3/mol/s, b: 2.7, Ea: 6260 cal/mol}}''',
                gas)

        R = ct.Reaction.fromCti('''reaction('O + H2 <=> H + OH',
                [3.87e1, 2.7, 2.619184e7])''')

        R = ct.Reaction.fromCti('''reaction('O + H2 <=> H + OH',
                        [(3.87e4, 'cm3/mol/s'), 2.7, (6260, 'cal/mol')])''')
    """
    reaction_type = 0

    def __cinit__(self, reactants='', products='', init=True, **kwargs):
        if init:
            self._reaction.reset(newReaction(self.reaction_type))
            self.reaction = self._reaction.get()
            if reactants:
                self.reactants = reactants
            if products:
                self.products = products

    cdef _assign(self, shared_ptr[CxxReaction] other):
        self._reaction = other
        self.reaction = self._reaction.get()

    @staticmethod
    def fromCti(text):
        """
        Create a Reaction object from its CTI string representation.

        .. deprecated:: 2.5

            The CTI input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFromString(stringify(text))))
        assert cxx_reactions.size() == 1, cxx_reactions.size()
        return wrapReaction(cxx_reactions[0])

    @staticmethod
    def fromXml(text):
        """
        Create a Reaction object from its XML string representation.

        .. deprecated:: 2.5

            The XML input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_reaction = CxxNewReaction(deref(CxxGetXmlFromString(stringify(text))))
        return wrapReaction(cxx_reaction)

    @staticmethod
    def fromYaml(text, Kinetics kinetics):
        """
        Create a `Reaction` object from its YAML string representation.

        :param text:
            The YAML reaction string
        :param kinetics:
            A `Kinetics` object whose associated phase(s) contain the species
            involved in the reaction.
        """
        cxx_reaction = CxxNewReaction(AnyMapFromYamlString(stringify(text)),
                                      deref(kinetics.kinetics))
        return wrapReaction(cxx_reaction)

    @staticmethod
    def listFromFile(filename, Kinetics kinetics=None, section='reactions'):
        """
        Create a list of Reaction objects from all of the reactions defined in a
        YAML, CTI, or XML file.

        For YAML input files, a `Kinetics` object is required as the second
        argument, and reactions from the section *section* will be returned.

        Directories on Cantera's input file path will be searched for the
        specified file.

        In the case of an XML file, the ``<reactions>`` nodes are assumed to be
        children of the ``<reactionsData>`` node in a document with a ``<ctml>``
        root node, as in the XML files produced by conversion from CTI files.

        .. deprecated:: 2.5

            The CTI and XML input formats are deprecated and will be removed in
            Cantera 3.0.
        """
        if filename.lower().split('.')[-1] in ('yml', 'yaml'):
            if kinetics is None:
                raise ValueError("A Kinetics object is required.")
            root = AnyMapFromYamlFile(stringify(filename))
            cxx_reactions = CxxGetReactions(root[stringify(section)],
                                            deref(kinetics.kinetics))
        else:
            cxx_reactions = CxxGetReactions(deref(CxxGetXmlFile(stringify(filename))))
        return [wrapReaction(r) for r in cxx_reactions]

    @staticmethod
    def listFromXml(text):
        """
        Create a list of Reaction objects from all the reaction defined in an
        XML string. The ``<reaction>`` nodes are assumed to be children of the
        ``<reactionData>`` node in a document with a ``<ctml>`` root node, as in
        the XML files produced by conversion from CTI files.

        .. deprecated:: 2.5

            The XML input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFromString(stringify(text))))
        return [wrapReaction(r) for r in cxx_reactions]

    @staticmethod
    def listFromCti(text):
        """
        Create a list of `Reaction` objects from all the reactions defined in a
        CTI string.

        .. deprecated:: 2.5

            The CTI input format is deprecated and will be removed in Cantera 3.0.
        """
        # Currently identical to listFromXml since get_XML_from_string is able
        # to distinguish between CTI and XML.
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFromString(stringify(text))))
        return [wrapReaction(r) for r in cxx_reactions]

    @staticmethod
    def listFromYaml(text, Kinetics kinetics):
        """
        Create a list of `Reaction` objects from all the reactions defined in a
        YAML string.
        """
        root = AnyMapFromYamlString(stringify(text))
        cxx_reactions = CxxGetReactions(root[stringify("items")],
                                        deref(kinetics.kinetics))
        return [wrapReaction(r) for r in cxx_reactions]

    property reactant_string:
        """
        A string representing the reactants side of the chemical equation for
        this reaction. Determined automatically based on `reactants`.
        """
        def __get__(self):
            return pystr(self.reaction.reactantString())

    property product_string:
        """
        A string representing the products side of the chemical equation for
        this reaction. Determined automatically based on `products`.
        """
        def __get__(self):
            return pystr(self.reaction.productString())

    property equation:
        """
        A string giving the chemical equation for this reaction. Determined
        automatically based on `reactants` and `products`.
        """
        def __get__(self):
            return pystr(self.reaction.equation())

    property reactants:
        """
        Get/Set the reactants in this reaction as a dict where the keys are
        species names and the values, are the stoichiometric coefficients, e.g.
        ``{'CH4':1, 'OH':1}``, or as a composition string, e.g.
        ``'CH4:1, OH:1'``.
        """
        def __get__(self):
            return comp_map_to_dict(self.reaction.reactants)
        def __set__(self, reactants):
            self.reaction.reactants = comp_map(reactants)

    property products:
        """
        Get/Set the products in this reaction as a dict where the keys are
        species names and the values, are the stoichiometric coefficients, e.g.
        ``{'CH3':1, 'H2O':1}``, or as a composition string, e.g.
        ``'CH3:1, H2O:1'``.
        """
        def __get__(self):
            return comp_map_to_dict(self.reaction.products)
        def __set__(self, products):
            self.reaction.products = comp_map(products)

    def __contains__(self, species):
        return species in self.reactants or species in self.products

    property orders:
        """
        Get/Set the reaction order with respect to specific species as a dict
        with species names as the keys and orders as the values, or as a
        composition string. By default, mass-action kinetics is assumed, with
        the reaction order for each reactant species equal to each its
        stoichiometric coefficient.
        """
        def __get__(self):
            return comp_map_to_dict(self.reaction.orders)
        def __set__(self, orders):
            self.reaction.orders = comp_map(orders)

    property ID:
        """
        Get/Set the identification string for the reaction, which can be used in
        filtering operations.
        """
        def __get__(self):
            return pystr(self.reaction.id)
        def __set__(self, ID):
            self.reaction.id = stringify(ID)

    property reversible:
        """
        Get/Set a flag which is `True` if this reaction is reversible or `False`
        otherwise.
        """
        def __get__(self):
            return self.reaction.reversible
        def __set__(self, reversible):
            self.reaction.reversible = reversible

    property duplicate:
        """
        Get/Set a flag which is `True` if this reaction is marked as a duplicate
        or `False` otherwise.
        """
        def __get__(self):
            return self.reaction.duplicate
        def __set__(self, duplicate):
             self.reaction.duplicate = duplicate

    property allow_nonreactant_orders:
        """
        Get/Set a flag which is `True` if reaction orders can be specified for
        non-reactant species. Default is `False`.
        """
        def __get__(self):
            return self.reaction.allow_nonreactant_orders
        def __set__(self, allow):
            self.reaction.allow_nonreactant_orders = allow

    property allow_negative_orders:
        """
        Get/Set a flag which is `True` if negative reaction orders are allowed.
        Default is `False`.
        """
        def __get__(self):
            return self.reaction.allow_negative_orders
        def __set__(self, allow):
            self.reaction.allow_negative_orders = allow

    def __repr__(self):
        return '<{}: {}>'.format(self.__class__.__name__, self.equation)

    def __str__(self):
        return self.equation


cdef class Arrhenius:
    r"""
    A reaction rate coefficient which depends on temperature only and follows
    the modified Arrhenius form:

    .. math::

        k_f = A T^b \exp{-\tfrac{E}{RT}}

    where *A* is the `pre_exponential_factor`, *b* is the `temperature_exponent`,
    and *E* is the `activation_energy`.
    """
    def __cinit__(self, A=0, b=0, E=0, init=True):
        if init:
            self.rate = new CxxArrhenius(A, b, E / gas_constant)
            self.reaction = None

    def __dealloc__(self):
        if self.reaction is None:
            del self.rate

    property pre_exponential_factor:
        """
        The pre-exponential factor *A* in units of m, kmol, and s raised to
        powers depending on the reaction order.
        """
        def __get__(self):
            return self.rate.preExponentialFactor()

    property temperature_exponent:
        """
        The temperature exponent *b*.
        """
        def __get__(self):
            return self.rate.temperatureExponent()

    property activation_energy:
        """
        The activation energy *E* [J/kmol].
        """
        def __get__(self):
            return self.rate.activationEnergy_R() * gas_constant

    def __repr__(self):
        return 'Arrhenius(A={:g}, b={:g}, E={:g})'.format(
            self.pre_exponential_factor, self.temperature_exponent,
            self.activation_energy)

    def __call__(self, float T):
        cdef double logT = np.log(T)
        cdef double recipT = 1/T
        return self.rate.updateRC(logT, recipT)


cdef wrapArrhenius(CxxArrhenius* rate, Reaction reaction):
    r = Arrhenius(init=False)
    r.rate = rate
    r.reaction = reaction
    return r

cdef copyArrhenius(CxxArrhenius* rate):
    r = Arrhenius(rate.preExponentialFactor(), rate.temperatureExponent(),
                  rate.activationEnergy_R() * gas_constant)
    return r


cdef class ElementaryReaction(Reaction):
    """
    A reaction which follows mass-action kinetics with a modified Arrhenius
    reaction rate.
    """
    reaction_type = ELEMENTARY_RXN

    property rate:
        """ Get/Set the `Arrhenius` rate coefficient for this reaction. """
        def __get__(self):
            cdef CxxElementaryReaction* r = <CxxElementaryReaction*>self.reaction
            return wrapArrhenius(&(r.rate), self)
        def __set__(self, Arrhenius rate):
            cdef CxxElementaryReaction* r = <CxxElementaryReaction*>self.reaction
            r.rate = deref(rate.rate)

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            cdef CxxElementaryReaction* r = <CxxElementaryReaction*>self.reaction
            return r.allow_negative_pre_exponential_factor
        def __set__(self, allow):
            cdef CxxElementaryReaction* r = <CxxElementaryReaction*>self.reaction
            r.allow_negative_pre_exponential_factor = allow


cdef class ThreeBodyReaction(ElementaryReaction):
    """
    A reaction with a non-reacting third body "M" that acts to add or remove
    energy from the reacting species.
    """
    reaction_type = THREE_BODY_RXN

    cdef CxxThreeBodyReaction* tbr(self):
        return <CxxThreeBodyReaction*>self.reaction

    property efficiencies:
        """
        Get/Set a `dict` defining non-default third-body efficiencies for this
        reaction, where the keys are the species names and the values are the
        efficiencies.
        """
        def __get__(self):
            return comp_map_to_dict(self.tbr().third_body.efficiencies)
        def __set__(self, eff):
            self.tbr().third_body.efficiencies = comp_map(eff)

    property default_efficiency:
        """
        Get/Set the default third-body efficiency for this reaction, used for
        species used for species not in `efficiencies`.
        """
        def __get__(self):
            return self.tbr().third_body.default_efficiency
        def __set__(self, default_eff):
            self.tbr().third_body.default_efficiency = default_eff

    def efficiency(self, species):
        """
        Get the efficiency of the third body named *species* considering both
        the default efficiency and species-specific efficiencies.
        """
        return self.tbr().third_body.efficiency(stringify(species))


cdef class Falloff:
    """
    A parameterization used to describe the fall-off in reaction rate constants
    due to intermolecular energy transfer. These functions are used by reactions
    defined using the `FalloffReaction` and `ChemicallyActivatedReaction`
    classes.

    This base class implements the simple falloff function
    :math:`F(T,P_r) = 1.0`.

    :param params:
        Not used for the "simple" falloff parameterization.
    :param init:
        Used internally when wrapping :ct:`Falloff` objects returned from C++.
    """
    falloff_type = "Lindemann"

    def __cinit__(self, params=(), init=True):
        if not init:
            return

        cdef vector[double] c
        for p in params:
            c.push_back(p)
        self._falloff = CxxNewFalloff(stringify(self.falloff_type), c)
        self.falloff = self._falloff.get()

    property type:
        """ A string defining the type of the falloff parameterization """
        def __get__(self):
            return pystr(self.falloff.type())

    property parameters:
        """ The array of parameters used to define this falloff function. """
        def __get__(self):
            N = self.falloff.nParameters()
            if N == 0:
                return np.empty(0)

            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(N)
            self.falloff.getParameters(&data[0])
            return data

    def __call__(self, float T, float Pr):
        """ Evaluate the falloff function :math:`F(T, P_r)` """
        N = max(self.falloff.workSize(), 1)
        cdef np.ndarray[np.double_t, ndim=1] work = np.empty(N)
        self.falloff.updateTemp(T, &work[0])
        return self.falloff.F(Pr, &work[0])


cdef class TroeFalloff(Falloff):
    """
    The 3- or 4-parameter Troe falloff function.

    :param params:
        An array of 3 or 4 parameters: :math:`[a, T^{***}, T^*, T^{**}]` where
        the final parameter is optional (with a default value of 0).
    """
    falloff_type = "Troe"


cdef class SriFalloff(Falloff):
    """
    The 3- or 5-parameter SRI falloff function.

    :param params:
        An array of 3 or 5 parameters: :math:`[a, b, c, d, e]` where the last
        two parameters are optional (with default values of 1 and 0,
        respectively).
    """
    falloff_type = "SRI"


cdef wrapFalloff(shared_ptr[CxxFalloff] falloff):
    falloff_type = pystr(falloff.get().type())
    if falloff_type in ["Lindemann", "Simple"]:
        f = Falloff(init=False)
    elif falloff_type == "Troe":
        f = TroeFalloff(init=False)
    elif falloff_type == "SRI":
        f = SriFalloff(init=False)
    else:
        warnings.warn('Unknown falloff type: {0}'.format(falloff_type))
        f = Falloff(init=False)
    f._falloff = falloff
    f.falloff = f._falloff.get()
    return f


cdef class FalloffReaction(Reaction):
    """
    A reaction that is first-order in [M] at low pressure, like a third-body
    reaction, but zeroth-order in [M] as pressure increases.
    """
    reaction_type = FALLOFF_RXN

    cdef CxxFalloffReaction* frxn(self):
        return <CxxFalloffReaction*>self.reaction

    property low_rate:
        """ Get/Set the `Arrhenius` rate constant in the low-pressure limit """
        def __get__(self):
            return wrapArrhenius(&(self.frxn().low_rate), self)
        def __set__(self, Arrhenius rate):
            self.frxn().low_rate = deref(rate.rate)

    property high_rate:
        """ Get/Set the `Arrhenius` rate constant in the high-pressure limit """
        def __get__(self):
            return wrapArrhenius(&(self.frxn().high_rate), self)
        def __set__(self, Arrhenius rate):
            self.frxn().high_rate = deref(rate.rate)

    property falloff:
        """
        Get/Set the `Falloff` function used to blend the high- and low-pressure
        rate coefficients
        """
        def __get__(self):
            return wrapFalloff(self.frxn().falloff)
        def __set__(self, Falloff f):
            self.frxn().falloff = f._falloff

    property efficiencies:
        """
        Get/Set a `dict` defining non-default third-body efficiencies for this
        reaction, where the keys are the species names and the values are the
        efficiencies.
        """
        def __get__(self):
            return comp_map_to_dict(self.frxn().third_body.efficiencies)
        def __set__(self, eff):
            self.frxn().third_body.efficiencies = comp_map(eff)

    property default_efficiency:
        """
        Get/Set the default third-body efficiency for this reaction, used for
        species used for species not in `efficiencies`.
        """
        def __get__(self):
            return self.frxn().third_body.default_efficiency
        def __set__(self, default_eff):
            self.frxn().third_body.default_efficiency = default_eff

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            cdef CxxFalloffReaction* r = <CxxFalloffReaction*>self.reaction
            return r.allow_negative_pre_exponential_factor
        def __set__(self, allow):
            cdef CxxFalloffReaction* r = <CxxFalloffReaction*>self.reaction
            r.allow_negative_pre_exponential_factor = allow

    def efficiency(self, species):
        """
        Get the efficiency of the third body named *species* considering both
        the default efficiency and species-specific efficiencies.
        """
        return self.frxn().third_body.efficiency(stringify(species))


cdef class ChemicallyActivatedReaction(FalloffReaction):
    """
    A reaction where the rate decreases as pressure increases due to collisional
    stabilization of a reaction intermediate. Like a `FalloffReaction`, except
    that the forward rate constant is written as being proportional to the low-
    pressure rate constant.
    """
    reaction_type = CHEMACT_RXN


cdef class PlogReaction(Reaction):
    """
    A pressure-dependent reaction parameterized by logarithmically interpolating
    between Arrhenius rate expressions at various pressures.
    """
    reaction_type = PLOG_RXN

    property rates:
        """
        Get/Set the rate coefficients for this reaction, which are given as a
        list of (pressure, `Arrhenius`) tuples.
        """
        def __get__(self):
            cdef CxxPlogReaction* r = <CxxPlogReaction*>self.reaction
            rates = []
            cdef vector[pair[double,CxxArrhenius]] cxxrates = r.rate.rates()
            cdef pair[double,CxxArrhenius] p_rate
            for p_rate in cxxrates:
                rates.append((p_rate.first,copyArrhenius(&p_rate.second)))
            return rates

        def __set__(self, rates):
            cdef multimap[double,CxxArrhenius] ratemap
            cdef Arrhenius rate
            cdef pair[double,CxxArrhenius] item
            for p,rate in rates:
                item.first = p
                item.second = deref(rate.rate)
                ratemap.insert(item)

            cdef CxxPlogReaction* r = <CxxPlogReaction*>self.reaction
            r.rate = CxxPlog(ratemap)

    def __call__(self, float T, float P):
        cdef CxxPlogReaction* r = <CxxPlogReaction*>self.reaction
        cdef double logT = np.log(T)
        cdef double recipT = 1/T
        cdef double logP = np.log(P)

        r.rate.update_C(&logP)
        return r.rate.updateRC(logT, recipT)


cdef class ChebyshevReaction(Reaction):
    """
    A pressure-dependent reaction parameterized by a bivariate Chebyshev
    polynomial in temperature and pressure.
    """
    reaction_type = CHEBYSHEV_RXN

    property Tmin:
        """ Minimum temperature [K] for the Chebyshev fit """
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Tmin()

    property Tmax:
        """ Maximum temperature [K] for the Chebyshev fit """
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Tmax()

    property Pmin:
        """ Minimum pressure [Pa] for the Chebyshev fit """
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Pmin()

    property Pmax:
        """ Maximum pressure [Pa] for the Chebyshev fit """
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.Pmax()

    property nPressure:
        """ Number of pressures over which the Chebyshev fit is computed """
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.nPressure()

    property nTemperature:
        """ Number of temperatures over which the Chebyshev fit is computed """
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            return r.rate.nTemperature()

    property coeffs:
        """
        2D array of Chebyshev coefficients of size `(nTemperature, nPressure)`.
        """
        def __get__(self):
            cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
            c = np.fromiter(r.rate.coeffs(), np.double)
            return c.reshape((r.rate.nTemperature(), r.rate.nPressure()))

    def set_parameters(self, Tmin, Tmax, Pmin, Pmax, coeffs):
        """
        Simultaneously set values for `Tmin`, `Tmax`, `Pmin`, `Pmax`, and
        `coeffs`.
        """
        cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction

        cdef CxxArray2D data
        data.resize(len(coeffs), len(coeffs[0]))
        cdef double value
        cdef int i
        cdef int j
        for i,row in enumerate(coeffs):
            for j,value in enumerate(row):
                CxxArray2D_set(data, i, j, value)

        r.rate = CxxChebyshevRate(Tmin, Tmax, Pmin, Pmax, data)

    def __call__(self, float T, float P):
        cdef CxxChebyshevReaction* r = <CxxChebyshevReaction*>self.reaction
        cdef double logT = np.log(T)
        cdef double recipT = 1/T
        cdef double logP = np.log10(P)

        r.rate.update_C(&logP)
        return r.rate.updateRC(logT, recipT)


cdef class InterfaceReaction(ElementaryReaction):
    """ A reaction occurring on an `Interface` (i.e. a surface or an edge) """
    reaction_type = INTERFACE_RXN

    property coverage_deps:
        """
        Get/Set a dict containing adjustments to the Arrhenius rate expression
        dependent on surface species coverages. The keys of the dict are species
        names, and the values are tuples specifying the three coverage
        parameters ``(a, m, E)`` which are the modifiers for the pre-exponential
        factor [m, kmol, s units], the temperature exponent [nondimensional],
        and the activation energy [J/kmol], respectively.
        """
        def __get__(self):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            deps = {}
            cdef pair[string,CxxCoverageDependency] item
            for item in r.coverage_deps:
                deps[pystr(item.first)] = (item.second.a, item.second.m,
                                           item.second.E * gas_constant)
            return deps
        def __set__(self, deps):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            r.coverage_deps.clear()
            cdef str species
            for species, D in deps.items():
                r.coverage_deps[stringify(species)] = CxxCoverageDependency(
                    D[0], D[2] / gas_constant, D[1])

    property is_sticking_coefficient:
        """
        Get/Set a boolean indicating if the rate coefficient for this reaction
        is expressed as a sticking coefficient rather than the forward rate
        constant.
        """
        def __get__(self):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            return r.is_sticking_coefficient
        def __set__(self, stick):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            r.is_sticking_coefficient = stick

    property use_motz_wise_correction:
        """
        Get/Set a boolean indicating whether to use the correction factor
        developed by Motz & Wise for reactions with high (near-unity) sticking
        coefficients when converting the sticking coefficient to a rate
        coefficient.
        """
        def __get__(self):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            return r.use_motz_wise_correction
        def __set__(self, mw):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            r.use_motz_wise_correction = mw

    property sticking_species:
        """
        The name of the sticking species. Needed only for reactions with
        multiple non-surface reactant species, where the sticking species is
        ambiguous.
        """
        def __get__(self):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            return pystr(r.sticking_species)
        def __set__(self, species):
            cdef CxxInterfaceReaction* r = <CxxInterfaceReaction*>self.reaction
            r.sticking_species = stringify(species)


cdef Reaction wrapReaction(shared_ptr[CxxReaction] reaction):
    """
    Wrap a C++ Reaction object with a Python object of the correct derived type.
    """
    cdef int reaction_type = reaction.get().reaction_type

    if reaction_type == ELEMENTARY_RXN:
        R = ElementaryReaction(init=False)
    elif reaction_type == THREE_BODY_RXN:
        R = ThreeBodyReaction(init=False)
    elif reaction_type == FALLOFF_RXN:
        R = FalloffReaction(init=False)
    elif reaction_type == CHEMACT_RXN:
        R = ChemicallyActivatedReaction(init=False)
    elif reaction_type == PLOG_RXN:
        R = PlogReaction(init=False)
    elif reaction_type == CHEBYSHEV_RXN:
        R = ChebyshevReaction(init=False)
    elif reaction_type == INTERFACE_RXN:
        R = InterfaceReaction(init=False)
    else:
        R = Reaction(init=False)

    R._assign(reaction)
    return R

cdef CxxReaction* newReaction(int reaction_type):
    """
    Create a new C++ Reaction object of the specified type
    """
    if reaction_type == ELEMENTARY_RXN:
        return new CxxElementaryReaction()
    elif reaction_type == THREE_BODY_RXN:
        return new CxxThreeBodyReaction()
    elif reaction_type == FALLOFF_RXN:
        return new CxxFalloffReaction()
    elif reaction_type == CHEMACT_RXN:
        return new CxxChemicallyActivatedReaction()
    elif reaction_type == PLOG_RXN:
        return new CxxPlogReaction()
    elif reaction_type == CHEBYSHEV_RXN:
        return new CxxChebyshevReaction()
    elif reaction_type == INTERFACE_RXN:
        return new CxxInterfaceReaction()
    else:
        return new CxxReaction(0)
