# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.


# dictionary to store reaction classes
cdef dict _reaction_class_registry = {}


cdef class _ReactionRate:
    """
    Base class for ReactionRate objects.

    ReactionRate objects are used to calculate reaction rates and are associated
    with a Reaction object.
    """

    def __repr__(self):
        return f"<{pystr(self.base.type())} at {id(self):0x}>"

    property _linked:
        """
        Check whether reaction rate was defined as stand-alone or is linked to
        a Kinetics object.
        """
        def __get__(self):
            return self.base.linked()

    def _release(self):
        """
        Break link between ReactionRate and Kinetics objects
        (for testing purposes)
        """
        self.base.releaseEvaluator()

    def __call__(self, double temperature, pressure=None):
        """
        Evaluate rate expression based on temperature and pressure. For rate
        expressions that are dependent of pressure, an omission of pressure
        will raise an exception.
        """
        if pressure:
            self.base.update(temperature, pressure)
            return self.base.eval(temperature, pressure)
        else:
            self.base.update(temperature)
            return self.base.eval(temperature)

    def ddT(self, double temperature, pressure=None):
        """
        Evaluate derivative of rate expression with respect to temperature.
        """
        if pressure:
            self.base.update(temperature, pressure)
            return self.base.ddT(temperature, pressure)
        else:
            self.base.update(temperature)
            return self.base.ddT(temperature)

    property input_data:
        """
        Get input data for this reaction rate with its current parameter values.
        """
        def __get__(self):
            return anymap_to_dict(self.base.parameters())


cdef class ArrheniusRate(_ReactionRate):
    r"""
    A reaction rate coefficient which depends on temperature only and follows
    the modified Arrhenius form:

    .. math::

        k_f = A T^b \exp{-\tfrac{E}{RT}}

    where *A* is the `pre_exponential_factor`, *b* is the `temperature_exponent`,
    and *Ea* is the `activation_energy`.
    """
    def __cinit__(self, A=None, b=None, Ea=None, input_data=None, init=True):

        if init:
            if isinstance(input_data, dict):
                self._base.reset(new CxxArrheniusRate(dict_to_anymap(input_data)))
            elif all([arg is not None for arg in [A, b, Ea]]):
                self._base.reset(new CxxArrheniusRate(A, b, Ea))
            elif all([arg is None for arg in [A, b, Ea, input_data]]):
                self._base.reset(new CxxArrheniusRate(dict_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameters 'A', 'b' or 'Ea'")
            self.base = self._base.get()
            self.rate = <CxxArrheniusRate*>(self.base)

    @staticmethod
    cdef wrap(shared_ptr[CxxReactionRateBase] rate):
        """
        Wrap a C++ ReactionRateBase object with a Python object.
        """
        # wrap C++ reaction
        cdef ArrheniusRate arr
        arr = ArrheniusRate(init=False)
        arr._base = rate
        arr.base = arr._base.get()
        arr.rate = <CxxArrheniusRate*>(arr.base)
        return arr

    property pre_exponential_factor:
        """
        The pre-exponential factor *A* in units of m, kmol, and s raised to
        powers depending on the reaction order.
        """
        def __set__(self, float A):
            self.rate.setPreExponentialFactor(A)
        def __get__(self):
            return self.rate.preExponentialFactor()

    property temperature_exponent:
        """
        The temperature exponent *b*.
        """
        def __set__(self, float b):
            self.rate.setTemperatureExponent(b)
        def __get__(self):
            return self.rate.temperatureExponent()

    property activation_energy:
        """
        The activation energy *Ea* [J/kmol].
        """
        def __set__(self, float Ea):
            self.rate.setActivationEnergy(Ea)
        def __get__(self):
            return self.rate.activationEnergy()

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            return self.rate.allow_negative_pre_exponential_factor
        def __set__(self, allow):
            self.rate.allow_negative_pre_exponential_factor = allow


cdef class PlogRate(_ReactionRate):
    r"""
    A pressure-dependent reaction rate parameterized by logarithmically
    interpolating between Arrhenius rate expressions at various pressures.
    """
    def __cinit__(self, rates=None, input_data=None, init=True):

        if init and isinstance(rates, list):
            self._base.reset(new CxxPlogRate(self._cxxrates(rates)))

        elif init:
            if isinstance(input_data, dict):
                self._base.reset(new CxxPlogRate(dict_to_anymap(input_data)))
            elif rates is None:
                self._base.reset(new CxxPlogRate(dict_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameter 'rates'")

        self.base = self._base.get()
        self.rate = <CxxPlogRate*>(self.base)

    @staticmethod
    cdef wrap(shared_ptr[CxxReactionRateBase] rate):
        """
        Wrap a C++ ReactionRateBase object with a Python object.
        """
        # wrap C++ reaction
        cdef PlogRate arr
        arr = PlogRate(init=False)
        arr._base = rate
        arr.base = arr._base.get()
        arr.rate = <CxxPlogRate*>(arr.base)
        return arr

    cdef vector[pair[double, CxxArrhenius]] _cxxrates(self, rates):
        cdef vector[pair[double, CxxArrhenius]] cxxrates
        cdef Arrhenius rate
        cdef pair[double, CxxArrhenius] item
        for p, rate in rates:
            item.first = p
            item.second = deref(rate.rate)
            cxxrates.push_back(item)
        return cxxrates

    property rates:
        """
        Get/Set the rate coefficients for this reaction, which are given as a
        list of (pressure, `Arrhenius`) tuples.
        """
        def __get__(self):
            rates = []
            cdef vector[pair[double, CxxArrhenius]] cxxrates = self.rate.rates()
            cdef pair[double, CxxArrhenius] p_rate
            for p_rate in cxxrates:
                rates.append((p_rate.first, copyArrhenius(&p_rate.second)))
            return rates

        def __set__(self, rates):
            self.rate.setRates(self._cxxrates(rates))


cdef class ChebyshevRate(_ReactionRate):
    r"""
    A pressure-dependent reaction rate parameterized by a bivariate Chebyshev
    polynomial in temperature and pressure.
    """
    def __cinit__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, data=None,
                  input_data=None, init=True):

        if init:
            if isinstance(input_data, dict):
                self._base.reset(new CxxChebyshevRate3(dict_to_anymap(input_data)))
            elif all([arg is not None for arg in [Tmin, Tmax, Pmin, Pmax, data]]):
                self._base.reset(new CxxChebyshevRate3(Tmin, Tmax, Pmin, Pmax,
                                                       self._cxxarray2d(data)))
            elif all([arg is None
                    for arg in [Tmin, Tmax, Pmin, Pmax, data, input_data]]):
                self._base.reset(new CxxChebyshevRate3(dict_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameters")
            self.base = self._base.get()
            self.rate = <CxxChebyshevRate3*>(self.base)

    cdef CxxArray2D _cxxarray2d(self, coeffs):
        cdef CxxArray2D data
        if isinstance(coeffs, np.ndarray):
            coeffs = coeffs.tolist()
        data.resize(len(coeffs), len(coeffs[0]))
        cdef double value
        cdef int i
        cdef int j
        for i,row in enumerate(coeffs):
            for j,value in enumerate(row):
                CxxArray2D_set(data, i, j, value)
        return data

    @staticmethod
    cdef wrap(shared_ptr[CxxReactionRateBase] rate):
        """
        Wrap a C++ ReactionRateBase object with a Python object.
        """
        # wrap C++ reaction
        cdef ChebyshevRate arr
        arr = ChebyshevRate(init=False)
        arr._base = rate
        arr.base = arr._base.get()
        arr.rate = <CxxChebyshevRate3*>(arr.base)
        return arr

    property Tmin:
        """ Minimum temperature [K] for the Chebyshev fit """
        def __get__(self):
            return self.rate.Tmin()

    property Tmax:
        """ Maximum temperature [K] for the Chebyshev fit """
        def __get__(self):
            return self.rate.Tmax()

    property Pmin:
        """ Minimum pressure [Pa] for the Chebyshev fit """
        def __get__(self):
            return self.rate.Pmin()

    property Pmax:
        """ Maximum pressure [Pa] for the Chebyshev fit """
        def __get__(self):
            return self.rate.Pmax()

    property n_pressure:
        """ Number of pressures over which the Chebyshev fit is computed """
        def __get__(self):
            return self.rate.nPressure()

    property n_temperature:
        """ Number of temperatures over which the Chebyshev fit is computed """
        def __get__(self):
            return self.rate.nTemperature()

    property coeffs:
        """
        2D array of Chebyshev coefficients of size `(n_temperature, n_pressure)`.
        """
        def __get__(self):
            cdef CxxArray2D cxxcoeffs = self.rate.coeffs()
            c = np.fromiter(cxxcoeffs.data(), np.double)
            return c.reshape(cxxcoeffs.nRows(), cxxcoeffs.nColumns(), order="F")

        def __set__(self, coeffs):
            self.rate.setCoeffs(self._cxxarray2d(coeffs))


cdef class CustomRate(_ReactionRate):
    r"""
    A custom rate coefficient which depends on temperature only.

    The simplest way to create a `CustomRate` object is to use a lambda function,
    for example::

        rr = CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15/T))
    """
    def __cinit__(self, k=None, init=True):

        if init:
            self._base.reset(new CxxCustomFunc1Rate())
            self.base = self._base.get()
            self.rate = <CxxCustomFunc1Rate*>(self.base)
            self.set_rate_function(k)

    def set_rate_function(self, k):
        r"""
        Set the function describing a custom reaction rate::

            rr = CustomRate()
            rr.set_rate_function(lambda T: 38.7 * T**2.7 * exp(-3150.15/T))
        """
        if k is None:
            self._rate_func = None
            return

        if isinstance(k, Func1):
            self._rate_func = k
        else:
            self._rate_func = Func1(k)

        self.rate.setRateFunction(self._rate_func._func)


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

        R = ct.Reaction.listFromFile("gri30.yaml", gas)
        R = ct.Reaction.listFromCti(open("path/to/gri30.cti").read())
        R = ct.Reaction.listFromXml(open("path/to/gri30.xml").read())

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
    _reaction_type = ""
    _has_legacy = False
    _hybrid = False # indicate whether legacy implementations are separate or merged

    def __cinit__(self, reactants="", products="", init=True, legacy=False, **kwargs):

        if init:
            rxn_type = self._reaction_type
            if (not self._hybrid and self._has_legacy) or (self._hybrid and legacy):
                rxn_type += "-legacy"
            self._reaction = CxxNewReaction(stringify((rxn_type)))
            self.reaction = self._reaction.get()
            if reactants:
                self.reactants = reactants
            if products:
                self.products = products

    @staticmethod
    cdef wrap(shared_ptr[CxxReaction] reaction):
        """
        Wrap a C++ Reaction object with a Python object of the correct derived type.
        """
        # ensure all reaction types are registered
        if not(_reaction_class_registry):
            def register_subclasses(cls):
                for c in cls.__subclasses__():
                    rxn_type = getattr(c, "_reaction_type")
                    if getattr(c, "_hybrid"):
                        # registry needs to contain both updated and "-legacy" variants
                        _reaction_class_registry[rxn_type] = c
                    if getattr(c, "_has_legacy", False):
                        rxn_type += "-legacy"
                    _reaction_class_registry[rxn_type] = c
                    register_subclasses(c)

            # update global reaction class registry
            register_subclasses(Reaction)

        # identify class
        rxn_type = pystr(reaction.get().type())
        cls = _reaction_class_registry.get(rxn_type, Reaction)

        # wrap C++ reaction
        cdef Reaction R
        R = cls(init=False)
        R._reaction = reaction
        R.reaction = R._reaction.get()
        return R

    @staticmethod
    def fromCti(text):
        """
        Create a Reaction object from its CTI string representation.

        .. deprecated:: 2.5

            The CTI input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_reactions = CxxGetReactions(deref(CxxGetXmlFromString(stringify(text))))
        assert cxx_reactions.size() == 1, cxx_reactions.size()
        return Reaction.wrap(cxx_reactions[0])

    @staticmethod
    def fromXml(text):
        """
        Create a Reaction object from its XML string representation.

        .. deprecated:: 2.5

            The XML input format is deprecated and will be removed in Cantera 3.0.
        """
        cxx_reaction = CxxNewReaction(deref(CxxGetXmlFromString(stringify(text))))
        return Reaction.wrap(cxx_reaction)

    @classmethod
    def from_dict(cls, data, Kinetics kinetics=None):
        """
        Create a `Reaction` object from a dictionary corresponding to its YAML
        representation.

        An example for the creation of a Reaction from a dictionary is::

            rxn = Reaction.from_dict(
                {"equation": "O + H2 <=> H + OH",
                 "rate-constant": {"A": 38.7, "b": 2.7, "Ea": 26191840.0}},
                kinetics=gas)

        In the example, *gas* is a Kinetics (or Solution) object.

        :param data:
            A dictionary corresponding to the YAML representation.
        :param kinetics:
            A `Kinetics` object whose associated phase(s) contain the species
            involved in the reaction.
        """
        if cls._reaction_type != "":
            raise TypeError(
                "Class method 'from_dict' was invoked from '{}' but should "
                "be called from base class 'Reaction'".format(cls.__name__))
        if kinetics is None:
            raise ValueError("A Kinetics object is required.")

        cdef CxxAnyMap any_map = dict_to_anymap(data)
        cxx_reaction = CxxNewReaction(any_map, deref(kinetics.kinetics))
        return Reaction.wrap(cxx_reaction)

    @classmethod
    def fromYaml(cls, text, Kinetics kinetics=None):
        """
        Create a `Reaction` object from its YAML string representation.

        An example for the creation of a Reaction from a YAML string is::

            rxn = Reaction.fromYaml('''
                equation: O + H2 <=> H + OH
                rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
                ''', kinetics=gas)

        In the example, *gas* is a Kinetics (or Solution) object.

        :param text:
            The YAML reaction string.
        :param kinetics:
            A `Kinetics` object whose associated phase(s) contain the species
            involved in the reaction.
        """
        if cls._reaction_type != "":
            raise TypeError(
                "Class method 'fromYaml' was invoked from '{}' but should "
                "be called from base class 'Reaction'".format(cls.__name__))
        if kinetics is None:
            raise ValueError("A Kinetics object is required.")

        cdef CxxAnyMap any_map
        any_map = AnyMapFromYamlString(stringify(text))
        cxx_reaction = CxxNewReaction(any_map, deref(kinetics.kinetics))
        return Reaction.wrap(cxx_reaction)

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
        return [Reaction.wrap(r) for r in cxx_reactions]

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
        return [Reaction.wrap(r) for r in cxx_reactions]

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
        return [Reaction.wrap(r) for r in cxx_reactions]

    @staticmethod
    def listFromYaml(text, Kinetics kinetics):
        """
        Create a list of `Reaction` objects from all the reactions defined in a
        YAML string.
        """
        root = AnyMapFromYamlString(stringify(text))
        cxx_reactions = CxxGetReactions(root[stringify("items")],
                                        deref(kinetics.kinetics))
        return [Reaction.wrap(r) for r in cxx_reactions]

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

    property reaction_type:
        """
        Retrieve the native type name of the reaction.
        """
        def __get__(self):
            return pystr(self.reaction.type())

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

    property input_data:
        """
        Get input data for this reaction with its current parameter values,
        along with any user-specified data provided with its input (YAML)
        definition.
        """
        def __get__(self):
            return anymap_to_dict(self.reaction.parameters(True))

    def update_user_data(self, data):
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.reaction.input.update(dict_to_anymap(data), False)

    def clear_user_data(self):
        """
        Clear all saved input data, so that the data given by `input_data` or
        `Solution.write_yaml` will only include values generated by Cantera based on
        the current object state.
        """
        self.reaction.input.clear()

    def __repr__(self):
        return f"<{self.__class__.__name__}: {self.equation}>"

    def __str__(self):
        return self.equation

    def _deprecation_warning(self, attr, what="property"):
        return ("\n{} '{}' to be removed after Cantera 2.6.\nThis {} is moved to "
                "the {} object accessed via the 'rate' property."
                ).format(what.capitalize(), attr, what, type(self.rate).__name__)

    property uses_legacy:
        """Indicate whether reaction uses a legacy implementation"""
        def __get__(self):
            return self.reaction.usesLegacy()

    property _linked:
        """
        Check whether reaction was defined as stand-alone or is linked to
        a Kinetics object.
        """
        def __get__(self):
            return self.reaction.linked()

    property index:
        """
        Index of reaction within the Kinetics object (if applicable). Raises
        an exception if the reaction is not linked.
        """
        def __get__(self):
            return self.reaction.index()


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
            self.own_rate = True
            self.reaction = None
        else:
            self.own_rate = False

    def __dealloc__(self):
        if self.own_rate:
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

    An example for the definition of an `ElementaryReaction` object is given as::

        rxn = ElementaryReaction(
            equation="O + H2 <=> H + OH",
            rate={"A": 38.7, "b": 2.7, "Ea": 2.619184e+07},
            kinetics=gas)

    The YAML description corresponding to this reaction is::

        equation: O + H2 <=> H + OH
        rate-constant: {A: 3.87e+04 cm^3/mol/s, b: 2.7, Ea: 6260.0 cal/mol}
    """
    _reaction_type = "elementary"
    _has_legacy = True
    _hybrid = True

    cdef CxxElementaryReaction3* er(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxElementaryReaction3*>self.reaction

    cdef CxxElementaryReaction2* er2(self):
        if not self.uses_legacy:
            raise AttributeError("Incorrect accessor for legacy implementation")
        return <CxxElementaryReaction2*>self.reaction

    def __init__(self, equation=None, rate=None, Kinetics kinetics=None,
                 init=True, legacy=False, **kwargs):

        if init and equation and kinetics:

            rxn_type = self._reaction_type
            if legacy:
                rxn_type += "-legacy"
            spec = {"equation": equation, "type": rxn_type}
            if isinstance(rate, dict):
                spec["rate-constant"] = rate
            elif legacy and (isinstance(rate, Arrhenius) or rate is None):
                spec["rate-constant"] = dict.fromkeys(["A", "b", "Ea"], 0.)
            elif rate is None:
                pass
            elif not legacy and isinstance(rate, (Arrhenius, ArrheniusRate)):
                pass
            else:
                raise TypeError("Invalid rate definition")

            self._reaction = CxxNewReaction(dict_to_anymap(spec),
                                            deref(kinetics.kinetics))
            self.reaction = self._reaction.get()

            if legacy and isinstance(rate, Arrhenius):
                self.rate = rate
            elif not legacy and isinstance(rate, (ArrheniusRate, Arrhenius)):
                self.rate = rate

    cdef _legacy_set_rate(self, Arrhenius rate):
        cdef CxxElementaryReaction2* r = self.er2()
        r.rate = deref(rate.rate)

    property rate:
        """ Get/Set the `ArrheniusRate` rate coefficient for this reaction. """
        def __get__(self):
            if self.uses_legacy:
                return wrapArrhenius(&(self.er2().rate), self)

            return ArrheniusRate.wrap(self.er().rate())
        def __set__(self, rate):
            if self.uses_legacy:
                self._legacy_set_rate(rate)
                return

            cdef ArrheniusRate rate3
            if isinstance(rate, ArrheniusRate):
                rate3 = rate
            elif isinstance(rate, Arrhenius):
                warnings.warn("Setting the rate using an 'Arrhenius' object is "
                    "deprecated and will be removed after Cantera 2.6. The argument "
                    "type is replaceable by 'ArrheniusRate'.", DeprecationWarning)
                rate3 = ArrheniusRate(rate.pre_exponential_factor,
                                      rate.temperature_exponent,
                                      rate.activation_energy)
            else:
                raise TypeError("Invalid rate definition")
            self.er().setRate(rate3._base)

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ArrheniusRate.allow_negative_pre_exponential_factor`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.er2().allow_negative_pre_exponential_factor

            attr = "allow_negative_pre_exponential_factor"
            warnings.warn(self._deprecation_warning(attr), DeprecationWarning)
            return self.rate.allow_negative_pre_exponential_factor
        def __set__(self, allow):
            if self.uses_legacy:
                self.er2().allow_negative_pre_exponential_factor = allow
                return

            attr = "allow_negative_pre_exponential_factor"
            warnings.warn(self._deprecation_warning(attr), DeprecationWarning)
            self.rate.allow_negative_pre_exponential_factor = allow


cdef class ThreeBodyReaction(ElementaryReaction):
    """
    A reaction with a non-reacting third body "M" that acts to add or remove
    energy from the reacting species.

    An example for the definition of an `ThreeBodyReaction` object is given as::

        rxn = ThreeBodyReaction(
            equation="2 O + M <=> O2 + M",
            type="three-body",
            rate={"A": 1.2e+17, "b": -1.0, "Ea": 0.0},
            efficiencies={"H2": 2.4, "H2O": 15.4, "AR": 0.83},
            kinetics=gas)

    The YAML description corresponding to this reaction is::

        equation: 2 O + M <=> O2 + M
        type: three-body
        rate-constant: {A: 1.2e+17 cm^6/mol^2/s, b: -1.0, Ea: 0.0 cal/mol}
        efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}
    """
    _reaction_type = "three-body"
    _has_legacy = True
    _hybrid = True

    cdef CxxThreeBodyReaction3* tbr(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxThreeBodyReaction3*>self.reaction

    cdef CxxThreeBodyReaction2* tbr2(self):
        if not self.uses_legacy:
            raise AttributeError("Incorrect accessor for legacy implementation")
        return <CxxThreeBodyReaction2*>self.reaction

    cdef CxxThirdBody* thirdbody(self):
        if self.uses_legacy:
            return &(self.tbr2().third_body)
        return <CxxThirdBody*>(self.tbr().thirdBody().get())

    def __init__(self, equation=None, rate=None, efficiencies=None,
                 Kinetics kinetics=None, legacy=False, init=True, **kwargs):

        if init and equation and kinetics:

            rxn_type = self._reaction_type
            if legacy:
                rxn_type += "-legacy"
            spec = {"equation": equation, "type": rxn_type}
            if isinstance(rate, dict):
                spec["rate-constant"] = rate
            elif legacy and (isinstance(rate, Arrhenius) or rate is None):
                spec["rate-constant"] = dict.fromkeys(["A", "b", "Ea"], 0.)
            elif rate is None:
                pass
            elif not legacy and isinstance(rate, (Arrhenius, ArrheniusRate)):
                pass
            else:
                raise TypeError("Invalid rate definition")

            if isinstance(efficiencies, dict):
                spec["efficiencies"] = efficiencies

            self._reaction = CxxNewReaction(dict_to_anymap(spec),
                                            deref(kinetics.kinetics))
            self.reaction = self._reaction.get()

            if legacy and isinstance(rate, Arrhenius):
                self.rate = rate
            elif not legacy and isinstance(rate, (Arrhenius, ArrheniusRate)):
                self.rate = rate

    property efficiencies:
        """
        Get/Set a `dict` defining non-default third-body efficiencies for this
        reaction, where the keys are the species names and the values are the
        efficiencies.
        """
        def __get__(self):
            return comp_map_to_dict(self.thirdbody().efficiencies)
        def __set__(self, eff):
            self.thirdbody().efficiencies = comp_map(eff)

    property default_efficiency:
        """
        Get/Set the default third-body efficiency for this reaction, used for
        species used for species not in `efficiencies`.
        """
        def __get__(self):
            return self.thirdbody().default_efficiency
        def __set__(self, default_eff):
            self.thirdbody().default_efficiency = default_eff

    def efficiency(self, species):
        """
        Get the efficiency of the third body named *species* considering both
        the default efficiency and species-specific efficiencies.
        """
        return self.thirdbody().efficiency(stringify(species))


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
    _reaction_type = "falloff"

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
    _reaction_type = "chemically-activated"


cdef class PlogReaction(Reaction):
    """
    A pressure-dependent reaction parameterized by logarithmically interpolating
    between Arrhenius rate expressions at various pressures.

    An example for the definition of a `PlogReaction` object is given as::

        rxn = PlogReaction(
            equation="H2 + O2 <=> 2 OH",
            rate=[(1013.25, Arrhenius(1.2124e+16, -0.5779, 45491376.8)),
                  (101325., Arrhenius(4.9108e+31, -4.8507, 103649395.2)),
                  (1013250., Arrhenius(1.2866e+47, -9.0246, 166508556.0)),
                  (10132500., Arrhenius(5.9632e+56, -11.529, 220076726.4))],
            kinetics=gas)

    The YAML description corresponding to this reaction is::

        equation: H2 + O2 <=> 2 OH
        type: pressure-dependent-Arrhenius
        rate-constants:
        - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04 cal/mol}
        - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04 cal/mol}
        - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04 cal/mol}
        - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04 cal/mol}or.
    """
    _reaction_type = "pressure-dependent-Arrhenius"
    _has_legacy = True
    _hybrid = True

    cdef CxxPlogReaction3* pr(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxPlogReaction3*>self.reaction

    cdef CxxPlogReaction2* cp2(self):
        if not self.uses_legacy:
            raise AttributeError("Incorrect accessor for legacy implementation")
        return <CxxPlogReaction2*>self.reaction

    def __init__(self, equation=None, rate=None, Kinetics kinetics=None,
                 init=True, legacy=False, **kwargs):

        if init and equation and kinetics:

            rxn_type = self._reaction_type
            if legacy:
                rxn_type += "-legacy"
            spec = {"equation": equation, "type": rxn_type}
            if legacy and isinstance(rate, list):
                rates = []
                for r in rate:
                    rates.append({
                        "P": r[0],
                        "A": r[1].pre_exponential_factor,
                        "b": r[1].temperature_exponent,
                        "Ea": r[1].activation_energy,
                    })
                spec.update({'rate-constants': rates})
            elif not legacy and (isinstance(rate, (PlogRate, list)) or rate is None):
                pass
            else:
                raise TypeError("Invalid rate definition")

            self._reaction = CxxNewReaction(dict_to_anymap(spec),
                                            deref(kinetics.kinetics))
            self.reaction = self._reaction.get()

            if not legacy and isinstance(rate, PlogRate):
                self.rate = rate
            elif not legacy and isinstance(rate, list):
                self.rate = PlogRate(rate)

    property rate:
        """ Get/Set the `PlogRate` rate coefficients for this reaction. """
        def __get__(self):
            if self.uses_legacy:
                raise AttributeError("Legacy implementation does not use rate property.")
            return PlogRate.wrap(self.pr().rate())
        def __set__(self, PlogRate rate):
            if self.uses_legacy:
                raise AttributeError("Legacy implementation does not use rate property.")
            self.pr().setRate(rate._base)

    cdef list _legacy_get_rates(self):
        cdef CxxPlogReaction2* r = self.cp2()
        cdef vector[pair[double,CxxArrhenius]] cxxrates = r.rate.rates()
        cdef pair[double,CxxArrhenius] p_rate
        rates = []
        for p_rate in cxxrates:
            rates.append((p_rate.first,copyArrhenius(&p_rate.second)))
        return rates

    cdef _legacy_set_rates(self, list rates):
        cdef multimap[double,CxxArrhenius] ratemap
        cdef Arrhenius rate
        cdef pair[double,CxxArrhenius] item
        for p, rate in rates:
            item.first = p
            item.second = deref(rate.rate)
            ratemap.insert(item)

        cdef CxxPlogReaction2* r = self.cp2()
        r.rate = CxxPlog(ratemap)

    property rates:
        """
        Get/Set the rate coefficients for this reaction, which are given as a
        list of (pressure, `Arrhenius`) tuples.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `PlogRate.rates`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self._legacy_get_rates()

            warnings.warn(self._deprecation_warning("rates"), DeprecationWarning)
            return self.rate.rates

        def __set__(self, rates):
            if self.uses_legacy:
                self._legacy_set_rates(rates)
                return

            warnings.warn("Property 'rates' to be removed after Cantera 2.6. "
                "Setter is replaceable by assigning a new 'PlogRate' object created "
                "from rates to the rate property.", DeprecationWarning)
            rate_ = self.rate
            rate_.rates = rates
            self.rate = rate_

    cdef _legacy_call(self, float T, float P):
        cdef CxxPlogReaction2* r = self.cp2()
        cdef double logT = np.log(T)
        cdef double recipT = 1/T
        cdef double logP = np.log(P)

        r.rate.update_C(&logP)
        return r.rate.updateRC(logT, recipT)

    def __call__(self, float T, float P):
        """
        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaceable by call to `rate` property.
        """
        if self.uses_legacy:
            return self._legacy_call(T, P)

        warnings.warn(
            self._deprecation_warning("__call__", "method"), DeprecationWarning)
        return self.rate(T, P)


cdef class ChebyshevReaction(Reaction):
    """
    A pressure-dependent reaction parameterized by a bivariate Chebyshev
    polynomial in temperature and pressure.

    An example for the definition of a `ChebyshevReaction` object is given as::

        rxn = ChebyshevReaction(
            equation="HO2 <=> OH + O",
            rate={"Tmin": 290.0, "Tmax": 3000.0,
                  "Pmin": 1e3, "Pmax": 1e8,
                  "data": [[8.2883, -1.1397, -0.12059, 0.016034],
                           [1.9764, 1.0037, 7.2865e-03, -0.030432],
                           [0.3177, 0.26889, 0.094806, -7.6385e-03]]},
            kinetics=gas)

    The YAML description corresponding to this reaction is::

        equation: HO2 <=> OH + O
        type: Chebyshev
        temperature-range: [290.0, 3000.0]
        pressure-range: [1.e-03 bar, 10. bar]
        data:
        - [8.2883, -1.1397, -0.12059, 0.016034]
        - [1.9764, 1.0037, 7.2865e-03, -0.030432]
        - [0.3177, 0.26889, 0.094806, -7.6385e-03]
    """
    _reaction_type = "Chebyshev"
    _has_legacy = True
    _hybrid = True

    cdef CxxChebyshevReaction3* cr(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxChebyshevReaction3*>self.reaction

    cdef CxxChebyshevReaction2* cr2(self):
        if not self.uses_legacy:
            raise AttributeError("Incorrect accessor for legacy implementation")
        return <CxxChebyshevReaction2*>self.reaction

    def __init__(self, equation=None, rate=None, Kinetics kinetics=None,
                 init=True, legacy=False, **kwargs):

        if init and equation and kinetics:

            rxn_type = self._reaction_type
            if legacy:
                rxn_type += "-legacy"
            spec = {"equation": equation, "type": rxn_type}
            if isinstance(rate, dict):
                spec["temperature-range"] = [rate["Tmin"], rate["Tmax"]]
                spec["pressure-range"] = [rate["Pmin"], rate["Pmax"]]
                spec["data"] = rate["data"]
            elif not legacy and (isinstance(rate, ChebyshevRate) or rate is None):
                pass
            else:
                raise TypeError("Invalid rate definition")

            self._reaction = CxxNewReaction(dict_to_anymap(spec),
                                            deref(kinetics.kinetics))
            self.reaction = self._reaction.get()

            if not legacy and isinstance(rate, ChebyshevRate):
                self.rate = rate

    property rate:
        """ Get/Set the `ChebyshevRate` rate coefficients for this reaction. """
        def __get__(self):
            if self.uses_legacy:
                raise AttributeError("Legacy implementation does not use rate property.")
            return ChebyshevRate.wrap(self.cr().rate())
        def __set__(self, ChebyshevRate rate):
            if self.uses_legacy:
                raise AttributeError("Legacy implementation does not use rate property.")
            self.cr().setRate(rate._base)

    property Tmin:
        """
        Minimum temperature [K] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.Tmin`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cr2().rate.Tmin()

            warnings.warn(self._deprecation_warning("Tmin"), DeprecationWarning)
            return self.rate.Tmin

    property Tmax:
        """
        Maximum temperature [K] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.Tmax`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cr2().rate.Tmax()

            warnings.warn(self._deprecation_warning("Tmax"), DeprecationWarning)
            return self.rate.Tmax

    property Pmin:
        """
        Minimum pressure [Pa] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.Pmin`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cr2().rate.Pmin()

            warnings.warn(self._deprecation_warning("Pmin"), DeprecationWarning)
            return self.rate.Pmin

    property Pmax:
        """ Maximum pressure [Pa] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.Pmax`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cr2().rate.Pmax()

            warnings.warn(self._deprecation_warning("Pmax"), DeprecationWarning)
            return self.rate.Pmax

    property nPressure:
        """
        Number of pressures over which the Chebyshev fit is computed

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.nPressure`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cr2().rate.nPressure()

            warnings.warn(self._deprecation_warning("nPressure"), DeprecationWarning)
            return self.rate.n_pressure

    property nTemperature:
        """
        Number of temperatures over which the Chebyshev fit is computed

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.nTemperature`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cr2().rate.nTemperature()

            warnings.warn(
                self._deprecation_warning("nTemperature"), DeprecationWarning)
            return self.rate.n_temperature

    cdef _legacy_get_coeffs(self):
        cdef CxxChebyshevReaction2* r = self.cr2()
        c = np.fromiter(r.rate.coeffs(), np.double)
        return c.reshape((r.rate.nTemperature(), r.rate.nPressure()))

    property coeffs:
        """
        2D array of Chebyshev coefficients of size `(n_temperature, n_pressure)`.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.coeffs`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self._legacy_get_coeffs()

            warnings.warn(self._deprecation_warning("coeffs"), DeprecationWarning)
            return self.rate.coeffs

    cdef _legacy_set_parameters(self, Tmin, Tmax, Pmin, Pmax, coeffs):
        cdef CxxChebyshevReaction2* r = self.cr2()

        cdef CxxArray2D data
        data.resize(len(coeffs), len(coeffs[0]))
        cdef double value
        cdef int i
        cdef int j
        for i,row in enumerate(coeffs):
            for j, value in enumerate(row):
                CxxArray2D_set(data, i, j, value)

        r.rate = CxxChebyshev(Tmin, Tmax, Pmin, Pmax, data)

    def set_parameters(self, Tmin, Tmax, Pmin, Pmax, coeffs):
        """
        Simultaneously set values for `Tmin`, `Tmax`, `Pmin`, `Pmax`, and
        `coeffs`.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by `ChebyshevRate` constructor.
        """
        if self.uses_legacy:
            return self._legacy_set_parameters(Tmin, Tmax, Pmin, Pmax, coeffs)

        warnings.warn("Method 'set_parameters' to be removed after Cantera 2.6. "
            "Method is replaceable by assigning a new 'ChebyshevRate' object to the "
            "rate property.", DeprecationWarning)
        self.rate = ChebyshevRate(Tmin, Tmax, Pmin, Pmax, coeffs)

    cdef _legacy_call(self, float T, float P):
        cdef CxxChebyshevReaction2* r = self.cr2()
        cdef double logT = np.log(T)
        cdef double recipT = 1/T
        cdef double logP = np.log10(P)

        r.rate.update_C(&logP)
        return r.rate.updateRC(logT, recipT)

    def __call__(self, float T, float P):
        """
        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaceable by call to `rate` property.
        """
        if self.uses_legacy:
            return self._legacy_call(T, P)

        warnings.warn(
            self._deprecation_warning("__call__", "method"), DeprecationWarning)
        return self.rate(T, P)


cdef class BlowersMasel:
    """
    A reaction rate coefficient which depends on temperature and enthalpy change
    of the reaction follows Blowers-Masel approximation and modified Arrhenius form
    described in `Arrhenius`. The functions are used by reactions defined using
    `BlowersMaselReaction` and `BlowersMaselInterfaceReaction`.
    """
    def __cinit__(self, A=0, b=0, E0=0, w=0, init=True):
        if init:
            self.rate = new CxxBlowersMasel(A, b, E0 / gas_constant, w / gas_constant)
            self.own_rate = True
            self.reaction = None
        else:
            self.own_rate = False

    def __dealloc__(self):
        if self.own_rate:
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

    def activation_energy(self, float dH):
        """
        The activation energy *E* [J/kmol]

        :param dH: The enthalpy of reaction [J/kmol] at the current temperature
        """
        return self.rate.activationEnergy_R(dH) * gas_constant

    property bond_energy:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction *E* [J/kmol].
        """
        def __get__(self):
            return self.rate.bondEnergy() * gas_constant

    property intrinsic_activation_energy:
        """
        The intrinsic activation energy *E0* [J/kmol].
        """
        def __get__(self):
            return self.rate.activationEnergy_R0() * gas_constant

    def __repr__(self):
        return 'BlowersMasel(A={:g}, b={:g}, E0={:g}, w={:g})'.format(
            self.pre_exponential_factor, self.temperature_exponent,
            self.intrinsic_activation_energy, self.bond_energy)

    def __call__(self, float T, float dH):
        cdef double logT = np.log(T)
        cdef double recipT = 1/T
        return self.rate.updateRC(logT, recipT, dH)


cdef wrapBlowersMasel(CxxBlowersMasel* rate, Reaction reaction):
    r = BlowersMasel(init=False)
    r.rate = rate
    r.reaction = reaction
    return r


cdef class BlowersMaselReaction(Reaction):
    """
    A reaction which follows mass-action kinetics with Blowers Masel
    reaction rate.
    """
    _reaction_type = "Blowers-Masel"

    property rate:
        """ Get/Set the `Arrhenius` rate coefficient for this reaction. """
        def __get__(self):
            cdef CxxBlowersMaselReaction* r = <CxxBlowersMaselReaction*>self.reaction
            return wrapBlowersMasel(&(r.rate), self)
        def __set__(self, BlowersMasel rate):
            cdef CxxBlowersMaselReaction* r = <CxxBlowersMaselReaction*>self.reaction
            r.rate = deref(rate.rate)

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            cdef CxxBlowersMaselReaction* r = <CxxBlowersMaselReaction*>self.reaction
            return r.allow_negative_pre_exponential_factor
        def __set__(self, allow):
            cdef CxxBlowersMaselReaction* r = <CxxBlowersMaselReaction*>self.reaction
            r.allow_negative_pre_exponential_factor = allow


cdef class InterfaceReaction(ElementaryReaction):
    """ A reaction occurring on an `Interface` (i.e. a surface or an edge) """
    _reaction_type = "interface"
    _has_legacy = False
    _hybrid = False

    property uses_legacy:
        # legacy framework is implicitly used
        def __get__(self):
            return True

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


cdef class BlowersMaselInterfaceReaction(BlowersMaselReaction):
    """
    A reaction occurring on an `Interface` (i.e. a surface or an edge)
    with the rate parameterization of `BlowersMasel`.
    """
    _reaction_type = "surface-Blowers-Masel"

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
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
            deps = {}
            cdef pair[string,CxxCoverageDependency] item
            for item in r.coverage_deps:
                deps[pystr(item.first)] = (item.second.a, item.second.m,
                                           item.second.E * gas_constant)
            return deps
        def __set__(self, deps):
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
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
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
            return r.is_sticking_coefficient
        def __set__(self, stick):
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
            r.is_sticking_coefficient = stick

    property use_motz_wise_correction:
        """
        Get/Set a boolean indicating whether to use the correction factor
        developed by Motz & Wise for reactions with high (near-unity) sticking
        coefficients when converting the sticking coefficient to a rate
        coefficient.
        """
        def __get__(self):
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
            return r.use_motz_wise_correction
        def __set__(self, mw):
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
            r.use_motz_wise_correction = mw

    property sticking_species:
        """
        The name of the sticking species. Needed only for reactions with
        multiple non-surface reactant species, where the sticking species is
        ambiguous.
        """
        def __get__(self):
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
            return pystr(r.sticking_species)
        def __set__(self, species):
            cdef CxxBlowersMaselInterfaceReaction* r = <CxxBlowersMaselInterfaceReaction*>self.reaction
            r.sticking_species = stringify(species)


cdef class CustomReaction(Reaction):
    """
    A reaction which follows mass-action kinetics with a custom reaction rate.

    An example for the definition of a `CustomReaction` object is given as::

        rxn = CustomReaction(
            equation="H2 + O <=> H + OH",
            rate=lambda T: 38.7 * T**2.7 * exp(-3150.15428/T),
            kinetics=gas)

    Warning: this class is an experimental part of the Cantera API and
        may be changed or removed without notice.
    """
    _reaction_type = "custom-rate-function"

    def __init__(self, equation=None, rate=None, Kinetics kinetics=None,
                 init=True, **kwargs):

        if init and equation and kinetics:

            spec = {"equation": equation, "type": self._reaction_type}

            self._reaction = CxxNewReaction(dict_to_anymap(spec),
                                            deref(kinetics.kinetics))
            self.reaction = self._reaction.get()
            self.rate = CustomRate(rate)

    property rate:
        """ Get/Set the `CustomRate` object for this reaction. """
        def __get__(self):
            return self._rate
        def __set__(self, CustomRate rate):
            self._rate = rate
            cdef CxxCustomFunc1Reaction* r = <CxxCustomFunc1Reaction*>self.reaction
            r.setRate(self._rate._base)
