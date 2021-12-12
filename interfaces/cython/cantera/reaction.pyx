# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.


# dictionary to store reaction classes
cdef dict _reaction_class_registry = {}


# dictionary to store reaction rate classes
cdef dict _reaction_rate_class_registry = {}


cdef class ReactionRate:
    """
    Base class for ReactionRate objects.

    ReactionRate objects are used to calculate reaction rates and are associated
    with a Reaction object.
    """
    _reaction_rate_type = ""

    def __repr__(self):
        return f"<{type(self).__name__} at {id(self):0x}>"

    def __call__(self, double temperature):
        """
        Evaluate rate expression based on temperature.
        """
        return self.rate.eval(temperature)

    property type:
        """ Get the C++ ReactionRate type """
        def __get__(self):
            return pystr(self.rate.type())

    @staticmethod
    cdef wrap(shared_ptr[CxxReactionRate] rate):
        """
        Wrap a C++ Reaction object with a Python object of the correct derived type.
        """
        # ensure all reaction types are registered
        if not _reaction_rate_class_registry:
            def register_subclasses(cls):
                for c in cls.__subclasses__():
                    rate_type = getattr(c, "_reaction_rate_type")
                    _reaction_rate_class_registry[rate_type] = c
                    register_subclasses(c)

            # update global reaction class registry
            register_subclasses(ReactionRate)

        # identify class
        rate_type = pystr(rate.get().type())
        cls = _reaction_rate_class_registry.get(rate_type, ReactionRate)

        # wrap C++ reaction rate
        cdef ReactionRate rr
        rr = cls(init=False)
        rr._rate = rate
        rr.set_cxx_object()
        return rr

    cdef set_cxx_object(self):
        self.rate = self._rate.get()

    @classmethod
    def from_dict(cls, data):
        """
        Create a `ReactionRate` object from a dictionary corresponding to its YAML
        representation.

        An example for the creation of a `ReactionRate` from a dictionary is::

            rate = ReactionRate.from_dict(
                {"rate-constant": {"A": 38.7, "b": 2.7, "Ea": 26191840.0}})

        :param data:
            A dictionary corresponding to the YAML representation.
        """
        if cls._reaction_rate_type != "":
            raise TypeError(
                f"Class method 'from_dict' was invoked from '{cls.__name__}' but "
                "should be called from base class 'ReactionRate'")

        cdef CxxAnyMap any_map = dict_to_anymap(data)
        cxx_rate = CxxNewReactionRate(any_map)
        return ReactionRate.wrap(cxx_rate)

    @classmethod
    def from_yaml(cls, text):
        """
        Create a `ReactionRate` object from its YAML string representation.

        An example for the creation of a `ReactionRate` from a YAML string is::

            rate = ReactionRate.from_yaml(
                "rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}")

        Units for ``A`` require a unit system with length in ``m`` and quantity in
        ``kmol`` (standard Cantera units).

        :param text:
            The YAML reaction rate string.
        """
        if cls._reaction_rate_type != "":
            raise TypeError(
                f"Class method 'from_yaml' was invoked from '{cls.__name__}' but "
                "should be called from base class 'ReactionRate'")

        cdef CxxAnyMap any_map
        any_map = AnyMapFromYamlString(stringify(text))
        cxx_rate = CxxNewReactionRate(any_map)
        return ReactionRate.wrap(cxx_rate)

    property input_data:
        """
        Get input data for this reaction rate with its current parameter values.
        """
        def __get__(self):
            return anymap_to_dict(self.rate.parameters())


cdef class ArrheniusRate(ReactionRate):
    r"""
    A reaction rate coefficient which depends on temperature only and follows
    the modified Arrhenius form:

    .. math::

        k_f = A T^b \exp{-\tfrac{E}{RT}}

    where ``A`` is the `pre_exponential_factor`, ``b`` is the `temperature_exponent`,
    and ``Ea`` is the `activation_energy`.
    """
    _reaction_rate_type = "Arrhenius"

    def __cinit__(self, A=None, b=None, Ea=None, input_data=None, init=True):

        if init:
            if isinstance(input_data, dict):
                self._rate.reset(new CxxArrheniusRate(dict_to_anymap(input_data)))
            elif all([arg is not None for arg in [A, b, Ea]]):
                self._rate.reset(new CxxArrheniusRate(A, b, Ea))
            elif all([arg is None for arg in [A, b, Ea, input_data]]):
                self._rate.reset(new CxxArrheniusRate(dict_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameters 'A', 'b' or 'Ea'")
            self.set_cxx_object()

    cdef CxxArrheniusRate* cxx_object(self):
        return <CxxArrheniusRate*>self.rate

    property pre_exponential_factor:
        """
        The pre-exponential factor ``A`` in units of m, kmol, and s raised to
        powers depending on the reaction order.
        """
        def __get__(self):
            return self.cxx_object().preExponentialFactor()

    property temperature_exponent:
        """
        The temperature exponent ``b``.
        """
        def __get__(self):
            return self.cxx_object().temperatureExponent()

    property activation_energy:
        """
        The activation energy ``E`` [J/kmol].
        """
        def __get__(self):
            return self.cxx_object().activationEnergy()

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            return self.cxx_object().allowNegativePreExponentialFactor()
        def __set__(self, cbool allow):
            self.cxx_object().setAllowNegativePreExponentialFactor(allow)


cdef class BlowersMaselRate(ReactionRate):
    r"""
    A reaction rate coefficient which depends on temperature and enthalpy change
    of the reaction follows the Blowers-Masel approximation and modified Arrhenius form
    described in `ArrheniusRate`.
    """
    _reaction_rate_type = "Blowers-Masel"

    def __cinit__(self, A=None, b=None, Ea0=None, w=None, input_data=None, init=True):

        if init:
            if isinstance(input_data, dict):
                self._rate.reset(new CxxBlowersMaselRate(dict_to_anymap(input_data)))
            elif all([arg is not None for arg in [A, b, Ea0, w]]):
                self._rate.reset(new CxxBlowersMaselRate(A, b, Ea0, w))
            elif all([arg is None for arg in [A, b, Ea0, w, input_data]]):
                self._rate.reset(new CxxBlowersMaselRate(dict_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameters 'A', 'b', 'Ea0' or 'w'")
            self.set_cxx_object()

    def __call__(self, double temperature, double deltaH):
        """
        Evaluate rate expression based on temperature and enthalpy change of reaction.
        """
        return self.rate.eval(temperature, deltaH)

    cdef CxxBlowersMaselRate* cxx_object(self):
        return <CxxBlowersMaselRate*>self.rate

    property pre_exponential_factor:
        """
        The pre-exponential factor ``A`` in units of m, kmol, and s raised to
        powers depending on the reaction order.
        """
        def __get__(self):
            return self.cxx_object().preExponentialFactor()

    property temperature_exponent:
        """
        The temperature exponent ``b``.
        """
        def __get__(self):
            return self.cxx_object().temperatureExponent()

    def activation_energy(self, double deltaH):
        """
        The activation energy ``E`` [J/kmol], based on enthalpy change of reaction
        ``deltaH`` (in [J/kmol])
        """
        return self.cxx_object().activationEnergy(deltaH)

    property intrinsic_activation_energy:
        """
        The intrinsic activation energy ``E0`` [J/kmol].
        """
        def __get__(self):
            return self.cxx_object().activationEnergy0()

    property bond_energy:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        def __get__(self):
            return self.cxx_object().bondEnergy()

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            return self.cxx_object().allowNegativePreExponentialFactor()
        def __set__(self, cbool allow):
            self.cxx_object().setAllowNegativePreExponentialFactor(allow)


cdef class FalloffRate(ReactionRate):
    """
    Base class for parameterizations used to describe the fall-off in reaction rates
    due to intermolecular energy transfer. These objects are used by `FalloffReaction`.
    Note that `FalloffRate` is a base class for specialized fall-off parameterizations
    and cannot be instantiated by itself.
    """

    _reaction_rate_type = None

    def __cinit__(self, low=None, high=None, falloff_coeffs=None, input_data=None, init=True):

        if self._reaction_rate_type is None:
            raise TypeError(
                f"Base class '{self.__class__.__name__}' cannot be instantiated "
                "by itself; use specialized fall-off rate instead.")

        if init:
            if isinstance(input_data, dict):
                if input_data:
                    self._from_dict(input_data)
                else:
                    self._from_empty()
            elif input_data is None:
                self._from_empty()
            else:
                raise TypeError("Invalid input parameters")
            self.set_cxx_object()

            if (low is not None) and (high is not None):
                self.low_rate = low
                self.high_rate = high
                if falloff_coeffs is None:
                    falloff_coeffs = ()
                self.falloff_coeffs = falloff_coeffs

    def __call__(self, double temperature, double concm):
        """
        Evaluate rate expression based on temperature and third-body concentration.
        """
        return self.rate.eval(temperature, concm)

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = <CxxFalloffRate*>self.rate

    property low_rate:
        """ Get/Set the `Arrhenius` rate constant in the low-pressure limit """
        def __get__(self):
            return Arrhenius.wrap(&(self.falloff.lowRate()))
        def __set__(self, Arrhenius rate):
            self.falloff.setLowRate(deref(rate.rate))

    property high_rate:
        """ Get/Set the `Arrhenius` rate constant in the high-pressure limit """
        def __get__(self):
            return Arrhenius.wrap(&(self.falloff.highRate()))
        def __set__(self, Arrhenius rate):
            self.falloff.setHighRate(deref(rate.rate))

    property falloff_coeffs:
        """ The array of coefficients used to define this falloff function. """
        def __get__(self):
            cdef vector[double] cxxdata
            self.falloff.getFalloffCoeffs(cxxdata)
            return np.fromiter(cxxdata, np.double)
        def __set__(self, data):
            cdef vector[double] cxxdata
            for c in data:
                cxxdata.push_back(c)
            self.falloff.setFalloffCoeffs(cxxdata)

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            return self.falloff.allowNegativePreExponentialFactor()
        def __set__(self, cbool allow):
            self.falloff.setAllowNegativePreExponentialFactor(allow)

    property chemically_activated:
        """
        Get whether the object is a chemically-activated reaction rate.
        """
        def __get__(self):
            return self.falloff.chemicallyActivated()
        def __set__(self, cbool activated):
            self.falloff.setChemicallyActivated(activated)

    def falloff_function(self, double temperature, conc3b=None):
        """
        Evaluate the falloff function
        """
        if conc3b is None:
            conc3b = self.third_body_concentration
        return self.falloff.evalF(temperature, conc3b)


cdef class LindemannRate(FalloffRate):
    r"""
    The Lindemann falloff parameterization.

    This class implements the simple falloff function :math:`F(T,P_r) = 1.0`.
    """
    _reaction_rate_type = "Lindemann"

    def _from_empty(self):
        self._rate.reset(new CxxLindemannRate())

    def _from_dict(self, dict input_data):
        self._rate = CxxNewReactionRate(dict_to_anymap(input_data))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = <CxxFalloffRate*>self.rate


cdef class TroeRate(FalloffRate):
    r"""
    The 3- or 4-parameter Troe falloff function.

    :param falloff_coeffs:
        An array of 3 or 4 parameters: :math:`[a, T^{***}, T^*, T^{**}]` where
        the final parameter is optional (with a default value of 0).
    """
    _reaction_rate_type = "Troe"

    def _from_empty(self):
        self._rate.reset(new CxxTroeRate())

    def _from_dict(self, dict input_data):
        self._rate = CxxNewReactionRate(dict_to_anymap(input_data))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = <CxxFalloffRate*>self.rate


cdef class SriRate(FalloffRate):
    r"""
    The 3- or 5-parameter SRI falloff function.

    :param falloff_coeffs:
        An array of 3 or 5 parameters: :math:`[a, b, c, d, e]` where the last
        two parameters are optional (with default values of 1 and 0, respectively).
    """
    _reaction_rate_type = "SRI"

    def _from_empty(self):
        self._rate.reset(new CxxSriRate())

    def _from_dict(self, dict input_data):
        self._rate = CxxNewReactionRate(dict_to_anymap(input_data))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = <CxxFalloffRate*>self.rate


cdef class TsangRate(FalloffRate):
    r"""
    The Tsang falloff parameterization.
    """
    _reaction_rate_type = "Tsang"

    def _from_empty(self):
        self._rate.reset(new CxxTsangRate())

    def _from_dict(self, dict input_data):
        self._rate = CxxNewReactionRate(dict_to_anymap(input_data))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = <CxxFalloffRate*>self.rate


cdef class PlogRate(ReactionRate):
    r"""
    A pressure-dependent reaction rate parameterized by logarithmically
    interpolating between Arrhenius rate expressions at various pressures.
    """
    _reaction_rate_type = "pressure-dependent-Arrhenius"

    def __cinit__(self, rates=None, input_data=None, init=True):

        if init and isinstance(rates, list):
            self.rates = rates

        elif init:
            if isinstance(input_data, dict):
                self._rate.reset(new CxxPlogRate(dict_to_anymap(input_data)))
            elif rates is None:
                self._rate.reset(new CxxPlogRate(dict_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameter 'rates'")
            self.set_cxx_object()

    def __call__(self, double temperature, double pressure):
        """
        Evaluate rate expression based on temperature and pressure.
        """
        return self.rate.eval(temperature, pressure)

    cdef CxxPlogRate* cxx_object(self):
        return <CxxPlogRate*>self.rate

    property rates:
        """
        Get/Set the rate coefficients for this reaction, which are given as a
        list of (pressure, `Arrhenius`) tuples.
        """
        def __get__(self):
            rates = []
            cdef multimap[double, CxxArrhenius2] cxxrates
            cdef pair[double, CxxArrhenius2] p_rate
            cxxrates = self.cxx_object().getRates()
            for p_rate in cxxrates:
                rates.append((p_rate.first, copyArrhenius(&p_rate.second)))
            return rates

        def __set__(self, rates):
            cdef multimap[double, CxxArrhenius2] ratemap
            cdef Arrhenius rate
            cdef pair[double, CxxArrhenius2] item
            for p, rate in rates:
                item.first = p
                if rate.legacy is not NULL:
                    item.second = deref(rate.legacy)
                else:
                    item.second = CxxArrhenius2(
                        rate.rate.preExponentialFactor(),
                        rate.rate.temperatureExponent(),
                        rate.rate.activationEnergy() / gas_constant
                    )
                ratemap.insert(item)

            self._rate.reset(new CxxPlogRate(ratemap))
            self.rate = self._rate.get()


cdef class ChebyshevRate(ReactionRate):
    r"""
    A pressure-dependent reaction rate parameterized by a bivariate Chebyshev
    polynomial in temperature and pressure.
    """
    _reaction_rate_type = "Chebyshev"

    def __cinit__(self, temperature_range=None, pressure_range=None, data=None,
                  input_data=None, init=True):

        if init:
            if isinstance(input_data, dict):
                self._rate.reset(new CxxChebyshevRate3(dict_to_anymap(input_data)))
            elif all([arg is not None
                    for arg in [temperature_range, pressure_range, data]]):
                Tmin = temperature_range[0]
                Tmax = temperature_range[1]
                Pmin = pressure_range[0]
                Pmax = pressure_range[1]
                self._rate.reset(
                    new CxxChebyshevRate3(Tmin, Tmax, Pmin, Pmax, self._cxxarray2d(data)))
            elif all([arg is None
                    for arg in [temperature_range, pressure_range, data, input_data]]):
                self._rate.reset(new CxxChebyshevRate3(dict_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameters")
            self.set_cxx_object()

    def __call__(self, double temperature, double pressure):
        """
        Evaluate rate expression based on temperature and pressure.
        """
        return self.rate.eval(temperature, pressure)

    cdef CxxArray2D _cxxarray2d(self, coeffs):
        """ Internal function to assign coefficient matrix values """
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

    cdef CxxChebyshevRate3* cxx_object(self):
        return <CxxChebyshevRate3*>self.rate

    property temperature_range:
        """ Valid temperature range [K] for the Chebyshev fit """
        def __get__(self):
            return self.cxx_object().Tmin(), self.cxx_object().Tmax()

    property pressure_range:
        """ Valid pressure range [Pa] for the Chebyshev fit """
        def __get__(self):
            return self.cxx_object().Pmin(), self.cxx_object().Pmax()

    property n_temperature:
        """
        Number of temperatures over which the Chebyshev fit is computed.
        (same as number of rows of `data` property).
        """
        def __get__(self):
            return self.cxx_object().Pmin(), self.cxx_object().nTemperature()

    property n_pressure:
        """
        Number of pressures over which the Chebyshev fit is computed
        (same as number of columns of `data` property).
        """
        def __get__(self):
            return self.cxx_object().Pmin(), self.cxx_object().nPressure()

    property data:
        """
        2D array of Chebyshev coefficients where rows and columns correspond to
        temperature and pressure dimensions over which the Chebyshev fit is computed.
        """
        def __get__(self):
            cdef CxxArray2D cxxcoeffs = self.cxx_object().data()
            c = np.fromiter(cxxcoeffs.data(), np.double)
            return c.reshape(cxxcoeffs.nRows(), cxxcoeffs.nColumns(), order="F")


cdef class CustomRate(ReactionRate):
    r"""
    A custom rate coefficient which depends on temperature only.

    The simplest way to create a `CustomRate` object is to use a lambda function,
    for example::

        rr = CustomRate(lambda T: 38.7 * T**2.7 * exp(-3150.15/T))

    Warning: this class is an experimental part of the Cantera API and
        may be changed or removed without notice.
    """
    _reaction_rate_type = "custom-rate-function"

    def __cinit__(self, k=None, init=True):

        if init:
            self._rate.reset(new CxxCustomFunc1Rate())
            self.set_cxx_object()
            try:
                self.set_rate_function(k)
            except Exception:
                raise TypeError(f"Cannot convert input with type '{type(k)}' to rate expression.")

    cdef CxxCustomFunc1Rate* cxx_object(self):
        return <CxxCustomFunc1Rate*>self.rate

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

        self.cxx_object().setRateFunction(self._rate_func._func)


cdef class Reaction:
    """
    A class which stores data about a reaction and its rate parameterization so
    that it can be added to a `Kinetics` object.

    :param reactants:
        Value used to set `reactants`
    :param products:
        Value used to set `products`

    The static methods `listFromFile`, `list_from_yaml`, `listFromCti`, and
    `listFromXml` can be used to create lists of `Reaction` objects from
    existing definitions in the YAML, CTI, or XML formats. All of the following
    will produce a list of the 325 reactions which make up the GRI 3.0
    mechanism::

        R = ct.Reaction.listFromFile("gri30.yaml", gas)
        R = ct.Reaction.listFromCti(open("path/to/gri30.cti").read())
        R = ct.Reaction.listFromXml(open("path/to/gri30.xml").read())

    where `gas` is a `Solution` object with the appropriate thermodynamic model,
    which is the `ideal-gas` model in this case.

    The static method `list_from_yaml` can be used to create lists of `Reaction`
    objects from a YAML list::

        rxns = '''
          - equation: O + H2 <=> H + OH
            rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
          - equation: O + HO2 <=> OH + O2
            rate-constant: {A: 2.0e+13, b: 0.0, Ea: 0.0}
        '''
        R = ct.Reaction.list_from_yaml(rxns, gas)

    The methods `from_yaml`, `fromCti`, and `fromXml` can be used to create
    individual `Reaction` objects from definitions in these formats. In the case
    of using YAML or CTI definitions, it is important to verify that either the
    pre-exponential factor and activation energy are supplied in SI units, or
    that they have their units specified::

        R = ct.Reaction.from_yaml('''{equation: O + H2 <=> H + OH,
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
        if not _reaction_class_registry:
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
    def from_dict(cls, data, Kinetics kinetics):
        """
        Create a `Reaction` object from a dictionary corresponding to its YAML
        representation.

        An example for the creation of a Reaction from a dictionary is::

            rxn = Reaction.from_dict(
                {"equation": "O + H2 <=> H + OH",
                 "rate-constant": {"A": 38.7, "b": 2.7, "Ea": 26191840.0}},
                kinetics=gas)

        In the example, ``gas`` is a Kinetics (or Solution) object.

        :param data:
            A dictionary corresponding to the YAML representation.
        :param kinetics:
            A `Kinetics` object whose associated phase(s) contain the species
            involved in the reaction.
        """
        if cls._reaction_type != "":
            raise TypeError(
                f"Class method 'from_dict' was invoked from '{cls.__name__}' but "
                "should be called from base class 'Reaction'")

        cdef CxxAnyMap any_map = dict_to_anymap(data)
        cxx_reaction = CxxNewReaction(any_map, deref(kinetics.kinetics))
        return Reaction.wrap(cxx_reaction)

    @classmethod
    def fromYaml(cls, text, Kinetics kinetics=None):
        """
        Create a `Reaction` object from its YAML string representation.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by `Reaction.from_yaml`.
        """
        warnings.warn("Class method 'fromYaml' is renamed to 'from_yaml' "
            "and will be removed after Cantera 2.6.", DeprecationWarning)

        return cls.from_yaml(text, kinetics)

    @classmethod
    def from_yaml(cls, text, Kinetics kinetics):
        """
        Create a `Reaction` object from its YAML string representation.

        An example for the creation of a Reaction from a YAML string is::

            rxn = Reaction.from_yaml('''
                equation: O + H2 <=> H + OH
                rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0 cal/mol}
                ''', kinetics=gas)

        In the example, ``gas`` is a Kinetics (or Solution) object.

        :param text:
            The YAML reaction string.
        :param kinetics:
            A `Kinetics` object whose associated phase(s) contain the species
            involved in the reaction.
        """
        if cls._reaction_type != "":
            raise TypeError(
                f"Class method 'from_yaml' was invoked from '{cls.__name__}' but "
                "should be called from base class 'Reaction'")

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
        argument, and reactions from the section ``section`` will be returned.

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

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by `Reaction.list_from_yaml`.
        """
        warnings.warn("Class method 'listFromYaml' is renamed to 'list_from_yaml' "
            "and will be removed after Cantera 2.6.", DeprecationWarning)

        return Reaction.list_from_yaml(text, kinetics)

    @staticmethod
    def list_from_yaml(text, Kinetics kinetics):
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

    def _deprecation_warning(self, attr, what="property", new=None):
        if new:
            return (f"\n{what.capitalize()} '{attr}' to be removed after Cantera 2.6.\n"
                    f"This {what} is moved to the {type(self.rate).__name__} object "
                    f"accessed via the 'rate' property as '{new}'.")
        return (f"\n{what.capitalize()} '{attr}' to be removed after Cantera 2.6.\n"
                f"This {what} is moved to the {type(self.rate).__name__} object accessed "
                "via the 'rate' property.")

    property uses_legacy:
        """Indicate whether reaction uses a legacy implementation"""
        def __get__(self):
            return self.reaction.usesLegacy()

    property rate_coeff_units:
        """Get reaction rate coefficient units"""
        def __get__(self):
            cdef CxxUnits rate_units = self.reaction.rate_units
            return Units.copy(rate_units)


cdef class Arrhenius:
    r"""
    A reaction rate coefficient which depends on temperature only and follows
    the modified Arrhenius form:

    .. math::

        k_f = A T^b \exp{-\tfrac{E}{RT}}

    where ``A`` is the `pre_exponential_factor`, ``b`` is the `temperature_exponent`,
    and ``E`` is the `activation_energy`.
    """
    def __cinit__(self, A=0, b=0, E=0, init=True):
        if init:
            self.rate = new CxxArrheniusRate(A, b, E)
            self.base = self.rate
            self.own_rate = True
            self.reaction = None
        else:
            self.own_rate = False

    def __dealloc__(self):
        if self.own_rate:
            del self.base

    @staticmethod
    cdef wrap(CxxArrheniusRate* rate):
        r = Arrhenius(init=False)
        r.rate = rate
        r.base = rate
        r.reaction = None
        return r

    property pre_exponential_factor:
        """
        The pre-exponential factor ``A`` in units of m, kmol, and s raised to
        powers depending on the reaction order.
        """
        def __get__(self):
            return self.base.preExponentialFactor()

    property temperature_exponent:
        """
        The temperature exponent ``b``.
        """
        def __get__(self):
            return self.base.temperatureExponent()

    property activation_energy:
        """
        The activation energy ``E`` [J/kmol].
        """
        def __get__(self):
            return self.base.intrinsicActivationEnergy()

    def __repr__(self):
        return 'Arrhenius(A={:g}, b={:g}, E={:g})'.format(
            self.pre_exponential_factor, self.temperature_exponent,
            self.activation_energy)

    def __call__(self, float T):
        if self.rate != NULL:
            return self.base.eval(T)
        else:
            return self.legacy.updateRC(np.log(T), 1/T)


cdef wrapArrhenius(CxxArrheniusBase* rate, Reaction reaction):
    r = Arrhenius(init=False)
    r.base = rate
    if reaction.uses_legacy:
        r.legacy = <CxxArrhenius2*>r.base
        r.rate = NULL
    else:
        r.rate = <CxxArrheniusRate*>r.base
        r.legacy = NULL

    r.reaction = reaction
    return r

cdef copyArrhenius(CxxArrhenius2* rate):
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

    cdef CxxElementaryReaction3* cxx_object(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxElementaryReaction3*>self.reaction

    cdef CxxElementaryReaction2* cxx_object2(self):
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
        cdef CxxElementaryReaction2* r = self.cxx_object2()
        r.rate = deref(<CxxArrhenius2*>rate.base)

    property rate:
        """ Get/Set the `ArrheniusRate` rate coefficient for this reaction. """
        def __get__(self):
            if self.uses_legacy:
                return wrapArrhenius(&(self.cxx_object2().rate), self)

            return ArrheniusRate.wrap(self.cxx_object().rate())
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
            self.cxx_object().setRate(rate3._rate)

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
                return self.cxx_object2().allow_negative_pre_exponential_factor

            attr = "allow_negative_pre_exponential_factor"
            warnings.warn(self._deprecation_warning(attr), DeprecationWarning)
            return self.rate.allow_negative_pre_exponential_factor
        def __set__(self, allow):
            if self.uses_legacy:
                self.cxx_object2().allow_negative_pre_exponential_factor = allow
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

    cdef CxxThreeBodyReaction3* cxx_threebody(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxThreeBodyReaction3*>self.reaction

    cdef CxxThreeBodyReaction2* cxx_threebody2(self):
        if not self.uses_legacy:
            raise AttributeError("Incorrect accessor for legacy implementation")
        return <CxxThreeBodyReaction2*>self.reaction

    cdef CxxThirdBody* thirdbody(self):
        if self.uses_legacy:
            return &(self.cxx_threebody2().third_body)
        return <CxxThirdBody*>(self.cxx_threebody().thirdBody().get())

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
        Get the efficiency of the third body named ``species`` considering both
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
        Used internally when wrapping :ct:`FalloffRate` objects returned from C++.
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

    An example for the definition of a `FalloffReaction` object is given as::

        rxn = FalloffReaction(
            equation="2 OH (+ M) <=> H2O2 (+ M)",
            type="falloff",
            rate=ct.TroeRate(low=ct.Arrhenius(2.3e+12, -0.9, -7112800.0),
                             high=ct.Arrhenius(7.4e+10, -0.37, 0),
                             data=[0.7346, 94.0, 1756.0, 5182.0])
            efficiencies={"AR": 0.7, "H2": 2.0, "H2O": 6.0},
            kinetics=gas)

    The YAML description corresponding to this reaction is::

        equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 3
        type: falloff-legacy
        low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0 cal/mol}
        high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0 cal/mol}
        Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
        efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
    """
    _reaction_type = "falloff"
    _has_legacy = True
    _hybrid = True

    cdef CxxFalloffReaction3* cxx_object(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxFalloffReaction3*>self.reaction

    cdef CxxFalloffReaction2* cxx_object2(self):
        if not self.uses_legacy:
            raise AttributeError("Incorrect accessor for legacy implementation")
        return <CxxFalloffReaction2*>self.reaction

    cdef CxxThirdBody* thirdbody(self):
        if self.uses_legacy:
            return &(self.cxx_object2().third_body)
        return <CxxThirdBody*>(self.cxx_object().thirdBody().get())

    def __init__(
            self, equation=None, rate=None, efficiencies=None,
            Kinetics kinetics=None, init=True, legacy=False, **kwargs
        ):

        if init and equation and kinetics:

            rxn_type = self._reaction_type
            if legacy:
                rxn_type += "-legacy"
            spec = {"equation": equation, "type": rxn_type}
            if isinstance(rate, dict):
                for key, value in rate.items():
                    spec[key.replace("_", "-")] = value
            elif legacy and isinstance(rate, Falloff):
                raise NotImplementedError("Not implemented for legacy rates")
            elif not legacy and (isinstance(rate, FalloffRate) or rate is None):
                pass
            else:
                raise TypeError(
                    f"Invalid rate definition; type is '{type(rate).__name__}'")

            if isinstance(efficiencies, dict):
                spec["efficiencies"] = efficiencies

            self._reaction = CxxNewReaction(dict_to_anymap(spec),
                                            deref(kinetics.kinetics))
            self.reaction = self._reaction.get()

            if not legacy and isinstance(rate, FalloffRate):
                self.rate = rate

    property rate:
        """ Get/Set the `FalloffRate` rate coefficients for this reaction. """
        def __get__(self):
            if not self.uses_legacy:
                return FalloffRate.wrap(self.cxx_object().rate())
            raise AttributeError("Legacy implementation does not use rate property.")
        def __set__(self, FalloffRate rate):
            if not self.uses_legacy:
                self.cxx_object().setRate(rate._rate)
                return
            raise AttributeError("Legacy implementation does not use rate property.")

    property low_rate:
        """ Get/Set the `Arrhenius` rate constant in the low-pressure limit """
        def __get__(self):
            if self.uses_legacy:
                return wrapArrhenius(&(self.cxx_object2().low_rate), self)
            warnings.warn(self._deprecation_warning("low_rate"), DeprecationWarning)
            return self.rate.low_rate
        def __set__(self, Arrhenius rate):
            if self.uses_legacy:
                self.cxx_object2().low_rate = deref(rate.legacy)
            warnings.warn(self._deprecation_warning("low_rate"), DeprecationWarning)
            self.rate.low_rate = rate

    property high_rate:
        """ Get/Set the `Arrhenius` rate constant in the high-pressure limit """
        def __get__(self):
            if self.uses_legacy:
                return wrapArrhenius(&(self.cxx_object2().high_rate), self)
            warnings.warn(self._deprecation_warning("high_rate"), DeprecationWarning)
            return self.rate.high_rate
        def __set__(self, Arrhenius rate):
            if self.uses_legacy:
                self.cxx_object2().high_rate = deref(rate.legacy)
            warnings.warn(self._deprecation_warning("high_rate"), DeprecationWarning)
            self.rate.high_rate = rate

    property falloff:
        """
        Get/Set the `Falloff` function used to blend the high- and low-pressure
        rate coefficients
        """
        def __get__(self):
            if self.uses_legacy:
                return wrapFalloff(self.cxx_object2().falloff)
            raise AttributeError("New implementation uses 'FalloffRate' objects that "
                "are retrieved using the 'rate' property.")
        def __set__(self, Falloff f):
            if self.uses_legacy:
                self.cxx_object2().falloff = f._falloff
            raise TypeError("New implementation requires 'FalloffRate' objects that "
                "are set using the 'rate' property.")

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

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `FalloffRate.allow_negative_pre_exponential_factor`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cxx_object2().allow_negative_pre_exponential_factor

            attr = "allow_negative_pre_exponential_factor"
            warnings.warn(self._deprecation_warning(attr), DeprecationWarning)
            return self.rate.allow_negative_pre_exponential_factor
        def __set__(self, allow):
            if self.uses_legacy:
                self.cxx_object2().allow_negative_pre_exponential_factor = allow
                return

            attr = "allow_negative_pre_exponential_factor"
            warnings.warn(self._deprecation_warning(attr), DeprecationWarning)
            self.rate.allow_negative_pre_exponential_factor = allow

    def efficiency(self, species):
        """
        Get the efficiency of the third body named ``species`` considering both
        the default efficiency and species-specific efficiencies.
        """
        return self.thirdbody().efficiency(stringify(species))


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

    cdef CxxPlogReaction3* cxx_object(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxPlogReaction3*>self.reaction

    cdef CxxPlogReaction2* cxx_object2(self):
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
            return PlogRate.wrap(self.cxx_object().rate())
        def __set__(self, PlogRate rate):
            if self.uses_legacy:
                raise AttributeError("Legacy implementation does not use rate property.")
            self.cxx_object().setRate(rate._rate)

    cdef list _legacy_get_rates(self):
        cdef CxxPlogReaction2* r = self.cxx_object2()
        cdef multimap[double,CxxArrhenius2] cxxrates = r.rate.getRates()
        cdef pair[double,CxxArrhenius2] p_rate
        rates = []
        for p_rate in cxxrates:
            rates.append((p_rate.first,copyArrhenius(&p_rate.second)))
        return rates

    cdef _legacy_set_rates(self, list rates):
        cdef multimap[double,CxxArrhenius2] ratemap
        cdef Arrhenius rate
        cdef pair[double,CxxArrhenius2] item
        for p, rate in rates:
            item.first = p
            if rate.legacy is not NULL:
                item.second = deref(rate.legacy)
            else:
                item.second = CxxArrhenius2(
                    rate.rate.preExponentialFactor(),
                    rate.rate.temperatureExponent(),
                    rate.rate.activationEnergy() / gas_constant
                )
            ratemap.insert(item)

        cdef CxxPlogReaction2* r = self.cxx_object2()
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
        cdef CxxPlogReaction2* r = self.cxx_object2()
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
            rate={"temperature-range": [290.0, 3000.0],
                  "pressure-range": [1e3, 1e8],
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

    cdef CxxChebyshevReaction3* cxx_object(self):
        if self.uses_legacy:
            raise AttributeError("Incorrect accessor for updated implementation")
        return <CxxChebyshevReaction3*>self.reaction

    cdef CxxChebyshevReaction2* cxx_object2(self):
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
                for key, value in rate.items():
                    spec[key.replace("_", "-")] = value
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
            return ChebyshevRate.wrap(self.cxx_object().rate())
        def __set__(self, ChebyshevRate rate):
            if self.uses_legacy:
                raise AttributeError("Legacy implementation does not use rate property.")
            self.cxx_object().setRate(rate._rate)

    property Tmin:
        """
        Minimum temperature [K] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.temperature_range[0]`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cxx_object2().rate.Tmin()

            warnings.warn(
                self._deprecation_warning(
                    "Tmin", new="ChebyshevRate.temperature_range[0]"),
                DeprecationWarning)
            return self.rate.temperature_range[0]

    property Tmax:
        """
        Maximum temperature [K] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.temperature_range[1]`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cxx_object2().rate.Tmax()

            warnings.warn(
                self._deprecation_warning(
                    "Tmax", new="ChebyshevRate.temperature_range[1]"),
                DeprecationWarning)
            return self.rate.temperature_range[1]

    property Pmin:
        """
        Minimum pressure [Pa] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.pressure_range[0]`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cxx_object2().rate.Pmin()

            warnings.warn(
                self._deprecation_warning(
                    "Pmin", new="ChebyshevRate.pressure_range[0]"),
                DeprecationWarning)
            return self.rate.pressure_range[0]

    property Pmax:
        """ Maximum pressure [Pa] for the Chebyshev fit

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.pressure_range[1]`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cxx_object2().rate.Pmax()

            warnings.warn(
                self._deprecation_warning(
                    "Pmax", new="ChebyshevRate.pressure_range[0]"),
                DeprecationWarning)
            return self.rate.pressure_range[1]

    property nPressure:
        """
        Number of pressures over which the Chebyshev fit is computed

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.n_pressure`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cxx_object2().rate.nPressure()

            warnings.warn(
                self._deprecation_warning(
                    "nPressure", new="ChebyshevRate.data.n_pressure"),
                DeprecationWarning)
            return self.rate.data.shape[1]

    property nTemperature:
        """
        Number of temperatures over which the Chebyshev fit is computed

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.n_temperature`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self.cxx_object2().rate.nTemperature()

            warnings.warn(
                self._deprecation_warning(
                    "nTemperature", new="ChebyshevRate.data.n_temperature"),
                DeprecationWarning)
            return self.rate.data.shape[0]

    cdef _legacy_get_coeffs(self):
        cdef CxxChebyshevReaction2* r = self.cxx_object2()
        cdef CxxArray2D cxxcoeffs = r.rate.data()
        c = np.fromiter(cxxcoeffs.data(), np.double)
        return c.reshape(cxxcoeffs.nRows(), cxxcoeffs.nColumns(), order="F")

    property coeffs:
        """
        2D array of Chebyshev coefficients of size `(n_temperature, n_pressure)`.

        .. deprecated:: 2.6
             To be deprecated with version 2.6, and removed thereafter.
             Replaced by property `ChebyshevRate.data`.
        """
        def __get__(self):
            if self.uses_legacy:
                return self._legacy_get_coeffs()

            warnings.warn(
                self._deprecation_warning(
                    "coeffs", new="ChebyshevRate.data"),
                DeprecationWarning)
            return self.rate.data

    cdef _legacy_set_parameters(self, Tmin, Tmax, Pmin, Pmax, coeffs):
        cdef CxxChebyshevReaction2* r = self.cxx_object2()

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
        cdef pair[double,double] Trange
        Trange.first, Trange.second = Tmin, Tmax
        cdef pair[double,double] Prange
        Prange.first, Prange.second = Pmin, Pmax
        self.rate = ChebyshevRate(Trange, Prange, coeffs)

    cdef _legacy_call(self, float T, float P):
        cdef CxxChebyshevReaction2* r = self.cxx_object2()
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
            self.rate = new CxxBlowersMasel2(A, b, E0 / gas_constant, w / gas_constant)
            self.own_rate = True
            self.reaction = None
        else:
            self.own_rate = False

    def __dealloc__(self):
        if self.own_rate:
            del self.rate

    property pre_exponential_factor:
        """
        The pre-exponential factor ``A`` in units of m, kmol, and s raised to
        powers depending on the reaction order.
        """
        def __get__(self):
            return self.rate.preExponentialFactor()

    property temperature_exponent:
        """
        The temperature exponent ``b``.
        """
        def __get__(self):
            return self.rate.temperatureExponent()

    def activation_energy(self, float dH):
        """
        The activation energy ``E`` [J/kmol]

        :param dH: The enthalpy of reaction [J/kmol] at the current temperature
        """
        return self.rate.activationEnergy_R(dH) * gas_constant

    property bond_energy:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        def __get__(self):
            return self.rate.bondEnergy() * gas_constant

    property intrinsic_activation_energy:
        """
        The intrinsic activation energy ``E0`` [J/kmol].
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


cdef wrapBlowersMasel(CxxBlowersMasel2* rate, Reaction reaction):
    r = BlowersMasel(init=False)
    r.rate = rate
    r.reaction = reaction
    return r


cdef class BlowersMaselReaction(Reaction):
    """
    A reaction which follows mass-action kinetics with Blowers Masel
    reaction rate.

    An example for the definition of an `BlowersMaselReaction` object is given as::
        rxn = BlowersMaselReaction(
            equation="O + H2 <=> H + OH",
            rate={"A": 38.7, "b": 2.7, "Ea0": 1.0958665856e8, "w": 1.7505856e13},
            kinetics=gas)
    The YAML description corresponding to this reaction is::
        equation: O + H2 <=> H + OH
        type: Blowers-Masel
        rate-constant: {A: 38700 cm^3/mol/s, b: 2.7, Ea0: 2.619184e4 cal/mol, w: 4.184e9 cal/mol}
    """
    _reaction_type = "Blowers-Masel"

    def __init__(self, equation=None, rate=None, Kinetics kinetics=None,
                 init=True, **kwargs):

        if init and equation and kinetics:

            rxn_type = self._reaction_type
            spec = {"equation": equation, "type": rxn_type}
            if isinstance(rate, dict):
                spec["rate-constant"] = rate
            elif rate is None or isinstance(rate, BlowersMaselRate):
                pass
            else:
                raise TypeError("Invalid rate definition")

            self._reaction = CxxNewReaction(dict_to_anymap(spec),
                                            deref(kinetics.kinetics))
            self.reaction = self._reaction.get()

            if isinstance(rate, BlowersMaselRate):
                self.rate = rate

    property rate:
        """ Get/Set the `BlowersMaselRate` rate coefficients for this reaction. """
        def __get__(self):
            cdef CxxBlowersMaselReaction* r = <CxxBlowersMaselReaction*>self.reaction
            return BlowersMaselRate.wrap(r.rate())
        def __set__(self, rate):
            cdef CxxBlowersMaselReaction* r = <CxxBlowersMaselReaction*>self.reaction
            cdef BlowersMaselRate rate_ = rate
            r.setRate(rate_._rate)


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


cdef class BlowersMaselInterfaceReaction(Reaction):
    """
    A reaction occurring on an `Interface` (i.e. a surface or an edge)
    with the rate parameterization of `BlowersMasel`.
    """
    _reaction_type = "surface-Blowers-Masel"

    cdef CxxBlowersMaselInterfaceReaction* cxx_object2(self):
        return <CxxBlowersMaselInterfaceReaction*>self.reaction

    property rate:
        """ Get/Set the `BlowersMaselRate` rate coefficients for this reaction. """
        def __get__(self):
            return wrapBlowersMasel(&(self.cxx_object2().rate), self)
        def __set__(self, BlowersMasel rate):
            cdef CxxBlowersMaselInterfaceReaction* r = self.cxx_object2()
            r.rate = deref(rate.rate)

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            return self.cxx_object2().allow_negative_pre_exponential_factor
        def __set__(self, allow):
            self.cxx_object2().allow_negative_pre_exponential_factor = allow

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
            r.setRate(self._rate._rate)
