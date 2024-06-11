# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

cimport numpy as np
import numpy as np
import warnings
from cython.operator cimport dereference as deref

from .kinetics cimport Kinetics
from ._utils cimport *
from ._utils import CanteraError
from .units cimport *
from .delegator cimport *

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

    property sub_type:
        """ Get the C++ ReactionRate sub-type """
        def __get__(self):
            return pystr(self.rate.subType())

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

        # Check for delegated rate, which will already have a Python wrapper
        cdef CxxDelegator* drate = dynamic_cast[CxxDelegatorPtr](rate.get())
        cdef CxxPythonHandle* handle

        if drate != NULL:
            handle = dynamic_pointer_cast[CxxPythonHandle, CxxExternalHandle](
                drate.getExternalHandle(stringify("python"))).get()
            if handle != NULL and handle.get() != NULL:
                py_rate = <ReactionRate>(<PyObject*>handle.get())
                py_rate._rate = rate
                return py_rate
            else:
                raise CanteraError("Internal Error: Delegator does not have a "
                                   "corresponding Python ExtensibleRate object")

        # identify class (either subType or type)
        rate_type = pystr(rate.get().subType())
        if rate_type == "":
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
    def from_dict(cls, data, hyphenize=True):
        """
        Create a `ReactionRate` object from a dictionary corresponding to its YAML
        representation. By default, underscores in keys are replaced by hyphens.

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

        cdef CxxAnyMap any_map = py_to_anymap(data, hyphenize=hyphenize)
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
            return anymap_to_py(self.rate.parameters())

    @property
    def conversion_units(self) -> Units:
        """
        Get the units for converting the leading term in the reaction rate expression
        to different unit systems.
        """
        return Units.copy(self.rate.conversionUnits())


cdef class ArrheniusRateBase(ReactionRate):
    """
    Base class collecting commonly used features of Arrhenius-type rate objects.
    Objects should be instantiated by specialized classes, for example `ArrheniusRate`,
    `BlowersMaselRate` and `TwoTempPlasmaRate`.
    """
    _reaction_rate_type = None

    def _cinit(self, input_data, **kwargs):
        """
        Helper function called by __cinit__. The method is used as a uniform interface
        for object construction. A separate method is necessary as default argument
        values to __cinit__ defined in derived classes are not available in the base
        class's __cinit__ (which gets called first).
        """
        if self._reaction_rate_type is None:
            raise TypeError(
                f"Base class '{self.__class__.__name__}' cannot be instantiated "
                "by itself; use specialized rate constructors instead.")

        if isinstance(input_data, dict):
            self._from_dict(input_data)
        elif all([kwargs[k] is not None for k in kwargs]):
            self._from_parameters(*[kwargs[k] for k in kwargs])
        elif all([kwargs[k] is None for k in kwargs]) and input_data is None:
            self._from_dict({})
        elif input_data:
            raise TypeError("Invalid parameter 'input_data'")
        else:
            par_list = [f"'{k}'" for k in kwargs]
            par_string = ", ".join(par_list[:-1])
            par_string += f" or {par_list[-1]}"
            raise TypeError(f"Invalid parameters {par_string}")
        self.set_cxx_object()

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
            return self.base.activationEnergy()

    property allow_negative_pre_exponential_factor:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        def __get__(self):
            return self.base.allowNegativePreExponentialFactor()
        def __set__(self, cbool allow):
            self.base.setAllowNegativePreExponentialFactor(allow)


cdef class ArrheniusRate(ArrheniusRateBase):
    r"""
    A reaction rate coefficient which depends on temperature only and follows
    the modified Arrhenius form:

    .. math::

        k_f = A T^b \exp(-\tfrac{E_a}{RT})

    where ``A`` is the
    `pre_exponential_factor <ArrheniusRateBase.pre_exponential_factor>`, ``b`` is the
    `temperature_exponent <ArrheniusRateBase.temperature_exponent>`, and ``Ea`` is the
    `activation_energy <ArrheniusRateBase.activation_energy>`.
    """
    _reaction_rate_type = "Arrhenius"

    def __cinit__(self, A=None, b=None, Ea=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea=Ea)

    def _from_dict(self, input_data):
        self._rate.reset(new CxxArrheniusRate(py_to_anymap(input_data)))

    def _from_parameters(self, A, b, Ea):
        self._rate.reset(new CxxArrheniusRate(A, b, Ea))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxArrheniusRate*>self.rate

    cdef CxxArrheniusRate* cxx_object(self):
        return <CxxArrheniusRate*>self.rate


cdef class BlowersMaselRate(ArrheniusRateBase):
    r"""
    A reaction rate coefficient which depends on temperature and enthalpy change
    of the reaction follows the Blowers-Masel approximation and modified Arrhenius form
    described in `ArrheniusRate`.
    """
    _reaction_rate_type = "Blowers-Masel"

    def __cinit__(self, A=None, b=None, Ea0=None, w=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea0=Ea0, w=w)

    def _from_dict(self, input_data):
        self._rate.reset(new CxxBlowersMaselRate(py_to_anymap(input_data)))

    def _from_parameters(self, A, b, Ea0, w):
        self._rate.reset(new CxxBlowersMaselRate(A, b, Ea0, w))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxArrheniusBase*>self.rate

    cdef CxxBlowersMaselRate* cxx_object(self):
        return <CxxBlowersMaselRate*>self.rate

    property bond_energy:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        def __get__(self):
            return self.cxx_object().bondEnergy()

    property delta_enthalpy:
        """
        Enthalpy change of reaction ``deltaH`` [J/kmol]

        The enthalpy change of reaction is a function of temperature and thus not
        an independent property. Accordingly, the setter should be only used for
        testing purposes, as any value will be overwritten by an update of the
        thermodynamic state.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        def __get__(self):
            return self.cxx_object().deltaH()
        def __set__(self, double delta_H):
            self.cxx_object().setDeltaH(delta_H)


cdef class TwoTempPlasmaRate(ArrheniusRateBase):
    r"""
    A reaction rate coefficient which depends on both gas and electron temperature
    with the form similar to the modified Arrhenius form. Specifically, the temperature
    exponent (b) is applied to the electron temperature instead. In addition, the
    exponential term with activation energy for electron is included.

    .. math::

        k_f = A T_e^b \exp(-\tfrac{E_{a,g}}{RT}) \exp(\tfrac{E_{a,e}(T_e - T)}{R T T_e})

    where :math:`A` is the
    `pre_exponential_factor <ArrheniusRateBase.pre_exponential_factor>`,
    :math:`b` is the `temperature_exponent <ArrheniusRateBase.temperature_exponent>`,
    :math:`E_{a,g}` (``Ea_gas``) is the
    `activation_energy <ArrheniusRateBase.activation_energy>`, and
    :math:`E_{a,e}` (``Ea_electron``) is the `activation_electron_energy`.
    """
    _reaction_rate_type = "two-temperature-plasma"

    def __cinit__(self, A=None, b=None, Ea_gas=0.0, Ea_electron=0.0,
            input_data=None, init=True):
        if init:
            if A is None and b is None:
                Ea_gas = None
                Ea_electron = None
            self._cinit(input_data, A=A, b=b, Ea_gas=Ea_gas, Ea_electron=Ea_electron)

    def __call__(self, double temperature, double elec_temp):
        """
        Evaluate rate expression based on temperature and enthalpy change of reaction.
        """
        return self.rate.eval(temperature, elec_temp)

    def _from_dict(self, input_data):
        self._rate.reset(
            new CxxTwoTempPlasmaRate(py_to_anymap(input_data, hyphenize=True))
        )

    def _from_parameters(self, A, b, Ea_gas, Ea_electron):
        self._rate.reset(new CxxTwoTempPlasmaRate(A, b, Ea_gas, Ea_electron))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxArrheniusBase*>self.rate

    cdef CxxTwoTempPlasmaRate* cxx_object(self):
        return <CxxTwoTempPlasmaRate*>self.rate

    property activation_electron_energy:
        """
        The activation electron energy :math:`E_{a,e}` [J/kmol].
        """
        def __get__(self):
            return self.cxx_object().activationElectronEnergy()


cdef class ElectronCollisionPlasmaRate(ReactionRate):
    r"""
    A reaction rate coefficient which depends on electron collision cross section
    and electron energy distribution

    """
    _reaction_rate_type = "electron-collision-plasma"

    def __cinit__(self, energy_levels=None, cross_sections=None,
                  input_data=None, init=True):
        if init:
            if isinstance(input_data, dict):
                self._from_dict(input_data)
            elif energy_levels is not None and cross_sections is not None:
                self._from_parameters(energy_levels, cross_sections)
            elif input_data is None:
                self._from_dict({})
            else:
                raise TypeError("Invalid input parameters")
            self.set_cxx_object()

    def _from_dict(self, dict input_data):
        self._rate.reset(
            new CxxElectronCollisionPlasmaRate(py_to_anymap(input_data, hyphenize=True))
        )

    def _from_parameters(self, energy_levels, cross_sections):
        # check length
        if len(energy_levels) != len(cross_sections):
            raise ValueError('Length of energy levels and '
                             'cross sections are different')
        cdef np.ndarray[np.double_t, ndim=1] data_energy_levels = \
            np.ascontiguousarray(energy_levels, dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] data_cross_sections = \
            np.ascontiguousarray(cross_sections, dtype=np.double)
        input_data = {'type': 'electron-collision-plasma',
                      'energy-levels': data_energy_levels,
                      'cross-sections': data_cross_sections}
        self._from_dict(input_data)

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxElectronCollisionPlasmaRate*>self.rate

    property energy_levels:
        """
        The energy levels [eV]. Each level corresponds to a cross section
        of `cross_sections`.
        """
        def __get__(self):
            return np.fromiter(self.base.energyLevels(), np.double)

    property cross_sections:
        """
        The cross sections [m2]. Each cross section corresponds to a energy
        level of `energy_levels`.
        """
        def __get__(self):
            return np.fromiter(self.base.crossSections(), np.double)


cdef class FalloffRate(ReactionRate):
    """
    Base class for parameterizations used to describe the fall-off in reaction rates
    due to intermolecular energy transfer.

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
                self._from_dict(input_data)
            elif input_data is None:
                self._from_dict({})
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
            self.falloff.setLowRate(deref(rate.base))

    property high_rate:
        """ Get/Set the `Arrhenius` rate constant in the high-pressure limit """
        def __get__(self):
            return Arrhenius.wrap(&(self.falloff.highRate()))
        def __set__(self, Arrhenius rate):
            self.falloff.setHighRate(deref(rate.base))

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

    def falloff_function(self, double temperature, double conc3b):
        """
        Evaluate the falloff function based on temperature and third-body
        concentration.
        """
        return self.falloff.evalF(temperature, conc3b)


cdef class LindemannRate(FalloffRate):
    r"""
    The Lindemann falloff parameterization.

    This class implements the simple falloff function :math:`F(T,P_r) = 1.0`.
    """
    _reaction_rate_type = "Lindemann"

    def _from_dict(self, input_data):
        self._rate.reset(
            new CxxLindemannRate(py_to_anymap(input_data, hyphenize=True))
        )

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

    def _from_dict(self, input_data):
        self._rate.reset(
            new CxxTroeRate(py_to_anymap(input_data, hyphenize=True))
        )

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

    def _from_dict(self, input_data):
        self._rate.reset(
            new CxxSriRate(py_to_anymap(input_data, hyphenize=True))
        )

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = <CxxFalloffRate*>self.rate


cdef class TsangRate(FalloffRate):
    r"""
    The Tsang falloff parameterization.
    """
    _reaction_rate_type = "Tsang"

    def _from_dict(self, input_data):
        self._rate.reset(
            new CxxTsangRate(py_to_anymap(input_data, hyphenize=True))
        )

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
                self._rate.reset(new CxxPlogRate(py_to_anymap(input_data)))
            elif rates is None:
                self._rate.reset(new CxxPlogRate(py_to_anymap({})))
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
            cdef multimap[double, CxxArrheniusRate] cxxrates
            cdef pair[double, CxxArrheniusRate] p_rate
            cxxrates = self.cxx_object().getRates()
            for p_rate in cxxrates:
                rates.append((p_rate.first, copyArrhenius(&p_rate.second)))
            return rates

        def __set__(self, rates):
            cdef multimap[double, CxxArrheniusRate] ratemap
            cdef Arrhenius rate
            cdef pair[double, CxxArrheniusRate] item
            for p, rate in rates:
                item.first = p
                item.second = deref(rate.base)
                ratemap.insert(item)

            self._rate.reset(new CxxPlogRate(ratemap))
            self.rate = self._rate.get()

cdef class LmrRate(ReactionRate):
    r"""
    A pressure-dependent reaction rate parameterized by logarithmically
    interpolating between Arrhenius rate expressions at various pressures.
    """
    _reaction_rate_type = "LMR_R"

    def __cinit__(self, rates=None, input_data=None, init=True):

        if init and isinstance(rates, list):
            self.rates = rates

        elif init:
            if isinstance(input_data, dict):
                self._rate.reset(new CxxLmrRate(py_to_anymap(input_data)))
            elif rates is None:
                self._rate.reset(new CxxLmrRate(py_to_anymap({})))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameter 'rates'")
            self.set_cxx_object()

    cdef CxxLmrRate* cxx_object(self):
        return <CxxLmrRate*>self.rate


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
                self._rate.reset(new CxxChebyshevRate(py_to_anymap(input_data)))
            elif all([arg is not None
                    for arg in [temperature_range, pressure_range, data]]):
                Tmin = temperature_range[0]
                Tmax = temperature_range[1]
                Pmin = pressure_range[0]
                Pmax = pressure_range[1]
                self._rate.reset(
                    new CxxChebyshevRate(Tmin, Tmax, Pmin, Pmax, self._cxxarray2d(data)))
            elif all([arg is None
                    for arg in [temperature_range, pressure_range, data, input_data]]):
                self._rate.reset(new CxxChebyshevRate(py_to_anymap({})))
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

    cdef CxxChebyshevRate* cxx_object(self):
        return <CxxChebyshevRate*>self.rate

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
            return self.cxx_object().nTemperature()

    property n_pressure:
        """
        Number of pressures over which the Chebyshev fit is computed
        (same as number of columns of `data` property).
        """
        def __get__(self):
            return self.cxx_object().nPressure()

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

    .. warning::

        This class is an experimental part of the Cantera API and
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
                raise TypeError(
                    f"Cannot convert input with type '{type(k)}' to rate expression.")

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


cdef class ExtensibleRate(ReactionRate):
    """
    A base class for a user-defined reaction rate parameterization. Classes derived from
    this class should be decorated with the `extension` decorator to specify the name
    of the rate parameterization and its corresponding data class, and to make these
    rates constructible through factory functions and input files.

    Classes derived from `ExtensibleRate` should implement the `set_parameters`,
    `get_parameters`, `eval`, and (optionally) `validate` methods, which will be called
    as delegates from the C++ :ct:`ReactionRate` class.

    .. warning::

        The delegatable methods defined here are an experimental part of the
        Cantera API and may change without notice.

    .. versionadded:: 3.0
    """

    _reaction_rate_type = "extensible"

    delegatable_methods = {
        "eval": ("evalFromStruct", "double(void*)", "replace"),
        "set_parameters": ("setParameters", "void(AnyMap&, UnitStack&)", "after"),
        "get_parameters": ("getParameters", "void(AnyMap&)", "replace"),
        "validate": ("validate", "void(string, void*)", "replace")
    }

    def __cinit__(self, *args, init=True, **kwargs):
        if init:
            self.set_cxx_object()

    def __init__(self, *args, input_data=None, **kwargs):
        # ReactionRate does not define __init__, so it does not need to be called
        if input_data is not None:
            self.set_parameters(input_data, UnitStack())

    def set_parameters(self, params: AnyMap, rate_coeff_units: UnitStack) -> None:
        """
        Responsible for setting rate parameters based on the input data. For example,
        for reactions created from YAML, ``params`` is the YAML reaction entry converted
        to an ``AnyMap``. ``rate_coeff_units`` specifies the units of the rate
        coefficient.

        Input values contained in ``params`` may be in non-default unit systems,
        specified in the user-provided input file. To convert them to Cantera's native
        mks+kmol unit system, use the functions `AnyMap.convert`,
        `AnyMap.convert_activation_energy`, and `AnyMap.convert_rate_coeff` as needed.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.set_parameters")

    def get_parameters(self, params: AnyMap) -> None:
        """
        Responsible for serializing the state of the ExtensibleRate object, using the
        same format as a YAML reaction entry. This is the inverse of `set_parameters`.

        Serialization methods may request output in unit systems other than Cantera's
        native mks+kmol system. To enable conversions to the user-specified unit system,
        dimensional values should be added to ``params`` using the methods
        `AnyMap.set_quantity` and `AnyMap.set_activation_energy`.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.get_parameters")

    def eval(self, data: ExtensibleRateData) -> float:
        """
        Responsible for calculating the forward rate constant based on the current state
        of the phase, stored in an instance of a class derived from
        `ExtensibleRateData`.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.eval")

    def validate(self, equation: str, soln: "Solution") -> None:
        """
        Responsible for validating that the rate expression is configured with valid
        parameters. This may depend on properties of the Solution, for example
        temperature ranges over which the rate expression can be evaluated. Raises an
        exception if any validation fails.
        """
        pass

    cdef set_cxx_object(self, CxxReactionRate* rate=NULL):
        cdef CxxDelegator* drate
        cdef shared_ptr[CxxExternalHandle] handle

        if rate is NULL:
            # Started with Python object first. Create the C++ object and attach to it.
            # In this case, the Python object owns the C++ object, via self._rate
            self._rate.reset(new CxxReactionRateDelegator())
            self.rate = self._rate.get()
            drate = dynamic_cast[CxxDelegatorPtr](self.rate)
            handle.reset(new CxxPythonHandle(<PyObject*>self, True))
            drate.holdExternalHandle(stringify('python'), handle)
        else:
            # Set up Python object from a C++ object that was created first. In this
            # case, the C++ object owns the Python object, and self._rate is empty to
            # avoid creating a circular dependency.
            self._rate.reset()
            self.rate = rate

        assign_delegates(self, dynamic_cast[CxxDelegatorPtr](self.rate))
        (<CxxReactionRateDelegator*>self.rate).setType(
            stringify(self._reaction_rate_type))


cdef class ExtensibleRateData:
    """
    A base class for data used when evaluating instances of `ExtensibleRate`. Classes
    derived from `ExtensibleRateData` are used to store common data needed to evaluate
    all reactions of a particular type.

    Classes derived from `ExtensibleRateData` must implement the `update` method. After
    the `update` method has been called, instances of `ExtensibleRateData` are passed as
    the argument to `ExtensibleRate.eval`.

    .. versionadded:: 3.0
    """
    delegatable_methods = {
        "update": ("update", "double(void*)", "replace")
    }

    def update(self, soln):
        """
        This method takes a `Solution` object and stores any thermodynamic data (for
        example, temperature and pressure) needed to evaluate all reactions of the
        corresponding ReactionRate type.

        If this state data has changed since the last time `update` was called and the
        reaction rates need to be updated, this method should return `True`. Otherwise,
        it should return `False`.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.update")

    cdef set_cxx_object(self, CxxReactionDataDelegator* data):
        assign_delegates(self, dynamic_cast[CxxDelegatorPtr](data))


cdef class InterfaceRateBase(ArrheniusRateBase):
    """
    Base class collecting commonly used features of Arrhenius-type rate objects
    that include coverage dependencies.
    """

    def __call__(self, double temperature, np.ndarray coverages):
        """
        Evaluate rate expression based on temperature and surface coverages.

        .. warning::

            This method is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        cdef vector[double] cxxdata
        for c in coverages:
            cxxdata.push_back(c)
        return self.rate.eval(temperature, cxxdata)

    property coverage_dependencies:
        """
        Get/set a dictionary containing adjustments to the Arrhenius rate expression
        dependent on surface species coverages. The keys of the dictionary are species
        names, and the values are dictionaries specifying the three coverage
        parameters ``a``, ``m`` and ``E`` which are the modifiers for the pre-exponential
        factor [m, kmol, s units], the temperature exponent [nondimensional],
        and the activation energy [J/kmol], respectively.
        """
        def __get__(self):
            cdef CxxAnyMap cxx_deps
            self.interface.getCoverageDependencies(cxx_deps)
            return anymap_to_py(cxx_deps)
        def __set__(self, deps):
            cdef CxxAnyMap cxx_deps = py_to_anymap(deps)

            self.interface.setCoverageDependencies(cxx_deps)

    def set_species(self, species):
        """
        Set association with an ordered list of all species associated with an
        `InterfaceKinetics` object.
        """
        cdef vector[string] cxxvector
        for s in species:
            cxxvector.push_back(stringify(s))
        self.interface.setSpecies(cxxvector)

    property site_density:
        """
        Site density [kmol/m^2]

        The site density is not an independent property, as it is set by an associated
        `InterfaceKinetics` object. Accordingly, the setter should be only used for
        testing purposes, as the value will be overwritten by an update of the
        thermodynamic state.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        def __get__(self):
            return self.interface.siteDensity()
        def __set__(self, double site_density):
            self.interface.setSiteDensity(site_density)

    property uses_electrochemistry:
        """
        Return boolean flag indicating whether rate involves a charge transfer.
        """
        def __get__(self):
            return self.interface.usesElectrochemistry()

    property beta:
        """
        Return the charge transfer beta parameter
        """
        def __get__(self):
            return self.interface.beta()


cdef class InterfaceArrheniusRate(InterfaceRateBase):
    r"""
    A reaction rate coefficient which depends on temperature and interface coverage
    """
    _reaction_rate_type = "interface-Arrhenius"

    def __cinit__(self, A=None, b=None, Ea=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea=Ea)

    def _from_dict(self, input_data):
        self._rate.reset(new CxxInterfaceArrheniusRate(py_to_anymap(input_data)))

    def _from_parameters(self, A, b, Ea):
        self._rate.reset(new CxxInterfaceArrheniusRate(A, b, Ea))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxArrheniusRate*>self.rate
        self.interface = <CxxInterfaceRateBase*>self.cxx_object()

    cdef CxxInterfaceArrheniusRate* cxx_object(self):
        return <CxxInterfaceArrheniusRate*>self.rate


cdef class InterfaceBlowersMaselRate(InterfaceRateBase):
    r"""
    A reaction rate coefficient which depends on temperature and enthalpy change
    of the reaction follows the Blowers-Masel approximation and modified Arrhenius form
    described in `ArrheniusRate`.
    """
    _reaction_rate_type = "interface-Blowers-Masel"

    def __cinit__(self, A=None, b=None, Ea0=None, w=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea0=Ea0, w=w)

    def _from_dict(self, input_data):
        self._rate.reset(new CxxInterfaceBlowersMaselRate(py_to_anymap(input_data)))

    def _from_parameters(self, A, b, Ea0, w):
        self._rate.reset(new CxxInterfaceBlowersMaselRate(A, b, Ea0, w))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxArrheniusRate*>self.rate
        self.interface = <CxxInterfaceRateBase*>self.cxx_object()

    cdef CxxInterfaceBlowersMaselRate* cxx_object(self):
        return <CxxInterfaceBlowersMaselRate*>self.rate

    property bond_energy:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        def __get__(self):
            return self.cxx_object().bondEnergy()

    property delta_enthalpy:
        """
        Enthalpy change of reaction ``deltaH`` [J/kmol]

        The enthalpy change of reaction is a function of temperature and thus not
        an independent property. Accordingly, the setter should be only used for
        testing purposes, as any value will be overwritten by an update of the
        thermodynamic state.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        def __get__(self):
            return self.cxx_object().deltaH()
        def __set__(self, double delta_H):
            self.cxx_object().setDeltaH(delta_H)


cdef class StickRateBase(InterfaceRateBase):
    """
    Base class collecting commonly used features of Arrhenius-type sticking rate
    objects that include coverage dependencies.
    """

    property motz_wise_correction:
        """
        Get/Set a boolean indicating whether to use the correction factor developed by
        Motz & Wise for reactions with high (near-unity) sticking coefficients when
        converting the sticking coefficient to a rate coefficient.
        """
        def __get__(self):
            return self.stick.motzWiseCorrection()
        def __set__(self, cbool motz_wise):
            self.stick.setMotzWiseCorrection(motz_wise)

    property sticking_species:
        """
        The name of the sticking species. Needed only for reactions with
        multiple non-surface reactant species, where the sticking species is
        ambiguous.
        """
        def __get__(self):
            return pystr(self.stick.stickingSpecies())
        def __set__(self, species):
            self.stick.setStickingSpecies(stringify(species))

    property sticking_order:
        """
        The exponent applied to site density (sticking order).

        The sticking order is not an independent property and is detected automatically
        by Cantera. Accordingly, the setter should be only used for testing purposes.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        def __get__(self):
            return self.stick.stickingOrder()
        def __set__(self, double order):
            self.stick.setStickingOrder(order)

    property sticking_weight:
        """
        The molecular weight of the sticking species.

        The sticking weight is not an independent property and is detected automatically
        by Cantera. Accordingly, the setter should be only used for testing purposes.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        def __get__(self):
            return self.stick.stickingWeight()
        def __set__(self, double weight):
            self.stick.setStickingWeight(weight)


cdef class StickingArrheniusRate(StickRateBase):
    r"""
    A surface sticking rate expression based on the Arrhenius parameterization
    """
    _reaction_rate_type = "sticking-Arrhenius"

    def __cinit__(self, A=None, b=None, Ea=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea=Ea)

    def _from_dict(self, input_data):
        self._rate.reset(new CxxStickingArrheniusRate(py_to_anymap(input_data)))

    def _from_parameters(self, A, b, Ea):
        self._rate.reset(new CxxStickingArrheniusRate(A, b, Ea))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxArrheniusRate*>self.rate
        self.stick = <CxxStickingCoverage*>self.cxx_object()
        self.interface = <CxxInterfaceRateBase*>self.stick

    cdef CxxStickingArrheniusRate* cxx_object(self):
        return <CxxStickingArrheniusRate*>self.rate


cdef class StickingBlowersMaselRate(StickRateBase):
    r"""
    A surface sticking rate expression based on the Blowers-Masel parameterization
    """
    _reaction_rate_type = "sticking-Blowers-Masel"

    def __cinit__(self, A=None, b=None, Ea0=None, w=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea0=Ea0, w=w)

    def _from_dict(self, input_data):
        self._rate.reset(new CxxStickingBlowersMaselRate(py_to_anymap(input_data)))

    def _from_parameters(self, A, b, Ea0, w):
        self._rate.reset(new CxxStickingBlowersMaselRate(A, b, Ea0, w))

    cdef set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = <CxxArrheniusRate*>self.rate
        self.stick = <CxxStickingCoverage*>self.cxx_object()
        self.interface = <CxxInterfaceRateBase*>self.stick

    cdef CxxStickingBlowersMaselRate* cxx_object(self):
        return <CxxStickingBlowersMaselRate*>self.rate

    property bond_energy:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        def __get__(self):
            return self.cxx_object().bondEnergy()

    property delta_enthalpy:
        """
        Enthalpy change of reaction ``deltaH`` [J/kmol]

        The enthalpy change of reaction is a function of temperature and thus not
        an independent property. Accordingly, the setter should be only used for
        testing purposes, as any value will be overwritten by an update of the
        thermodynamic state.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        def __get__(self):
            return self.cxx_object().deltaH()
        def __set__(self, double delta_H):
            self.cxx_object().setDeltaH(delta_H)


cdef class ThirdBody:
    r"""
    Class representing third-body collision partners in three-body or falloff reactions.

    :param collider:
        Name of the third-body collider. If ``M`` (default), the `default_efficiency`
        is set to 1 and the collider is assumed to participate in a three-body reaction.
        If the collider includes parentheses, - for example ``(+M)``, - a falloff form
        is assumed, where the collider is not considered for the law of mass action.
        For species other than ``M``, the third-body collider represents a named species
        with collision efficiency 1, and the `default_efficiency` is set to zero.
    :param efficiencies:
        Non-default third-body efficiencies
    :param default_efficiency:
        Default third-body efficiency

    .. versionadded:: 3.0
    """
    def __cinit__(self, collider="M", *,
                  efficiencies=None, default_efficiency=None, init=True):
        if not init:
            return
        self._third_body.reset(new CxxThirdBody(stringify(collider)))
        self.third_body = self._third_body.get()

        if efficiencies is not None:
            self.efficiencies = efficiencies

        if default_efficiency is not None:
            self.default_efficiency = default_efficiency

    @staticmethod
    cdef wrap(shared_ptr[CxxThirdBody] third_body):
        tb = ThirdBody(init=False)
        tb._third_body = third_body
        tb.third_body = tb._third_body.get()
        return tb

    property name:
        """
        Get the name of the third-body collider used in the reaction equation.
        """
        def __get__(self):
            return pystr(self.third_body.name())

    property mass_action:
        """
        Retrieve flag indicating whether third-body collider participates
        in the law of mass action.
        """
        def __get__(self):
            return self.third_body.mass_action

    property efficiencies:
        """
        Get/Set a `dict` defining non-default third-body efficiencies for this reaction,
        where the keys are the species names and the values are the efficiencies.
        """
        def __get__(self):
            return comp_map_to_dict(self.third_body.efficiencies)
        def __set__(self, eff):
            self.third_body.efficiencies = comp_map(eff)

    property default_efficiency:
        """
        Get/Set the default third-body efficiency for this reaction, used for species
        not in `efficiencies`.
        """
        def __get__(self):
            return self.third_body.default_efficiency
        def __set__(self, default_eff):
            self.third_body.default_efficiency = default_eff

    def efficiency(self, species):
        """
        Get the efficiency of the third body named ``species`` considering both
        the default efficiency and species-specific efficiencies.
        """
        return self.third_body.efficiency(stringify(species))


cdef class Reaction:
    """
    A class which stores data about a reaction and its rate parameterization so
    that it can be added to a `Kinetics` object.

    :param reactants:
        Value used to set `reactants`
    :param products:
        Value used to set `products`
    :param rate:
        The rate parameterization for the reaction, given as one of the following:

           - a `ReactionRate` object
           - a `dict` containing the parameters needed to construct a `ReactionRate`
             object, with keys corresponding to the YAML format
           - a `dict` containing Arrhenius parameters (``A``, ``b``, and ``Ea``)
    :param equation:
        The reaction equation, used to set the reactants and products if values for
        those arguments are not provided.

    Examples::

        R = ct.Reaction({"O": 1, "H2": 1}, {"H": 1, "OH": 1},
                        ct.ArrheniusRate(38.7, 2.7, 26191840.0))
        R = ct.Reaction(equation="O + H2 <=> H + OH",
                        rate={"A": 38.7, "b": 2.7, "Ea": 26191840.0})
        R = ct.Reaction(equation="HO2 <=> OH + O", rate=ChebyshevRate(...))

    The static method `list_from_file` can be used to create a list of `Reaction`
    objects from existing definitions in the YAML format. The following will produce a
    list of the 325 reactions which make up the GRI 3.0 mechanism::

        R = ct.Reaction.list_from_file("gri30.yaml", gas)

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

    The method `from_yaml` can be used to create individual `Reaction` objects from
    definitions in the YAML format. It is important to verify that either the
    pre-exponential factor and activation energy are supplied in SI units, or
    that they have their units specified::

        R = ct.Reaction.from_yaml('''{equation: O + H2 <=> H + OH,
                rate-constant: {A: 3.87e+04 cm^3/mol/s, b: 2.7, Ea: 6260 cal/mol}}''',
                gas)
    """

    def __cinit__(self, reactants=None, products=None, rate=None, *,
                  equation=None, init=True, third_body=None):
        if not init:
            return

        cdef ReactionRate _rate
        if isinstance(rate, ReactionRate):
            _rate = rate
        elif rate is None:
            # default to Arrhenius expression
            raise ValueError("Missing reaction rate information.")
        elif isinstance(rate, dict):
            if {"A", "b"} - set(rate) == set():
                # Allow simple syntax for Arrhenius-type rates
                args = {"rate-constant": rate}
                if "Ea0" in rate:
                    args.update({"type": "Blowers-Masel"})
                elif "Ea_gas" in rate:
                    args.update({"type": "two-temperature-plasma"})
                _rate = ReactionRate.from_dict(args)
            else:
                _rate = ReactionRate.from_dict(rate)
        elif isinstance(rate, Arrhenius):
            _rate = ArrheniusRate(
                A=rate.pre_exponential_factor,
                b=rate.temperature_exponent,
                Ea=rate.activation_energy
            )
        elif callable(rate):
            _rate = CustomRate(rate)
        else:
            raise TypeError(f"Invalid rate definition with type '{type(rate)}'")

        cdef ThirdBody _third_body
        if isinstance(third_body, ThirdBody):
            _third_body = third_body
        elif isinstance(third_body, str):
            _third_body = ThirdBody(third_body)

        if reactants and products:
            # create from reactant and product compositions
            if third_body:
                self._reaction.reset(
                    new CxxReaction(comp_map(reactants), comp_map(products),
                    _rate._rate, _third_body._third_body)
                )
            else:
                self._reaction.reset(
                    new CxxReaction(comp_map(reactants), comp_map(products),
                    _rate._rate)
                )
        elif equation:
            # create from reaction equation
            if third_body:
                self._reaction.reset(
                    new CxxReaction(stringify(equation),
                    _rate._rate, _third_body._third_body)
                )
            else:
                self._reaction.reset(
                    new CxxReaction(stringify(equation),
                    _rate._rate)
                )
        else:
            # create default object
            raise ValueError("Missing reactant and/or product information.")

        self.reaction = self._reaction.get()
        self._rate = _rate

    @staticmethod
    cdef wrap(shared_ptr[CxxReaction] reaction):
        """
        Wrap a C++ Reaction object with a Python object of the correct derived type.
        """
        # wrap C++ reaction
        cdef Reaction R
        R = Reaction(init=False)
        R._reaction = reaction
        R.reaction = R._reaction.get()
        R._rate = ReactionRate.wrap(R.reaction.rate())
        return R

    @classmethod
    def from_dict(cls, data, Kinetics kinetics, hyphenize=True):
        """
        Create a `Reaction` object from a dictionary corresponding to its YAML
        representation. By default, underscores in keys are replaced by hyphens.

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
        cdef CxxAnyMap any_map = py_to_anymap(data, hyphenize=hyphenize)
        cxx_reaction = CxxNewReaction(any_map, deref(kinetics.kinetics))
        return Reaction.wrap(cxx_reaction)

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
        cdef CxxAnyMap any_map = AnyMapFromYamlString(stringify(text))
        cxx_reaction = CxxNewReaction(any_map, deref(kinetics.kinetics))
        return Reaction.wrap(cxx_reaction)

    @staticmethod
    def list_from_file(filename, Kinetics kinetics, section="reactions"):
        """
        Create a list of Reaction objects from all of the reactions defined in a
        YAML file. Reactions from the section ``section`` will be returned.

        Directories on Cantera's input file path will be searched for the
        specified file.
        """
        root = AnyMapFromYamlFile(stringify(str(filename)))
        cxx_reactions = CxxGetReactions(root[stringify(section)],
                                        deref(kinetics.kinetics))
        return [Reaction.wrap(r) for r in cxx_reactions]

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
        species names and the values, are the stoichiometric coefficients, for example
        ``{'CH4':1, 'OH':1}``, or as a composition string, for example
        ``'CH4:1, OH:1'``.

        .. versionchanged:: 3.1  This is a read-only property
        """
        def __get__(self):
            return comp_map_to_dict(self.reaction.reactants)

    property products:
        """
        Get/Set the products in this reaction as a dict where the keys are
        species names and the values, are the stoichiometric coefficients, for example
        ``{'CH3':1, 'H2O':1}``, or as a composition string, for example
        ``'CH3:1, H2O:1'``.

        .. versionchanged:: 3.1  This is a read-only property
        """
        def __get__(self):
            return comp_map_to_dict(self.reaction.products)

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

    property rate:
        """ Get/Set the reaction rate evaluator for this reaction. """
        def __get__(self):
            return self._rate

        def __set__(self, rate):
            if isinstance(rate, ReactionRate):
                self._rate = rate
            elif callable(rate):
                self._rate = CustomRate(rate)
            else:
                raise TypeError(f"Invalid rate definition with type '{type(rate)}'")
            self.reaction.setRate(self._rate._rate)

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
            return anymap_to_py(self.reaction.parameters(True))

    def update_user_data(self, data):
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.reaction.input.update(py_to_anymap(data), False)

    def clear_user_data(self):
        """
        Clear all saved input data, so that the data given by `input_data` or
        `Solution.write_yaml` will only include values generated by Cantera based on
        the current object state.
        """
        self.reaction.input.clear()

    def __repr__(self):
        return f"{self.equation}    <{self.reaction_type}>"

    def __str__(self):
        return self.equation

    property rate_coeff_units:
        """Get reaction rate coefficient units"""
        def __get__(self):
            cdef CxxUnits rate_units = self.reaction.rate_units
            return Units.copy(rate_units)

    property third_body:
        """
        Returns a `ThirdBody` object if `Reaction` uses a third body collider, and
        ``None`` otherwise.

        .. versionadded:: 3.0
        """
        def __get__(self):
            if self.reaction.usesThirdBody():
                return ThirdBody.wrap(self.reaction.thirdBody())

    @property
    def third_body_name(self):
        """
        Returns name of `ThirdBody` collider if `Reaction` uses a third body collider,
        and ``None`` otherwise.

        .. versionadded:: 3.0
        """
        if self.reaction.usesThirdBody():
            return ThirdBody.wrap(self.reaction.thirdBody()).name
        return None


cdef class Arrhenius:
    r"""
    A reaction rate coefficient which depends on temperature only and follows
    the modified Arrhenius form:

    .. math::

        k_f = A T^b \exp(-\tfrac{E}{RT})

    where ``A`` is the `pre_exponential_factor`, ``b`` is the `temperature_exponent`,
    and ``E`` is the `activation_energy`.
    """
    def __cinit__(self, A=0, b=0, E=0, init=True):
        if init:
            self.base = new CxxArrheniusRate(A, b, E)
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
            return self.base.activationEnergy()

    def __repr__(self):
        return 'Arrhenius(A={:g}, b={:g}, E={:g})'.format(
            self.pre_exponential_factor, self.temperature_exponent,
            self.activation_energy)

    def __call__(self, float T):
        return self.base.evalRate(np.log(T), 1/T)


cdef copyArrhenius(CxxArrheniusRate* rate):
    r = Arrhenius(rate.preExponentialFactor(), rate.temperatureExponent(),
                  rate.activationEnergy())
    return r
