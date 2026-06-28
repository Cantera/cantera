# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3

import numpy as np
from collections.abc import (Callable as _Callable, Iterable as _Iterable,
                             Sequence as _Sequence)
from typing import (Any as _Any, ClassVar as _ClassVar, Generic as _Generic,
                    TypeAlias as _TypeAlias, TYPE_CHECKING,
                    TypedDict as _TypedDict, TypeVar as _TypeVar)
from typing_extensions import NotRequired as _NotRequired
import warnings

if TYPE_CHECKING:
    from .composite import Solution as _Solution
    from .kinetics import Kinetics as _Kinetics

import cython
import cython.cimports.numpy as cnp  # Required: triggers import_array() for numpy C-API

from cython.cimports.cantera.kinetics import Kinetics
from cython.cimports.cantera._utils import (
    stringify, pystr, anymap_to_py, py_to_anymap, comp_map, comp_map_to_dict,
    AnyMapFromYamlString, AnyMapFromYamlFile)
from cython.cimports.cantera.units import Units, UnitStack
from cython.cimports.cantera.delegator import (
    CxxDelegatorPtr, assign_delegates,
    CxxPythonHandle, CxxExternalHandle)
from ._types import Array as _Array, ArrayLike as _ArrayLike
from ._utils import AnyMap as _AnyMap, CanteraError
from .func1 import _Func1Like
from .units import Units as _Units, UnitStack as _UnitStack


class _ArrheniusParameters(_TypedDict):
    A: float
    b: float
    Ea: float


class _BlowersMaselParameters(_TypedDict):
    A: float
    b: float
    Ea0: float
    w: float


class _TwoTempPlasmaParameters(_TypedDict):
    A: float
    b: float
    Ea_gas: float
    Ea_electron: float


class _ElectronCollisionPlasmaParameters(_TypedDict):
    energy_levels: list[float]
    cross_sections: list[float]


class _PlogParameters(_ArrheniusParameters):
    P: float


_ReactionRateParameters: _TypeAlias = (
    _ArrheniusParameters
    | _BlowersMaselParameters
    | _TwoTempPlasmaParameters
    | _ElectronCollisionPlasmaParameters
    | _PlogParameters
)


class _TroeParameters(_TypedDict):
    A: float
    T3: float
    T1: float
    T2: _NotRequired[float]


class _CoverageParameters(_TypedDict):
    a: float
    m: float
    E: float


_T = _TypeVar("_T")

# Published as ``str`` while preserving the pre-merge runtime behavior of
# accepting string subclasses and path-like objects at selected call sites.
_Str: _TypeAlias = str


class _ReactionRateInput(_TypedDict, _Generic[_T], total=False):
    type: str
    rate_constant: _T
    efficiencies: dict[str, float]
    Troe: _TroeParameters
    coverage_dependencies: dict[str, _CoverageParameters]


class _ReactionInput(_ReactionRateInput[_T]):
    equation: str


class _FalloffRateInput(_TypedDict, total=False):
    type: str
    low_P_rate_constant: _ArrheniusParameters
    high_P_rate_constant: _ArrheniusParameters
    Troe: _TroeParameters


class _PlogRateInput(_TypedDict, total=False):
    type: str
    rate_constants: _Sequence[_PlogParameters]


class _ChebyshevRateInput(_TypedDict, total=False):
    type: str
    pressure_range: _Sequence[float]
    temperature_range: _Sequence[float]
    data: _Array

# dictionary to store reaction rate classes
_reaction_rate_class_registry: dict = {}


@cython.cclass
class ReactionRate:
    """
    Base class for ReactionRate objects.

    ReactionRate objects are used to calculate reaction rates and are associated
    with a Reaction object.
    """
    _reaction_rate_type: _ClassVar[str] = ""

    def __repr__(self):
        return f"<{type(self).__name__} at {id(self):0x}>"

    def __call__(self, temperature: cython.double) -> float:
        """
        Evaluate rate expression based on temperature.
        """
        return self.rate.eval(temperature)

    @property
    def type(self) -> str:
        """ Get the C++ ReactionRate type """
        return pystr(self.rate.type())

    @property
    def sub_type(self) -> str:
        """ Get the C++ ReactionRate sub-type """
        return pystr(self.rate.subType())

    @cython.cfunc
    @staticmethod
    def wrap(rate: shared_ptr[CxxReactionRate]):
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
        drate: cython.pointer(CxxDelegator) = dynamic_cast[CxxDelegatorPtr](rate.get())
        handle: cython.pointer(CxxPythonHandle)

        if drate != cython.NULL:
            handle = dynamic_pointer_cast[CxxPythonHandle, CxxExternalHandle](
                drate.getExternalHandle(stringify("python"))).get()
            if handle != cython.NULL and handle.get() != cython.NULL:
                py_rate = cython.cast(ReactionRate, cython.cast(cython.p_void, handle.get()))
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
        rr: ReactionRate
        rr = cls(init=False)
        rr._rate = rate
        rr.set_cxx_object()
        return rr

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()

    @classmethod
    def from_dict(
        cls,
        data: _ReactionRateInput[_ReactionRateParameters],
        hyphenize: bool = True,
    ) -> "ReactionRate":
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

        any_map: CxxAnyMap = py_to_anymap(data, hyphenize=hyphenize)
        cxx_rate = CxxNewReactionRate(any_map)
        return ReactionRate.wrap(cxx_rate)

    @classmethod
    def from_yaml(cls, text: str) -> "ReactionRate":
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

        any_map: CxxAnyMap
        any_map = AnyMapFromYamlString(stringify(text))
        cxx_rate = CxxNewReactionRate(any_map)
        return ReactionRate.wrap(cxx_rate)

    @property
    def input_data(self) -> _ReactionRateInput[_ReactionRateParameters]:
        """
        Get input data for this reaction rate with its current parameter values.
        """
        return anymap_to_py(self.rate.parameters())

    @property
    def conversion_units(self) -> _Units:
        """
        Get the units for converting the leading term in the reaction rate expression
        to different unit systems.
        """
        return Units.copy(self.rate.conversionUnits())


@cython.cclass
class ArrheniusRateBase(ReactionRate):
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

    @property
    def pre_exponential_factor(self) -> float:
        """
        The pre-exponential factor ``A`` in units of m, kmol, and s raised to
        powers depending on the reaction order.
        """
        return self.base.preExponentialFactor()

    @property
    def temperature_exponent(self) -> float:
        """
        The temperature exponent ``b``.
        """
        return self.base.temperatureExponent()

    @property
    def activation_energy(self) -> float:
        """
        The activation energy ``E`` [J/kmol].
        """
        return self.base.activationEnergy()

    @property
    def allow_negative_pre_exponential_factor(self) -> bool:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        return self.base.allowNegativePreExponentialFactor()

    @allow_negative_pre_exponential_factor.setter
    def allow_negative_pre_exponential_factor(self, allow: bool) -> None:
        self.base.setAllowNegativePreExponentialFactor(allow)


@cython.cclass
class ArrheniusRate(ArrheniusRateBase):
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

    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea: float | None = None,
        input_data: _ReactionRateInput[_ArrheniusParameters] | None = None,
        init: bool = True,
    ) -> None:
        pass

    def _from_dict(self, input_data: _ReactionRateInput[_ArrheniusParameters]) -> None:
        self._rate = make_shared[CxxArrheniusRate](py_to_anymap(input_data))

    def _from_parameters(self, A: cython.double, b: cython.double,
                         Ea: cython.double) -> None:
        self._rate = make_shared[CxxArrheniusRate](A, b, Ea)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxArrheniusRate), self.rate)

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxArrheniusRate):
        return cython.cast(cython.pointer(CxxArrheniusRate), self.rate)


@cython.cclass
class BlowersMaselRate(ArrheniusRateBase):
    r"""
    A reaction rate coefficient which depends on temperature and enthalpy change
    of the reaction follows the Blowers-Masel approximation and modified Arrhenius form
    described in `ArrheniusRate`.
    """
    _reaction_rate_type = "Blowers-Masel"

    def __cinit__(self, A=None, b=None, Ea0=None, w=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea0=Ea0, w=w)

    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea0: float | None = None,
        w: float | None = None,
        input_data: _ReactionRateInput[_BlowersMaselParameters] | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def _from_dict(self, input_data: _ReactionRateInput[_BlowersMaselParameters]) -> None:
        self._rate = make_shared[CxxBlowersMaselRate](py_to_anymap(input_data))

    def _from_parameters(self, A: cython.double, b: cython.double,
                         Ea0: cython.double, w: cython.double) -> None:
        self._rate = make_shared[CxxBlowersMaselRate](A, b, Ea0, w)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxArrheniusBase), self.rate)

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxBlowersMaselRate):
        return cython.cast(cython.pointer(CxxBlowersMaselRate), self.rate)

    @property
    def bond_energy(self) -> float:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        return self.cxx_object().bondEnergy()

    @property
    def delta_enthalpy(self) -> float:
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
        return self.cxx_object().deltaH()

    @delta_enthalpy.setter
    def delta_enthalpy(self, delta_H: float) -> None:
        self.cxx_object().setDeltaH(delta_H)


@cython.cclass
class TwoTempPlasmaRate(ArrheniusRateBase):
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

    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea_gas: float = 0.0,
        Ea_electron: float = 0.0,
        input_data: _ReactionRateInput[_TwoTempPlasmaParameters] | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def __call__(self, temperature: cython.double, elec_temp: cython.double) -> float:
        """
        Evaluate rate expression based on temperature and enthalpy change of reaction.
        """
        return self.rate.eval(temperature, elec_temp)

    def _from_dict(
        self, input_data: _ReactionRateInput[_TwoTempPlasmaParameters]
    ) -> None:
        self._rate = make_shared[CxxTwoTempPlasmaRate](
            py_to_anymap(input_data, hyphenize=True)
        )

    def _from_parameters(self, A: cython.double, b: cython.double,
                         Ea_gas: cython.double,
                         Ea_electron: cython.double) -> None:
        self._rate = make_shared[CxxTwoTempPlasmaRate](A, b, Ea_gas, Ea_electron)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxArrheniusBase), self.rate)

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxTwoTempPlasmaRate):
        return cython.cast(cython.pointer(CxxTwoTempPlasmaRate), self.rate)

    @property
    def activation_electron_energy(self) -> float:
        """
        The activation electron energy :math:`E_{a,e}` [J/kmol].
        """
        return self.cxx_object().activationElectronEnergy()


@cython.cclass
class ElectronCollisionPlasmaRate(ReactionRate):
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

    def __init__(
        self,
        energy_levels: _ArrayLike | None = None,
        cross_sections: _ArrayLike | None = None,
        input_data: _ReactionRateInput[_ElectronCollisionPlasmaParameters] | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def _from_dict(
        self, input_data: _ReactionRateInput[_ElectronCollisionPlasmaParameters]
    ) -> None:
        self._rate = make_shared[CxxElectronCollisionPlasmaRate](
            py_to_anymap(input_data, hyphenize=True)
        )

    def _from_parameters(self, energy_levels: _ArrayLike,
                         cross_sections: _ArrayLike) -> None:
        # check length
        if len(energy_levels) != len(cross_sections):
            raise ValueError('Length of energy levels and '
                             'cross sections are different')
        data_energy_levels = np.ascontiguousarray(energy_levels, dtype=np.double)
        data_cross_sections = np.ascontiguousarray(cross_sections, dtype=np.double)
        input_data = {'type': 'electron-collision-plasma',
                      'energy-levels': data_energy_levels,
                      'cross-sections': data_cross_sections}
        self._from_dict(input_data)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxElectronCollisionPlasmaRate), self.rate)

    @property
    def energy_levels(self) -> _Array:
        """
        The energy levels [eV]. Each level corresponds to a cross section
        of `cross_sections`.
        """
        levels: span[const_double] = self.base.energyLevels()
        data = np.empty(levels.size())
        i: cython.size_t
        for i in range(levels.size()):
            data[i] = levels[i]
        return data

    @property
    def cross_sections(self) -> _Array:
        """
        The cross sections [m2]. Each cross section corresponds to a energy
        level of `energy_levels`.
        """
        cxsections: span[const_double] = self.base.crossSections()
        data = np.empty(cxsections.size())
        i: cython.size_t
        for i in range(cxsections.size()):
            data[i] = cxsections[i]
        return data


@cython.cclass
class FalloffRate(ReactionRate):
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

    def __init__(
        self,
        low: ArrheniusRate | None = None,
        high: ArrheniusRate | None = None,
        falloff_coeffs: _Sequence[float] | None = None,
        input_data: _ReactionRateInput[_FalloffRateInput] | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def __call__(self, temperature: cython.double, concm: cython.double) -> float:
        """
        Evaluate rate expression based on temperature and third-body concentration.
        """
        return self.rate.eval(temperature, concm)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = cython.cast(cython.pointer(CxxFalloffRate), self.rate)

    @property
    def low_rate(self) -> ArrheniusRate:
        """ Get/Set the `ArrheniusRate` rate constant in the low-pressure limit """
        # Copy the C++ ArrheniusRate object that is returned by reference
        low: cython.pointer(CxxArrheniusRate) = cython.address(self.falloff.lowRate())
        rate_ptr: shared_ptr[CxxReactionRate]
        rate_ptr = make_shared[CxxArrheniusRate](low[0])
        return ArrheniusRate.wrap(rate_ptr)

    @low_rate.setter
    def low_rate(self, rate: ArrheniusRate) -> None:
        if isinstance(rate, ArrheniusRate):
            self.falloff.setLowRate(
                cython.cast(cython.pointer(CxxArrheniusRate),
                    cython.cast(ArrheniusRate, rate).rate)[0])
        else:
            raise TypeError("FalloffRate.low_rate setter: Expected 'ArrheniusRate'"
                f"but got '{type(rate).__name__}'")

    @property
    def high_rate(self) -> ArrheniusRate:
        """ Get/Set the `ArrheniusRate` rate constant in the high-pressure limit """
        # Copy the C++ ArrheniusRate object that is returned by reference
        high: cython.pointer(CxxArrheniusRate) = cython.address(self.falloff.highRate())
        rate_ptr: shared_ptr[CxxReactionRate]
        rate_ptr = make_shared[CxxArrheniusRate](high[0])
        return ArrheniusRate.wrap(rate_ptr)

    @high_rate.setter
    def high_rate(self, rate: ArrheniusRate) -> None:
        if isinstance(rate, ArrheniusRate):
            self.falloff.setHighRate(
                cython.cast(cython.pointer(CxxArrheniusRate),
                    cython.cast(ArrheniusRate, rate).rate)[0])
        else:
            raise TypeError("FalloffRate.high_rate setter: Expected 'ArrheniusRate'"
                f"but got '{type(rate).__name__}'")

    @property
    def falloff_coeffs(self) -> _Array:
        """ The array of coefficients used to define this falloff function. """
        n: cython.size_t = self.falloff.nParameters()
        data = np.empty(n)
        if n:
            cdata: cython.double[::1] = data
            self.falloff.getFalloffCoeffs(span[cython.double](
                cython.address(cdata[0]), cython.cast(cython.size_t, n)))
        return data

    @falloff_coeffs.setter
    def falloff_coeffs(self, data: _Iterable[float]) -> None:
        cxxdata: vector[cython.double]
        for c in data:
            cxxdata.push_back(c)
        if cxxdata.size():
            self.falloff.setFalloffCoeffs(
                span[cython.double](cython.address(cxxdata[0]), cxxdata.size()))
        else:
            self.falloff.setFalloffCoeffs(span[cython.double]())

    @property
    def allow_negative_pre_exponential_factor(self) -> bool:
        """
        Get/Set whether the rate coefficient is allowed to have a negative
        pre-exponential factor.
        """
        return self.falloff.allowNegativePreExponentialFactor()

    @allow_negative_pre_exponential_factor.setter
    def allow_negative_pre_exponential_factor(self, allow: bool) -> None:
        self.falloff.setAllowNegativePreExponentialFactor(allow)

    @property
    def chemically_activated(self) -> bool:
        """
        Get whether the object is a chemically-activated reaction rate.
        """
        return self.falloff.chemicallyActivated()

    @chemically_activated.setter
    def chemically_activated(self, activated: bool) -> None:
        self.falloff.setChemicallyActivated(activated)

    def falloff_function(self, temperature: cython.double,
                         conc3b: cython.double) -> float:
        """
        Evaluate the falloff function based on temperature and third-body
        concentration.
        """
        return self.falloff.evalF(temperature, conc3b)


@cython.cclass
class LindemannRate(FalloffRate):
    r"""
    The Lindemann falloff parameterization.

    This class implements the simple falloff function :math:`F(T,P_r) = 1.0`.
    """
    _reaction_rate_type = "Lindemann"

    def _from_dict(self, input_data):
        self._rate = make_shared[CxxLindemannRate](
            py_to_anymap(input_data, hyphenize=True)
        )

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = cython.cast(cython.pointer(CxxFalloffRate), self.rate)


@cython.cclass
class TroeRate(FalloffRate):
    r"""
    The 3- or 4-parameter Troe falloff function.

    :param falloff_coeffs:
        An array of 3 or 4 parameters: :math:`[a, T^{***}, T^*, T^{**}]` where
        the final parameter is optional (with a default value of 0).
    """
    _reaction_rate_type = "Troe"

    def _from_dict(self, input_data):
        self._rate = make_shared[CxxTroeRate](
            py_to_anymap(input_data, hyphenize=True)
        )

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = cython.cast(cython.pointer(CxxFalloffRate), self.rate)


@cython.cclass
class SriRate(FalloffRate):
    r"""
    The 3- or 5-parameter SRI falloff function.

    :param falloff_coeffs:
        An array of 3 or 5 parameters: :math:`[a, b, c, d, e]` where the last
        two parameters are optional (with default values of 1 and 0, respectively).
    """
    _reaction_rate_type = "SRI"

    def _from_dict(self, input_data):
        self._rate = make_shared[CxxSriRate](
            py_to_anymap(input_data, hyphenize=True)
        )

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = cython.cast(cython.pointer(CxxFalloffRate), self.rate)


@cython.cclass
class TsangRate(FalloffRate):
    r"""
    The Tsang falloff parameterization.
    """
    _reaction_rate_type = "Tsang"

    def _from_dict(self, input_data):
        self._rate = make_shared[CxxTsangRate](
            py_to_anymap(input_data, hyphenize=True)
        )

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.falloff = cython.cast(cython.pointer(CxxFalloffRate), self.rate)


@cython.cclass
class PlogRate(ReactionRate):
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
                self._rate = make_shared[CxxPlogRate](py_to_anymap(input_data))
            elif rates is None:
                self._rate = make_shared[CxxPlogRate](py_to_anymap({}))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameter 'rates'")
            self.set_cxx_object()

    def __init__(
        self,
        rates: list[tuple[float, ArrheniusRate]] | None = None,
        input_data: _PlogRateInput | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def __call__(self, temperature: cython.double, pressure: cython.double) -> float:
        """
        Evaluate rate expression based on temperature and pressure.
        """
        return self.rate.eval(temperature, pressure)

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxPlogRate):
        return cython.cast(cython.pointer(CxxPlogRate), self.rate)

    @property
    def rates(self) -> list[tuple[float, ArrheniusRate]]:
        """
        Get/Set the rate coefficients for this reaction, which are given as a
        list of (pressure, `ArrheniusRate`) tuples.
        """
        rates = []
        cxxrates: multimap[cython.double, CxxArrheniusRate]
        p_rate: pair[cython.double, CxxArrheniusRate]
        rate_ptr: shared_ptr[CxxReactionRate]
        cxxrates = self.cxx_object().getRates()
        for p_rate in cxxrates:
            rate_ptr = make_shared[CxxArrheniusRate](p_rate.second)
            rates.append((p_rate.first, ArrheniusRate.wrap(rate_ptr)))
        return rates

    @rates.setter
    def rates(self, rates: _Iterable[tuple[float, ArrheniusRate]]) -> None:
        ratemap: multimap[cython.double, CxxArrheniusRate]
        item: pair[cython.double, CxxArrheniusRate]
        for p, rate in rates:
            item.first = p
            item.second = cython.cast(cython.pointer(CxxArrheniusRate),
                cython.cast(ArrheniusRate, rate).rate)[0]
            ratemap.insert(item)

        self._rate = make_shared[CxxPlogRate](ratemap)
        self.rate = self._rate.get()


@cython.cclass
class LinearBurkeRate(ReactionRate):
    r"""
    A reaction rate dependent on both pressure and mixture composition that accounts for
    collisions between reactants and bath gas species.
    """
    _reaction_rate_type = "linear-Burke"

    def __cinit__(self, input_data=None, init=True):
        self.set_cxx_object()

        if init:
            if isinstance(input_data, dict):
                self._rate = make_shared[CxxLinearBurkeRate](py_to_anymap(input_data))
            elif input_data:
                raise TypeError("'input_data' must be a dict")
            else:
                raise ValueError("No input data provided")
            self.set_cxx_object()

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxLinearBurkeRate):
        return cython.cast(cython.pointer(CxxLinearBurkeRate), self.rate)


@cython.cclass
class ChebyshevRate(ReactionRate):
    r"""
    A pressure-dependent reaction rate parameterized by a bivariate Chebyshev polynomial
    in temperature and pressure.
    """
    _reaction_rate_type = "Chebyshev"

    def __cinit__(self, temperature_range=None, pressure_range=None, data=None,
                  input_data=None, init=True):

        if init:
            if isinstance(input_data, dict):
                self._rate = make_shared[CxxChebyshevRate](py_to_anymap(input_data))
            elif all([arg is not None
                    for arg in [temperature_range, pressure_range, data]]):
                Tmin: cython.double = temperature_range[0]
                Tmax: cython.double = temperature_range[1]
                Pmin: cython.double = pressure_range[0]
                Pmax: cython.double = pressure_range[1]
                cxxdata: CxxArray2D = self._cxxarray2d(data)
                self._rate = make_shared[CxxChebyshevRate](
                    Tmin, Tmax, Pmin, Pmax, cxxdata
                    )
            elif all([arg is None
                    for arg in [temperature_range, pressure_range, data, input_data]]):
                self._rate = make_shared[CxxChebyshevRate](py_to_anymap({}))
            elif input_data:
                raise TypeError("Invalid parameter 'input_data'")
            else:
                raise TypeError("Invalid parameters")
            self.set_cxx_object()

    def __init__(
        self,
        temperature_range: _Sequence[float] | None = None,
        pressure_range: _Sequence[float] | None = None,
        data: _ArrayLike | None = None,
        input_data: _ChebyshevRateInput | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def __call__(self, temperature: cython.double, pressure: cython.double) -> float:
        """
        Evaluate rate expression based on temperature and pressure.
        """
        return self.rate.eval(temperature, pressure)

    @cython.cfunc
    def _cxxarray2d(self, coeffs) -> CxxArray2D:
        """ Internal function to assign coefficient matrix values """
        data: CxxArray2D
        if isinstance(coeffs, np.ndarray):
            coeffs = coeffs.tolist()
        data.resize(len(coeffs), len(coeffs[0]))
        value: cython.double
        i: cython.int
        j: cython.int
        for i, row in enumerate(coeffs):
            for j, value in enumerate(row):
                CxxArray2D_set(data, i, j, value)

        return data

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxChebyshevRate):
        return cython.cast(cython.pointer(CxxChebyshevRate), self.rate)

    @property
    def temperature_range(self) -> tuple[float, float]:
        """ Valid temperature range [K] for the Chebyshev fit """
        return self.cxx_object().Tmin(), self.cxx_object().Tmax()

    @property
    def pressure_range(self) -> tuple[float, float]:
        """ Valid pressure range [Pa] for the Chebyshev fit """
        return self.cxx_object().Pmin(), self.cxx_object().Pmax()

    @property
    def n_temperature(self) -> int:
        """
        Number of temperatures over which the Chebyshev fit is computed.
        (same as number of rows of `data` property).
        """
        return self.cxx_object().nTemperature()

    @property
    def n_pressure(self) -> int:
        """
        Number of pressures over which the Chebyshev fit is computed
        (same as number of columns of `data` property).
        """
        return self.cxx_object().nPressure()

    @property
    def data(self) -> _Array:
        """
        2D array of Chebyshev coefficients where rows and columns correspond to
        temperature and pressure dimensions over which the Chebyshev fit is computed.
        """
        # Must split declaration and assignment: 'data()' returns CxxArray2D&
        # and assigning it to a typed local copies the value into a CxxArray2D
        cxxcoeffs: CxxArray2D
        cxxcoeffs = self.cxx_object().data()
        c = np.fromiter(cxxcoeffs.data(), np.double)
        return c.reshape(cxxcoeffs.nRows(), cxxcoeffs.nColumns(), order="F")


@cython.cclass
class CustomRate(ReactionRate):
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
            self._rate = make_shared[CxxCustomFunc1Rate]()
            self.set_cxx_object()
            try:
                self.set_rate_function(k)
            except Exception:
                raise TypeError(
                    f"Cannot convert input with type '{type(k)}' to rate expression.")

    def __init__(self, k: _Func1Like | None = None, init: bool = True) -> None:
        """Published constructor signature."""

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxCustomFunc1Rate):
        return cython.cast(cython.pointer(CxxCustomFunc1Rate), self.rate)

    def set_rate_function(self, k: _Func1Like) -> None:
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


@cython.cclass
class ExtensibleRate(ReactionRate):
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

    delegatable_methods: dict[str, tuple[str, str, str]] = {
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

    def set_parameters(self, params: _AnyMap, rate_coeff_units: _UnitStack) -> None:
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

    def get_parameters(self, params: _AnyMap) -> None:
        """
        Responsible for serializing the state of the ExtensibleRate object, using the
        same format as a YAML reaction entry. This is the inverse of `set_parameters`.

        Serialization methods may request output in unit systems other than Cantera's
        native mks+kmol system. To enable conversions to the user-specified unit system,
        dimensional values should be added to ``params`` using the methods
        `AnyMap.set_quantity` and `AnyMap.set_activation_energy`.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.get_parameters")

    def eval(self, data: "ExtensibleRateData") -> float:
        """
        Responsible for calculating the forward rate constant based on the current state
        of the phase, stored in an instance of a class derived from
        `ExtensibleRateData`.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.eval")

    def validate(self, equation: str, soln: "_Solution") -> None:
        """
        Responsible for validating that the rate expression is configured with valid
        parameters. This may depend on properties of the Solution, for example
        temperature ranges over which the rate expression can be evaluated. Raises an
        exception if any validation fails.
        """
        pass

    @cython.cfunc
    def set_cxx_object(self, rate: cython.pointer(CxxReactionRate) = cython.NULL):
        drate: cython.pointer(CxxDelegator)
        handle: shared_ptr[CxxExternalHandle]

        if rate is cython.NULL:
            # Started with Python object first. Create the C++ object and attach to it.
            # In this case, the Python object owns the C++ object, via self._rate
            self._rate = make_shared[CxxReactionRateDelegator]()
            self.rate = self._rate.get()
            drate = dynamic_cast[CxxDelegatorPtr](self.rate)
            handle = make_shared[CxxPythonHandle](
                cython.cast(cython.pointer(PyObject), self), True)
            drate.holdExternalHandle(stringify('python'), handle)
        else:
            # Set up Python object from a C++ object that was created first. In this
            # case, the C++ object owns the Python object, and self._rate is empty to
            # avoid creating a circular dependency.
            self._rate.reset()
            self.rate = rate

        assign_delegates(self, dynamic_cast[CxxDelegatorPtr](self.rate))
        cython.cast(cython.pointer(CxxReactionRateDelegator), self.rate).setType(
            stringify(self._reaction_rate_type))


@cython.cclass
class ExtensibleRateData:
    """
    A base class for data used when evaluating instances of `ExtensibleRate`. Classes
    derived from `ExtensibleRateData` are used to store common data needed to evaluate
    all reactions of a particular type.

    Classes derived from `ExtensibleRateData` must implement the `update` method. After
    the `update` method has been called, instances of `ExtensibleRateData` are passed as
    the argument to `ExtensibleRate.eval`.

    .. versionadded:: 3.0
    """
    delegatable_methods: dict[str, tuple[str, str, str]] = {
        "update": ("update", "double(void*)", "replace")
    }

    def update(self, soln: "_Solution") -> bool:
        """
        This method takes a `Solution` object and stores any thermodynamic data (for
        example, temperature and pressure) needed to evaluate all reactions of the
        corresponding ReactionRate type.

        If this state data has changed since the last time `update` was called and the
        reaction rates need to be updated, this method should return `True`. Otherwise,
        it should return `False`.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.update")

    @cython.cfunc
    def set_cxx_object(self, data: cython.pointer(CxxReactionDataDelegator)):
        assign_delegates(self, dynamic_cast[CxxDelegatorPtr](data))


@cython.cclass
class InterfaceRateBase(ArrheniusRateBase):
    """
    Base class collecting commonly used features of Arrhenius-type rate objects
    that include coverage dependencies.
    """

    def __call__(self, temperature: cython.double, coverages: _Array) -> float:
        """
        Evaluate rate expression based on temperature and surface coverages.

        .. warning::

            This method is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        cxxdata: vector[cython.double]
        for c in coverages:
            cxxdata.push_back(c)
        return self.rate.eval(temperature, cxxdata)

    @property
    def coverage_dependencies(self) -> dict[str, _CoverageParameters]:
        """
        Get/set a dictionary containing adjustments to the Arrhenius rate expression
        dependent on surface species coverages. The keys of the dictionary are species
        names, and the values are dictionaries specifying the three coverage
        parameters ``a``, ``m`` and ``E`` which are the modifiers for the pre-exponential
        factor [m, kmol, s units], the temperature exponent [nondimensional],
        and the activation energy [J/kmol], respectively.
        """
        cxx_deps: CxxAnyMap
        self.interface.getCoverageDependencies(cxx_deps)
        return anymap_to_py(cxx_deps)

    @coverage_dependencies.setter
    def coverage_dependencies(self, deps: dict[str, _CoverageParameters]) -> None:
        cxx_deps: CxxAnyMap = py_to_anymap(deps)
        self.interface.setCoverageDependencies(cxx_deps)

    def set_species(self, species: _Iterable[str]) -> None:
        """
        Set association with an ordered list of all species associated with an
        `InterfaceKinetics` object.
        """
        cxxvector: vector[string]
        for s in species:
            cxxvector.push_back(stringify(s))
        self.interface.setSpecies(cxxvector)

    @property
    def site_density(self) -> float:
        """
        Site density [kmol/m²].

        The site density is not an independent property, as it is set by an associated
        `InterfaceKinetics` object. Accordingly, the setter should be only used for
        testing purposes, as the value will be overwritten by an update of the
        thermodynamic state.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return self.interface.siteDensity()

    @site_density.setter
    def site_density(self, site_density: float) -> None:
        self.interface.setSiteDensity(site_density)

    @property
    def uses_electrochemistry(self) -> bool:
        """
        Return boolean flag indicating whether rate involves a charge transfer.
        """
        return self.interface.usesElectrochemistry()

    @property
    def beta(self) -> float:
        """
        Return the charge transfer beta parameter
        """
        return self.interface.beta()


@cython.cclass
class InterfaceArrheniusRate(InterfaceRateBase):
    r"""
    A reaction rate coefficient which depends on temperature and interface coverage
    """
    _reaction_rate_type = "interface-Arrhenius"

    def __cinit__(self, A=None, b=None, Ea=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea=Ea)

    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea: float | None = None,
        input_data: _ReactionRateInput[_ArrheniusParameters] | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def _from_dict(self, input_data: _ReactionRateInput[_ArrheniusParameters]) -> None:
        self._rate = make_shared[CxxInterfaceArrheniusRate](py_to_anymap(input_data))

    def _from_parameters(self, A: cython.double, b: cython.double,
                         Ea: cython.double) -> None:
        self._rate = make_shared[CxxInterfaceArrheniusRate](A, b, Ea)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxArrheniusRate), self.rate)
        self.interface = cython.cast(cython.pointer(CxxInterfaceRateBase), self.cxx_object())

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxInterfaceArrheniusRate):
        return cython.cast(cython.pointer(CxxInterfaceArrheniusRate), self.rate)


@cython.cclass
class InterfaceBlowersMaselRate(InterfaceRateBase):
    r"""
    A reaction rate coefficient which depends on temperature and enthalpy change
    of the reaction follows the Blowers-Masel approximation and modified Arrhenius form
    described in `ArrheniusRate`.
    """
    _reaction_rate_type = "interface-Blowers-Masel"

    def __cinit__(self, A=None, b=None, Ea0=None, w=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea0=Ea0, w=w)

    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea0: float | None = None,
        w: float | None = None,
        input_data: _ReactionRateInput[_BlowersMaselParameters] | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def _from_dict(
        self, input_data: _ReactionRateInput[_BlowersMaselParameters]
    ) -> None:
        self._rate = make_shared[CxxInterfaceBlowersMaselRate](py_to_anymap(input_data))

    def _from_parameters(self, A: cython.double, b: cython.double,
                         Ea0: cython.double, w: cython.double) -> None:
        self._rate = make_shared[CxxInterfaceBlowersMaselRate](A, b, Ea0, w)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxArrheniusRate), self.rate)
        self.interface = cython.cast(cython.pointer(CxxInterfaceRateBase), self.cxx_object())

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxInterfaceBlowersMaselRate):
        return cython.cast(cython.pointer(CxxInterfaceBlowersMaselRate), self.rate)

    @property
    def bond_energy(self) -> float:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        return self.cxx_object().bondEnergy()

    @property
    def delta_enthalpy(self) -> float:
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
        return self.cxx_object().deltaH()

    @delta_enthalpy.setter
    def delta_enthalpy(self, delta_H: float) -> None:
        self.cxx_object().setDeltaH(delta_H)


@cython.cclass
class StickRateBase(InterfaceRateBase):
    """
    Base class collecting commonly used features of Arrhenius-type sticking rate
    objects that include coverage dependencies.
    """

    @property
    def motz_wise_correction(self) -> bool:
        """
        Get/Set a boolean indicating whether to use the correction factor developed by
        Motz & Wise for reactions with high (near-unity) sticking coefficients when
        converting the sticking coefficient to a rate coefficient.
        """
        return self.stick.motzWiseCorrection()

    @motz_wise_correction.setter
    def motz_wise_correction(self, motz_wise: bool) -> None:
        self.stick.setMotzWiseCorrection(motz_wise)

    @property
    def sticking_species(self) -> str:
        """
        The name of the sticking species. Needed only for reactions with
        multiple non-surface reactant species, where the sticking species is
        ambiguous.
        """
        return pystr(self.stick.stickingSpecies())

    @sticking_species.setter
    def sticking_species(self, species: str) -> None:
        self.stick.setStickingSpecies(stringify(species))

    @property
    def sticking_order(self) -> float:
        """
        The exponent applied to site density (sticking order).

        The sticking order is not an independent property and is detected automatically
        by Cantera. Accordingly, the setter should be only used for testing purposes.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return self.stick.stickingOrder()

    @sticking_order.setter
    def sticking_order(self, order: float) -> None:
        self.stick.setStickingOrder(order)

    @property
    def sticking_weight(self) -> float:
        """
        The molecular weight of the sticking species.

        The sticking weight is not an independent property and is detected automatically
        by Cantera. Accordingly, the setter should be only used for testing purposes.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return self.stick.stickingWeight()

    @sticking_weight.setter
    def sticking_weight(self, weight: float) -> None:
        self.stick.setStickingWeight(weight)


@cython.cclass
class StickingArrheniusRate(StickRateBase):
    r"""
    A surface sticking rate expression based on the Arrhenius parameterization
    """
    _reaction_rate_type = "sticking-Arrhenius"

    def __cinit__(self, A=None, b=None, Ea=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea=Ea)

    def _from_dict(self, input_data: _ReactionRateInput[_ArrheniusParameters]) -> None:
        self._rate = make_shared[CxxStickingArrheniusRate](py_to_anymap(input_data))

    def _from_parameters(self, A: cython.double, b: cython.double,
                         Ea: cython.double) -> None:
        self._rate = make_shared[CxxStickingArrheniusRate](A, b, Ea)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxArrheniusRate), self.rate)
        self.stick = cython.cast(cython.pointer(CxxStickingCoverage), self.cxx_object())
        self.interface = cython.cast(cython.pointer(CxxInterfaceRateBase), self.stick)

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxStickingArrheniusRate):
        return cython.cast(cython.pointer(CxxStickingArrheniusRate), self.rate)


@cython.cclass
class StickingBlowersMaselRate(StickRateBase):
    r"""
    A surface sticking rate expression based on the Blowers-Masel parameterization
    """
    _reaction_rate_type = "sticking-Blowers-Masel"

    def __cinit__(self, A=None, b=None, Ea0=None, w=None, input_data=None, init=True):

        if init:
            self._cinit(input_data, A=A, b=b, Ea0=Ea0, w=w)

    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea0: float | None = None,
        w: float | None = None,
        input_data: _ReactionRateInput[_BlowersMaselParameters] | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    def _from_dict(
        self, input_data: _ReactionRateInput[_BlowersMaselParameters]
    ) -> None:
        self._rate = make_shared[CxxStickingBlowersMaselRate](py_to_anymap(input_data))

    def _from_parameters(self, A: cython.double, b: cython.double,
                         Ea0: cython.double, w: cython.double) -> None:
        self._rate = make_shared[CxxStickingBlowersMaselRate](A, b, Ea0, w)

    @cython.cfunc
    def set_cxx_object(self):
        self.rate = self._rate.get()
        self.base = cython.cast(cython.pointer(CxxArrheniusRate), self.rate)
        self.stick = cython.cast(cython.pointer(CxxStickingCoverage), self.cxx_object())
        self.interface = cython.cast(cython.pointer(CxxInterfaceRateBase), self.stick)

    @cython.cfunc
    def cxx_object(self) -> cython.pointer(CxxStickingBlowersMaselRate):
        return cython.cast(cython.pointer(CxxStickingBlowersMaselRate), self.rate)

    @property
    def bond_energy(self) -> float:
        """
        Average bond dissociation energy of the bond being formed and broken
        in the reaction ``E`` [J/kmol].
        """
        return self.cxx_object().bondEnergy()

    @property
    def delta_enthalpy(self) -> float:
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
        return self.cxx_object().deltaH()

    @delta_enthalpy.setter
    def delta_enthalpy(self, delta_H: float) -> None:
        self.cxx_object().setDeltaH(delta_H)


@cython.cclass
class ThirdBody:
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
        self._third_body = make_shared[CxxThirdBody](stringify(collider))
        self.third_body = self._third_body.get()

        if efficiencies is not None:
            self.efficiencies = efficiencies

        if default_efficiency is not None:
            self.default_efficiency = default_efficiency

    def __init__(
        self,
        collider: _Str = "M",
        *,
        efficiencies: dict[str, float] | None = None,
        default_efficiency: float | None = None,
        init: bool = True,
    ) -> None:
        """Published constructor signature."""

    @cython.cfunc
    @staticmethod
    def wrap(third_body: shared_ptr[CxxThirdBody]):
        tb = ThirdBody(init=False)
        tb._third_body = third_body
        tb.third_body = tb._third_body.get()
        return tb

    @property
    def name(self) -> str:
        """
        Get the name of the third-body collider used in the reaction equation.
        """
        return pystr(self.third_body.name())

    @property
    def mass_action(self) -> bool:
        """
        Retrieve flag indicating whether third-body collider participates
        in the law of mass action.
        """
        return self.third_body.mass_action

    @property
    def efficiencies(self) -> dict[str, float]:
        """
        Get/Set a `dict` defining non-default third-body efficiencies for this reaction,
        where the keys are the species names and the values are the efficiencies.
        """
        return comp_map_to_dict(self.third_body.efficiencies)

    @efficiencies.setter
    def efficiencies(self, eff: dict[str, float]) -> None:
        self.third_body.efficiencies = comp_map(eff)

    @property
    def default_efficiency(self) -> float:
        """
        Get/Set the default third-body efficiency for this reaction, used for species
        not in `efficiencies`.
        """
        return self.third_body.default_efficiency

    @default_efficiency.setter
    def default_efficiency(self, default_eff: float) -> None:
        self.third_body.default_efficiency = default_eff

    def efficiency(self, species: _Str) -> float:
        """
        Get the efficiency of the third body named ``species`` considering both
        the default efficiency and species-specific efficiencies.
        """
        return self.third_body.efficiency(stringify(species))


@cython.cclass
class Reaction:
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

        _rate: ReactionRate
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
        elif callable(rate):
            _rate = CustomRate(rate)
        else:
            raise TypeError(f"Invalid rate definition with type '{type(rate)}'")

        _third_body: ThirdBody
        if isinstance(third_body, ThirdBody):
            _third_body = third_body
        elif isinstance(third_body, str):
            _third_body = ThirdBody(third_body)

        if reactants and products:
            # create from reactant and product compositions
            if third_body:
                self._reaction = make_shared[CxxReaction](
                    comp_map(reactants), comp_map(products),
                    _rate._rate, _third_body._third_body
                )
            else:
                self._reaction = make_shared[CxxReaction](
                    comp_map(reactants), comp_map(products),
                    _rate._rate
                )
        elif equation:
            # create from reaction equation
            if third_body:
                self._reaction = make_shared[CxxReaction](
                    stringify(equation),
                    _rate._rate, _third_body._third_body
                )
            else:
                self._reaction = make_shared[CxxReaction](
                    stringify(equation),
                    _rate._rate
                )
        else:
            # create default object
            raise ValueError("Missing reactant and/or product information.")

        self.reaction = self._reaction.get()
        self._rate = _rate

    def __init__(
        self,
        reactants: dict[str, float] | _Str | None = None,
        products: dict[str, float] | _Str | None = None,
        rate: ReactionRate
        | _ReactionRateInput[_ReactionRateParameters]
        | _ArrheniusParameters
        | _Callable[[float], float]
        | None = None,
        *,
        equation: _Str | None = None,
        init: bool = True,
        third_body: ThirdBody | _Str | None = None,
    ) -> None:
        """Published constructor signature."""

    @cython.cfunc
    @staticmethod
    def wrap(reaction: shared_ptr[CxxReaction]):
        """
        Wrap a C++ Reaction object with a Python object of the correct derived type.
        """
        # wrap C++ reaction
        R: Reaction
        R = Reaction(init=False)
        R._reaction = reaction
        R.reaction = R._reaction.get()
        R._rate = ReactionRate.wrap(R.reaction.rate())
        return R

    @classmethod
    def from_dict(
        cls,
        data: _ReactionRateInput[_ReactionRateParameters],
        kinetics: "_Kinetics",
        hyphenize: bool = True,
    ) -> "Reaction":
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
        any_map: CxxAnyMap = py_to_anymap(data, hyphenize=hyphenize)
        cxx_reaction = CxxNewReaction(any_map, cython.cast(Kinetics, kinetics).kinetics[0])
        return Reaction.wrap(cxx_reaction)

    @classmethod
    def from_yaml(cls, text: _Str, kinetics: "_Kinetics") -> "Reaction":
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
        any_map: CxxAnyMap = AnyMapFromYamlString(stringify(text))
        cxx_reaction = CxxNewReaction(any_map, cython.cast(Kinetics, kinetics).kinetics[0])
        return Reaction.wrap(cxx_reaction)

    @staticmethod
    def list_from_file(filename: _Str, kinetics: "_Kinetics",
                       section: _Str = "reactions") -> list["Reaction"]:
        """
        Create a list of Reaction objects from all of the reactions defined in a
        YAML file. Reactions from the section ``section`` will be returned.

        Directories on Cantera's input file path will be searched for the
        specified file.
        """
        root = AnyMapFromYamlFile(stringify(str(filename)))
        cxx_reactions = CxxGetReactions(root[stringify(section)],
                                        cython.cast(Kinetics, kinetics).kinetics[0])
        return [Reaction.wrap(r) for r in cxx_reactions]

    @staticmethod
    def list_from_yaml(text: _Str, kinetics: "_Kinetics") -> list["Reaction"]:
        """
        Create a list of `Reaction` objects from all the reactions defined in a
        YAML string.
        """
        root = AnyMapFromYamlString(stringify(text))
        cxx_reactions = CxxGetReactions(root[stringify("items")],
                                        cython.cast(Kinetics, kinetics).kinetics[0])
        return [Reaction.wrap(r) for r in cxx_reactions]

    @property
    def reactant_string(self) -> str:
        """
        A string representing the reactants side of the chemical equation for
        this reaction. Determined automatically based on `reactants`.
        """
        return pystr(self.reaction.reactantString())

    @property
    def product_string(self) -> str:
        """
        A string representing the products side of the chemical equation for
        this reaction. Determined automatically based on `products`.
        """
        return pystr(self.reaction.productString())

    @property
    def equation(self) -> str:
        """
        A string giving the chemical equation for this reaction. Determined
        automatically based on `reactants` and `products`.
        """
        return pystr(self.reaction.equation())

    @property
    def reactants(self) -> dict[str, float]:
        """
        Get/Set the reactants in this reaction as a dict where the keys are
        species names and the values, are the stoichiometric coefficients, for example
        ``{'CH4':1, 'OH':1}``, or as a composition string, for example
        ``'CH4:1, OH:1'``.

        .. versionchanged:: 3.1  This is a read-only property
        """
        return comp_map_to_dict(self.reaction.reactants)

    @property
    def products(self) -> dict[str, float]:
        """
        Get/Set the products in this reaction as a dict where the keys are
        species names and the values, are the stoichiometric coefficients, for example
        ``{'CH3':1, 'H2O':1}``, or as a composition string, for example
        ``'CH3:1, H2O:1'``.

        .. versionchanged:: 3.1  This is a read-only property
        """
        return comp_map_to_dict(self.reaction.products)

    def __contains__(self, species: _Str, /) -> bool:
        return species in self.reactants or species in self.products

    @property
    def orders(self) -> dict[str, float]:
        """
        Get/Set the reaction order with respect to specific species as a dict
        with species names as the keys and orders as the values, or as a
        composition string. By default, mass-action kinetics is assumed, with
        the reaction order for each reactant species equal to each its
        stoichiometric coefficient.
        """
        return comp_map_to_dict(self.reaction.orders)

    @orders.setter
    def orders(self, orders: dict[str, float]) -> None:
        self.reaction.orders = comp_map(orders)

    @property
    def ID(self) -> str:
        """
        Get/Set the identification string for the reaction, which can be used in
        filtering operations.
        """
        return pystr(self.reaction.id)

    @ID.setter
    def ID(self, ID: _Str) -> None:
        self.reaction.id = stringify(ID)

    @property
    def reaction_type(self) -> str:
        """
        Retrieve the native type name of the reaction.
        """
        return pystr(self.reaction.type())

    @property
    def rate(self) -> ReactionRate:
        """ Get/Set the reaction rate evaluator for this reaction. """
        return self._rate

    @rate.setter
    def rate(self, rate: ReactionRate | _Callable[[float], float]) -> None:
        if isinstance(rate, ReactionRate):
            self._rate = rate
        elif callable(rate):
            self._rate = CustomRate(rate)
        else:
            raise TypeError(f"Invalid rate definition with type '{type(rate)}'")
        self.reaction.setRate(self._rate._rate)

    @property
    def reversible(self) -> bool:
        """
        Get/Set a flag which is `True` if this reaction is reversible or `False`
        otherwise.
        """
        return self.reaction.reversible

    @reversible.setter
    def reversible(self, reversible: bool) -> None:
        self.reaction.reversible = reversible

    @property
    def duplicate(self) -> bool:
        """
        Get/Set a flag which is `True` if this reaction is marked as a duplicate
        or `False` otherwise.
        """
        return self.reaction.duplicate

    @duplicate.setter
    def duplicate(self, duplicate: bool) -> None:
        self.reaction.duplicate = duplicate

    @property
    def allow_nonreactant_orders(self) -> bool:
        """
        Get/Set a flag which is `True` if reaction orders can be specified for
        non-reactant species. Default is `False`.
        """
        return self.reaction.allow_nonreactant_orders

    @allow_nonreactant_orders.setter
    def allow_nonreactant_orders(self, allow: bool) -> None:
        self.reaction.allow_nonreactant_orders = allow

    @property
    def allow_negative_orders(self) -> bool:
        """
        Get/Set a flag which is `True` if negative reaction orders are allowed.
        Default is `False`.
        """
        return self.reaction.allow_negative_orders

    @allow_negative_orders.setter
    def allow_negative_orders(self, allow: bool) -> None:
        self.reaction.allow_negative_orders = allow

    @property
    def input_data(self) -> _ReactionRateInput[_ReactionRateParameters]:
        """
        Get input data for this reaction with its current parameter values,
        along with any user-specified data provided with its input (YAML)
        definition.
        """
        return anymap_to_py(self.reaction.parameters(True))

    def update_user_data(self, data: dict[str, str | float]) -> None:
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.reaction.input.update(py_to_anymap(data), False)

    def clear_user_data(self) -> None:
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

    @property
    def rate_coeff_units(self) -> _Units:
        """Get reaction rate coefficient units"""
        rate_units: CxxUnits = self.reaction.rate_units
        return Units.copy(rate_units)

    @property
    def third_body(self) -> ThirdBody | None:
        """
        Returns a `ThirdBody` object if `Reaction` uses a third body collider, and
        ``None`` otherwise.

        .. versionadded:: 3.0
        """
        if self.reaction.usesThirdBody():
            return ThirdBody.wrap(self.reaction.thirdBody())

    @property
    def third_body_name(self) -> str | None:
        """
        Returns name of `ThirdBody` collider if `Reaction` uses a third body collider,
        and ``None`` otherwise.

        .. versionadded:: 3.0
        """
        if self.reaction.usesThirdBody():
            return ThirdBody.wrap(self.reaction.thirdBody()).name
        return None
