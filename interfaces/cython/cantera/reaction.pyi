# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections.abc import Callable, Iterable, Sequence
from typing import (
    Generic,
    TypeAlias,
    TypedDict,
    TypeVar,
)

from typing_extensions import NotRequired

from ._types import Array, ArrayLike
from .composite import Solution
from .kinetics import Kinetics
from .units import Units

class ArrheniusParameters(TypedDict):
    A: float
    b: float
    Ea: float

class BlowersMaselParameters(TypedDict):
    A: float
    b: float
    Ea0: float
    w: float

class TwoTempPlasmaParameters(TypedDict):
    A: float
    b: float
    Ea_gas: float
    Ea_electron: float

class ElectronCollisionPlasmaParameters(TypedDict):
    energy_levels: list[float]
    cross_sections: list[float]

class PlogParameters(ArrheniusParameters):
    P: float

ReactionRateParameters: TypeAlias = (
    ArrheniusParameters
    | BlowersMaselParameters
    | TwoTempPlasmaParameters
    | ElectronCollisionPlasmaParameters
    | PlogParameters
)

class TroeParameters(TypedDict):
    A: float
    T3: float
    T1: float
    T2: NotRequired[float]

class CoverageParameters(TypedDict):
    a: float
    m: float
    E: float

T = TypeVar("T")

class ReactionRateInput(TypedDict, Generic[T], total=False):
    type: str
    rate_constant: T
    efficiencies: dict[str, float]
    Troe: TroeParameters
    coverage_dependencies: dict[str, CoverageParameters]

class ReactionInput(ReactionRateInput[T]):
    equation: str

class FalloffRateInput(TypedDict, total=False):
    type: str
    low_P_rate_constant: ArrheniusParameters
    high_P_rate_constant: ArrheniusParameters
    Troe: TroeParameters

class PlogRateInput(TypedDict, total=False):
    type: str
    rate_constants: Sequence[PlogParameters]

class ReactionRate:
    def __call__(self, temperature: float) -> float: ...
    @property
    def type(self) -> str: ...
    @property
    def sub_type(self) -> str: ...
    @classmethod
    def from_dict(
        cls, data: ReactionRateInput[ReactionRateParameters], hyphenize: bool = True
    ) -> ReactionRate: ...
    @classmethod
    def from_yaml(cls, data: str, hyphenize: bool = True) -> ReactionRate: ...
    @property
    def input_data(self) -> ReactionRateInput[ReactionRateParameters]: ...
    @property
    def conversion_units(self) -> Units: ...

class ArrheniusRateBase(ReactionRate):
    @property
    def pre_exponential_factor(self) -> float: ...
    @property
    def temperature_exponent(self) -> float: ...
    @property
    def activation_energy(self) -> float: ...
    @property
    def allow_negative_pre_exponential_factor(self) -> bool: ...
    @allow_negative_pre_exponential_factor.setter
    def allow_negative_pre_exponential_factor(self, allow: bool) -> None: ...

class ArrheniusRate(ArrheniusRateBase):
    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea: float | None = None,
        input_data: ReactionRateInput[ArrheniusParameters] | None = None,
        init: bool = True,
    ) -> None: ...
    def _from_parameters(self, A: float, b: float, Ea: float) -> None: ...

class BlowersMaselRate(ArrheniusRateBase):
    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea0: float | None = None,
        w: float | None = None,
        input_data: ReactionRateInput[BlowersMaselParameters] | None = None,
        init: bool = True,
    ) -> None: ...
    def _from_parameters(self, A: float, b: float, Ea0: float, w: float) -> None: ...
    @property
    def bond_energy(self) -> float: ...
    @property
    def delta_enthalpy(self) -> float: ...
    @delta_enthalpy.setter
    def delta_enthalpy(self, delta_H: float) -> None: ...

class TwoTempPlasmaRate(ArrheniusRateBase):
    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea_gas: float = 0.0,
        Ea_electron: float = 0.0,
        input_data: ReactionRateInput[TwoTempPlasmaParameters] | None = None,
        init: bool = True,
    ) -> None: ...
    def _from_parameters(
        self, A: float, b: float, Ea_gas: float, Ea_electron: float
    ) -> None: ...
    @property
    def activation_electron_energy(self) -> float: ...

class ElectronCollisionPlasmaRate(ReactionRate):
    def __init__(
        self,
        energy_levels: ArrayLike | None = None,
        cross_sections: ArrayLike | None = None,
        input_data: ReactionRateInput[ElectronCollisionPlasmaParameters] | None = None,
        init: bool = True,
    ) -> None: ...
    def _from_parameters(
        self, energy_levels: ArrayLike, cross_sections: ArrayLike
    ) -> None: ...
    @property
    def energy_levels(self) -> Array: ...
    @property
    def cross_sections(self) -> Array: ...

class FalloffRate(ReactionRate):
    def __init__(
        self,
        low: ArrheniusParameters | None = None,
        high: ArrheniusParameters | None = None,
        falloff_coeffs: Sequence[float] | None = None,
        input_data: ReactionRateInput[FalloffRateInput] | None = None,
        init: bool = True,
    ) -> None: ...
    def __call__(self, temperature: float, concm: float) -> float: ...  # type: ignore[override]
    @property
    def low_rate(self) -> Arrhenius: ...
    @low_rate.setter
    def low_rate(self, rate: Arrhenius) -> None: ...
    @property
    def high_rate(self) -> Arrhenius: ...
    @high_rate.setter
    def high_rate(self, rate: Arrhenius) -> None: ...
    @property
    def falloff_coeffs(self) -> Array: ...
    @falloff_coeffs.setter
    def falloff_coeffs(self, data: Iterable[float]) -> None: ...
    @property
    def allow_negative_pre_exponential_factor(self) -> bool: ...
    @allow_negative_pre_exponential_factor.setter
    def allow_negative_pre_exponential_factor(self, allow: bool) -> None: ...
    @property
    def chemically_activated(self) -> bool: ...
    @chemically_activated.setter
    def chemically_activated(self, activated: bool) -> None: ...
    def falloff_function(self, temperature: float, conc3b: float) -> float: ...

class LindemannRate(FalloffRate): ...
class TroeRate(FalloffRate): ...
class SriRate(FalloffRate): ...
class TsangRate(FalloffRate): ...

class PlogRate(ReactionRate):
    def __init__(
        self,
        rates: list[tuple[float, ArrheniusRate]] | None = None,
        input_data: PlogRateInput | None = None,
        init: bool = True,
    ) -> None: ...
    def __call__(self, temperature: float, pressure: float) -> float: ...  # type: ignore[override]
    @property
    def rates(self) -> list[tuple[float, ArrheniusRate]]: ...
    @rates.setter
    def rates(self, data: Iterable[tuple[float, ArrheniusRate]]) -> None: ...

class LinearBurkeRate(ReactionRate): ...

class ChebyshevRate(ReactionRate):
    @property
    def temperature_range(self) -> tuple[float, float]: ...
    @property
    def pressure_range(self) -> tuple[float, float]: ...
    @property
    def n_temperature(self) -> int: ...
    @property
    def n_pressure(self) -> int: ...
    @property
    def data(self) -> Array: ...

class CustomRate(ReactionRate): ...
class ExtensibleRate(ReactionRate): ...

class ExtensibleRateData:
    def update(self, soln: Solution) -> bool: ...

class InterfaceRateBase(ArrheniusRateBase):
    def __call__(self, temperature: float, coverages: Array) -> float: ...  # type: ignore[override]
    @property
    def coverage_dependencies(self) -> dict[str, CoverageParameters]: ...
    @coverage_dependencies.setter
    def coverage_dependencies(self, deps: dict[str, CoverageParameters]) -> None: ...
    def set_species(self, species: Iterable[str]) -> None: ...
    @property
    def site_density(self) -> float: ...
    @site_density.setter
    def site_density(self, site_density: float) -> None: ...
    @property
    def uses_electrochemistry(self) -> bool: ...
    @property
    def beta(self) -> float: ...

class InterfaceArrheniusRate(InterfaceRateBase):
    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea: float | None = None,
        input_data: ReactionRateInput[ArrheniusParameters] | None = None,
        init: bool = True,
    ) -> None: ...
    def _from_parameters(self, A: float, b: float, Ea: float) -> None: ...

class InterfaceBlowersMaselRate(InterfaceRateBase):
    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea0: float | None = None,
        w: float | None = None,
        input_data: ReactionRateInput[BlowersMaselParameters] | None = None,
        init: bool = True,
    ) -> None: ...
    def _from_dict(
        self, input_data: ReactionRateInput[BlowersMaselParameters]
    ) -> None: ...
    def _from_parameters(self, A: float, b: float, Ea0: float, w: float) -> None: ...
    @property
    def bond_energy(self) -> float: ...
    @property
    def delta_enthalpy(self) -> float: ...
    @delta_enthalpy.setter
    def delta_enthalpy(self, delta_H: float) -> None: ...

class StickRateBase(InterfaceRateBase): ...
class StickingArrheniusRate(StickRateBase): ...

class StickingBlowersMaselRate(StickRateBase):
    def __init__(
        self,
        A: float | None = None,
        b: float | None = None,
        Ea0: float | None = None,
        w: float | None = None,
        input_data: ReactionRateInput[BlowersMaselParameters] | None = None,
        init: bool = True,
    ) -> None: ...

class ThirdBody:
    def __init__(
        self,
        collider: str = "M",
        *,
        efficiencies: dict[str, float] | None = None,
        default_efficiency: float | None = None,
        init: bool = True,
    ) -> None: ...
    @property
    def name(self) -> str: ...
    @property
    def mass_action(self) -> bool: ...
    @property
    def efficiencies(self) -> dict[str, float]: ...
    @efficiencies.setter
    def efficiencies(self, eff: dict[str, float]) -> None: ...
    @property
    def default_efficiency(self) -> float: ...
    @default_efficiency.setter
    def default_efficiency(self, default_eff: float) -> None: ...
    def efficiency(self, species: str) -> float: ...

class Reaction:
    def __init__(
        self,
        reactants: dict[str, float] | None = None,
        products: dict[str, float] | None = None,
        rate: ReactionRate
        | ReactionRateInput[ReactionRateParameters]
        | ArrheniusParameters
        | Callable[[float], float]
        | None = None,
        *,
        equation: str | None = None,
        init: bool = True,
        third_body: ThirdBody | str | None = None,
    ) -> None: ...
    @classmethod
    def from_dict(
        cls,
        data: ReactionRateInput[ReactionRateParameters],
        kinetics: Kinetics,
        hyphenize: bool = True,
    ) -> Reaction: ...
    @classmethod
    def from_yaml(cls, data: str, kinetics: Kinetics) -> Reaction: ...
    @staticmethod
    def list_from_file(
        filename: str, kinetics: Kinetics, section: str = "reactions"
    ) -> list[Reaction]: ...
    @staticmethod
    def list_from_yaml(text: str, kinetics: Kinetics) -> list[Reaction]: ...
    @property
    def reactant_string(self) -> str: ...
    @property
    def product_string(self) -> str: ...
    @property
    def equation(self) -> str: ...
    @property
    def reactants(self) -> dict[str, float]: ...
    @property
    def products(self) -> dict[str, float]: ...
    def __contains__(self, species: str) -> bool: ...
    @property
    def orders(self) -> dict[str, float]: ...
    @orders.setter
    def orders(self, orders: dict[str, float]) -> None: ...
    @property
    def ID(self) -> str: ...
    @ID.setter
    def ID(self, ID: str) -> None: ...
    @property
    def reaction_type(self) -> str: ...
    @property
    def rate(self) -> ReactionRate: ...
    @rate.setter
    def rate(self, rate: ReactionRate | Callable[[float], float]) -> None: ...
    @property
    def reversible(self) -> bool: ...
    @reversible.setter
    def reversible(self, reversible: bool) -> None: ...
    @property
    def duplicate(self) -> bool: ...
    @duplicate.setter
    def duplicate(self, duplicate: bool) -> None: ...
    @property
    def allow_nonreactant_orders(self) -> bool: ...
    @allow_nonreactant_orders.setter
    def allow_nonreactant_orders(self, allow: bool) -> None: ...
    @property
    def allow_negative_orders(self) -> bool: ...
    @allow_negative_orders.setter
    def allow_negative_orders(self, allow: bool) -> None: ...
    @property
    def input_data(self) -> ReactionRateInput[ReactionRateParameters]: ...
    def update_user_data(self, data: dict[str, str | float]) -> None: ...
    def clear_user_data(self) -> None: ...
    @property
    def third_body(self) -> ThirdBody | None: ...
    @property
    def third_body_name(self) -> str | None: ...

class Arrhenius:
    def __init__(
        self, A: float = 0, b: float = 0, E: float = 0, init: bool = True
    ) -> None: ...
    @property
    def pre_exponential_factor(self) -> float: ...
    @property
    def temperature_exponent(self) -> float: ...
    @property
    def activation_energy(self) -> float: ...
    def __call__(self, T: float) -> float: ...
