from typing import Sequence, TypedDict, overload

from cantera._types import Array

UnitDict = TypedDict(
    "UnitDict",
    {
        "activation-energy": str,
        "current": str,
        "energy": str,
        "length": str,
        "mass": str,
        "pressure": str,
        "quantity": str,
        "temperature": str,
        "time": str,
    },
    total=False,
)
UnitDictBytes = TypedDict(
    "UnitDictBytes",
    {
        "activation-energy": bytes,
        "current": bytes,
        "energy": bytes,
        "length": bytes,
        "mass": bytes,
        "pressure": bytes,
        "quantity": bytes,
        "temperature": bytes,
        "time": bytes,
    },
    total=False,
)

class Units:
    def dimension(self, primary: str) -> float: ...
    @property
    def dimensions(self) -> dict[str, float]: ...
    @property
    def factor(self) -> float: ...

class UnitStack:
    def product(self) -> Units: ...
    def join(self, exponent: float) -> None: ...

class UnitSystem:
    def defaults(self) -> UnitDictBytes: ...
    @property
    def units(self) -> UnitDict: ...
    @units.setter
    def units(self, units: UnitDict | UnitDictBytes) -> None: ...
    @overload
    def convert_to(self, quantity: str | float, dest: str | Units) -> float: ...
    @overload
    def convert_to(self, quantity: Array, dest: str | Units) -> Array: ...
    @overload
    def convert_to(
        self, quantity: Sequence[str | float], dest: str | Units
    ) -> list[float]: ...
    @overload
    def convert_activation_energy_to(
        self, quantity: str | float, dest: str | Units
    ) -> float: ...
    @overload
    def convert_activation_energy_to(
        self, quantity: Array, dest: str | Units
    ) -> Array: ...
    @overload
    def convert_activation_energy_to(
        self, quantity: Sequence[str | float], dest: str | Units
    ) -> list[float]: ...
    @overload
    def convert_rate_coeff_to(
        self, quantity: str | float, dest: str | Units | UnitStack
    ) -> float: ...
    @overload
    def convert_rate_coeff_to(
        self, quantity: Array, dest: str | Units | UnitStack
    ) -> Array: ...
    @overload
    def convert_rate_coeff_to(
        self, quantity: Sequence[str | float], dest: str | Units | UnitStack
    ) -> list[float]: ...
