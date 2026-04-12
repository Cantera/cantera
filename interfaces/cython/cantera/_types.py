# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections.abc import Callable, Mapping
from types import EllipsisType
from typing import (
    Any,
    Concatenate,
    Literal,
    ParamSpec,
    TypeAlias,
    TypedDict,
    TypeGuard,
    TypeVar,
    get_args,
    overload,
)

import numpy as np
from numpy.typing import ArrayLike as ArrayLike
from numpy.typing import NDArray
try:
    # Requires typing_extensions >= 4.13, or possibly Python >= 3.15
    from typing_extensions import TypeForm as TypeForm
except ImportError:
    # Wrong, but better than crashing with an ImportError at runtime
    from typing_extensions import Type as TypeForm


Array: TypeAlias = NDArray[np.float64]
Index: TypeAlias = EllipsisType | int | slice | tuple[EllipsisType | int | slice, ...]

Basis: TypeAlias = Literal["mole", "molar", "mass"]
EquilibriumSolver: TypeAlias = Literal["element_potential", "gibbs", "vcs", "auto"]
LogLevel: TypeAlias = int  # Literal[0, 1, 2, 3, 4, 5]
CompressionLevel: TypeAlias = Literal[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

# ThermoPhase state definitions
StateVariable: TypeAlias = Literal["T", "D", "P", "U", "H", "S", "V"]
PropertyPair: TypeAlias = Literal["TP", "TV", "HP", "SP", "SV", "UV"]
CompositionVariable: TypeAlias = Literal["X", "Y"]
FullState: TypeAlias = Literal[
    "TDX",
    "TDY",
    "TPX",
    "TPY",
    "UVX",
    "UVY",
    "DPX",
    "DPY",
    "HPX",
    "HPY",
    "SPX",
    "SPY",
    "SVX",
    "SVY",
]
CompositionLike: TypeAlias = str | Mapping[str, float | int] | ArrayLike
State2Setter: TypeAlias = tuple[float | None, float | None]
StateSetter: TypeAlias = tuple[float | None, float | None, CompositionLike]
ArrayCompositionLike: TypeAlias = str | dict[str, ArrayLike] | ArrayLike
ArrayState2Setter: TypeAlias = tuple[ArrayLike | None, ArrayLike | None]
ArrayStateSetter: TypeAlias = tuple[ArrayLike | None, ArrayLike | None, ArrayCompositionLike]

# PureFluid state definitions
PureFluidStateVariable: TypeAlias = Literal[StateVariable, "Q"]
PureFluidPropertyPair: TypeAlias = Literal[
    PropertyPair, "TQ", "PQ", "ST", "TV", "PV", "UP", "VH", "TH", "SH"
]
PureFluidFullState: TypeAlias = Literal[
    FullState, "TDQ", "TPQ", "UVQ", "DPQ", "HPQ", "SPQ", "SVQ"
]
PureFluidStateSetter: TypeAlias = tuple[float, float, float]
ArrayPureFluidStateSetter: TypeAlias = tuple[ArrayLike, ArrayLike, ArrayLike]

RefineCriteria = TypedDict(
    "RefineCriteria",
    {
        "ratio": float,
        "slope": float,
        "curve": float,
        "prune": float,
        "grid-min": float,
        "max-points": int,
    },
    total=False,
)


class StateDefinition(TypedDict, total=False):
    D: float
    H: float
    P: float
    S: float
    T: float
    V: float
    X: CompositionLike
    Y: CompositionLike

    DPX: StateSetter
    DPY: StateSetter
    HPX: StateSetter
    HPY: StateSetter
    SPX: StateSetter
    SPY: StateSetter
    SVX: StateSetter
    SVY: StateSetter
    TDX: StateSetter
    TDY: StateSetter
    TPX: StateSetter
    TPY: StateSetter
    UVX: StateSetter
    UVY: StateSetter


_T0 = TypeVar("_T0")


def literal_type_guard(tag: str, literal: TypeForm[_T0]) -> TypeGuard[_T0]:
    """Utility function for narrowing strings to specified literals.

    Typically used to check a string against the permissible keys of a
    TypedDict before accessing the corresponding value. For example,
    running pyright against the following code gives:

    .. code-block:: python

        class Inputs(TypedDict):
            A: float
            B: int

        inputs: Inputs = {"A": 1.0, "B": 2}
        option: str = "test"
        reveal_type(option)  # Type of "option" is "Literal['test']"

        if option in ["A", "B"]:
            reveal_type(option)  # Type of "option" is "Literal['test']"
            value = inputs[option]  # error: Could not access item in TypedDict
            reveal_type(value)  # Type of "value" is "Unknown"

        if literal_type_guard(option, Literal["A", "B"]):
            reveal_type(option)  # Type of "option" is "Literal['A', 'B']"
            value = inputs[option]
            reveal_type(value)  # Type of "value" is "float | int"

    :param tag: Input string to be checked and narrowed
    :param literal: String literal to check the tag against

    .. versionadded:: 3.2
    """
    return tag in get_args(literal)
