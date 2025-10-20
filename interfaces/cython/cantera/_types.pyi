# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections.abc import Callable
from types import EllipsisType
from typing import (
    Concatenate,
    Literal,
    TypeAlias,
    TypedDict,
    TypeGuard,
    TypeVar,
    overload,
)

import numpy as np
from numpy.typing import ArrayLike as ArrayLike
from numpy.typing import NDArray
from typing_extensions import ParamSpec, TypeForm

Array: TypeAlias = NDArray[np.float64]
Index: TypeAlias = EllipsisType | int | slice | tuple[EllipsisType | int | slice, ...]

Basis: TypeAlias = Literal["mole", "molar", "mass"]
EquilibriumSolver: TypeAlias = Literal["element_potential", "gibbs", "vcs", "auto"]
LogLevel: TypeAlias = Literal[0, 1, 2, 3, 4, 5]
LogLevel7: TypeAlias = Literal[0, 1, 2, 3, 4, 5, 6, 7]
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
CompositionLike: TypeAlias = str | dict[str, float] | ArrayLike
StateSetter: TypeAlias = tuple[float, float, CompositionLike]
ArrayCompositionLike: TypeAlias = str | dict[str, ArrayLike] | ArrayLike
ArrayStateSetter: TypeAlias = tuple[ArrayLike, ArrayLike, ArrayCompositionLike]

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

# Convenience functions provided in _types.py:
_P = ParamSpec("_P")

_TSelf = TypeVar("_TSelf")
_TReturn = TypeVar("_TReturn")
_T0 = TypeVar("_T0")
_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")

@overload
def add_args_to_signature(
    to_signature: Callable[Concatenate[_TSelf, _P], _TReturn], new_arg_type: type[_T0]
) -> Callable[
    [Callable[..., _TReturn]], Callable[Concatenate[_TSelf, _T0, _P], _TReturn]
]: ...
@overload
def add_args_to_signature(
    to_signature: Callable[Concatenate[_TSelf, _P], _TReturn],
    new_arg_type0: type[_T0],
    new_arg_type1: type[_T1],
) -> Callable[
    [Callable[..., _TReturn]], Callable[Concatenate[_TSelf, _T0, _T1, _P], _TReturn]
]: ...
@overload
def add_args_to_signature(
    to_signature: Callable[Concatenate[_TSelf, _P], _TReturn],
    new_arg_type0: type[_T0],
    new_arg_type1: type[_T1],
    new_arg_type2: type[_T2],
) -> Callable[
    [Callable[..., _TReturn]], Callable[Concatenate[_TSelf, _T0, _T1, _P], _TReturn]
]: ...
def literal_type_guard(tag: str, literal: TypeForm[_T0]) -> TypeGuard[_T0]: ...  # type: ignore[valid-type]
