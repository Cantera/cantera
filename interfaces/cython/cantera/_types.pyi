# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from types import EllipsisType
from typing import Literal, TypeAlias, TypedDict

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__: list[str] = [
    "Array",
    "ArrayLike",
    "Basis",
    "CompositionLike",
    "CompositionVariable",
    "CompressionLevel",
    "EquilibriumSolver",
    "FullState",
    "Index",
    "LogLevel",
    "LogLevel7",
    "StateDefinition",
    "StateSetter",
    "StateVariable",
]

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
ArrayCompositionLike: TypeAlias = str | dict[str, Array | float] | ArrayLike
ArrayStateSetter: TypeAlias = tuple[Array | float, Array | float, ArrayCompositionLike]

# PureFluid state definitions
PureFluidStateVariable: TypeAlias = Literal[StateVariable, "Q"]
PureFluidPropertyPair: TypeAlias = Literal[
    PropertyPair, "TQ", "PQ", "ST", "TV", "PV", "UP", "VH", "TH", "SH"
]
PureFluidFullState: TypeAlias = Literal[
    FullState, "TDQ", "TPQ", "UVQ", "DPQ", "HPQ", "SPQ", "SVQ"
]
PureFluidStateSetter: TypeAlias = tuple[float, float, float]
ArrayPureFluidStateSetter: TypeAlias = tuple[
    Array | float, Array | float, Array | float
]

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
