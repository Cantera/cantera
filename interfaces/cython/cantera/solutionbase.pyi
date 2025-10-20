# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from collections.abc import Sequence
from pathlib import Path
from typing import (
    Any,
    Literal,
    TypeAlias,
    TypedDict,
    overload,
)

from typing_extensions import Never, Self

from ._types import Array, ArrayLike, Basis, CompressionLevel
from .kinetics import Kinetics, _KineticsType
from .reaction import Reaction
from .thermo import Species, ThermoPhase, _ThermoType
from .transport import _TransportModel
from .units import UnitSystem, _UnitDict, _UnitDictBytes

_SortingType: TypeAlias = Literal["alphabetical", "molar-mass"] | None

_YamlHeader = TypedDict(
    "_YamlHeader",
    {
        "description": str,
        "generator": str,
        "input-files": list[str],
        "cantera-version": str,
        "git-commit": str,
        "date": str,
    },
    total=False,
)

class _SolutionBase:
    def __init__(
        self,
        infile: Path | str = "",
        name: str = "",
        adjacent: Sequence[ThermoPhase] = (),
        *,
        origin: _SolutionBase | None = None,
        yaml: str | None = None,
        thermo: ThermoPhase | None = None,
        species: Sequence[Species] | None = None,
        kinetics: Kinetics | None = None,
        reactions: Sequence[Reaction] | None = (),
        init: bool = True,
        **kwargs: Any,
    ) -> None: ...
    @property
    def name(self) -> str: ...
    @name.setter
    def name(self, name: str) -> None: ...
    @property
    def source(self) -> str: ...
    @property
    def composite(
        self,
    ) -> tuple[_ThermoType | None, _KineticsType | None, _TransportModel | None]: ...
    @property
    def input_data(
        self,
    ) -> dict[str, str | list[str] | dict[str, float | dict[str, float]]]: ...
    @property
    def input_header(self) -> _YamlHeader: ...
    def update_user_data(self, data: dict[str, Any]) -> None: ...
    def clear_user_data(self) -> None: ...
    def update_user_header(self, data: dict[str, str | list[str]]) -> None: ...
    def clear_user_header(self) -> None: ...
    @overload
    def write_yaml(
        self,
        filename: None,
        phases: Sequence[ThermoPhase] | None = None,
        units: UnitSystem | _UnitDict | _UnitDictBytes | None = None,
        precision: int | None = None,
        skip_user_defined: bool | None = None,
        header: bool = True,
    ) -> str: ...
    @overload
    def write_yaml(
        self,
        filename: str | Path,
        phases: Sequence[ThermoPhase] | None = None,
        units: UnitSystem | _UnitDict | _UnitDictBytes | None = None,
        precision: int | None = None,
        skip_user_defined: bool | None = None,
        header: bool = True,
    ) -> None: ...
    @overload
    def write_yaml(
        self,
        filename: str | Path | None = None,
        phases: Sequence[ThermoPhase] | None = None,
        units: UnitSystem | _UnitDict | _UnitDictBytes | None = None,
        precision: int | None = None,
        skip_user_defined: bool | None = None,
        header: bool = True,
    ) -> str | None: ...
    def write_chemkin(
        self,
        mechanism_path: str | Path | None = None,
        thermo_path: str | Path | None = None,
        transport_path: str | Path | None = None,
        sort_species: _SortingType = None,
        sort_elements: _SortingType = None,
        overwrite: bool = False,
        quiet: bool = False,
    ) -> None: ...
    def __getitem__(self, selection: slice) -> Self: ...
    @property
    def selected_species(self) -> list[int]: ...
    @selected_species.setter
    def selected_species(
        self, species: str | int | Sequence[str] | Sequence[int]
    ) -> None: ...
    def __setstate__(self, pkl: str) -> None: ...
    def __copy__(self) -> Never: ...

class SolutionArrayBase:
    def __init__(
        self,
        phase: _SolutionBase,
        shape: int | tuple[int, ...] = (0,),
        states: ArrayLike | None = None,
        extra: str | Sequence[str] | dict[str, ArrayLike] | None = None,
        meta: dict[str, Any] = {},
        init: bool = True,
    ) -> None: ...
    @property
    def size(self) -> int: ...
    def _api_shape(self) -> tuple[int, ...]: ...
    def _set_api_shape(self, shape: Sequence[int]) -> None: ...
    @overload
    def info(
        self,
        keys: Sequence[str] | None = None,
        rows: int = 10,
        width: int | None = None,
        display: Literal[False] = False,
    ) -> str: ...
    @overload
    def info(
        self,
        keys: Sequence[str] | None = None,
        rows: int = 10,
        width: int | None = None,
        display: bool = True,
    ) -> None: ...
    @property
    def meta(self) -> dict[str, Any]: ...
    @meta.setter
    def meta(self, meta: dict[str, Any]) -> None: ...
    @property
    def extra(self) -> list[str]: ...
    @property
    def component_names(self) -> list[str]: ...
    def resize(self, size: int | tuple[int, ...]) -> None: ...
    def _has_component(self, name: str) -> bool: ...
    def _get_component(self, name: str) -> Array: ...
    def _set_component(self, name: str, data: Array) -> None: ...
    def _set_loc(self, loc: int) -> None: ...
    def _update_state(self, loc: int) -> None: ...
    def _get_state(self, loc: int) -> Array: ...
    def _set_state(self, loc: int, data: Array) -> None: ...
    def _has_extra(self, name: str) -> bool: ...
    def _add_extra(self, name: str, back: bool = True) -> None: ...
    def get_auxiliary(self, loc: int) -> dict[str, Any]: ...
    def set_auxiliary(self, loc: int, data: dict[str, Any]) -> None: ...
    def _append(self, state: Array, extra: dict[str, Any]) -> None: ...
    def _cxx_save(
        self,
        filename: str,
        name: str,
        sub: str,
        description: str,
        overwrite: bool,
        compression: CompressionLevel,
        basis: Basis,
    ) -> None: ...
    def _cxx_restore(self, filename: str, name: str, sub: str) -> dict[str, str]: ...
