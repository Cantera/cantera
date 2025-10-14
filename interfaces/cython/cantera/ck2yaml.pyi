# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import logging
from argparse import ArgumentParser
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, ClassVar, Concatenate, Literal

from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.nodes import MappingNode, ScalarNode
from ruamel.yaml.representer import SafeRepresenter
from typing_extensions import Never

from ._types import Array

yaml_version: tuple[int, int, int]
yaml_min_version: tuple[int, int, int]
BlockMap: type[CommentedMap]

class ErrorFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str: ...

logger: logging.Logger

def FlowMap(*args: Any, **kwargs: Any) -> CommentedMap: ...
def FlowList(*args: Any, **kwargs: Any) -> CommentedSeq: ...
def represent_float(self: SafeRepresenter, data: float) -> ScalarNode: ...

QUANTITY_UNITS: dict[str, str]
ENERGY_UNITS: dict[str, str]

def strip_nonascii(s: str) -> str: ...
def compatible_quantities(
    quantity_basis: Literal["mol", "molec"], units: str
) -> bool: ...

class InputError(Exception):
    def __init__(self, message: str) -> None: ...

class ParserLogger(logging.Handler):
    parser: Parser
    errors: list[str]
    def __init__(self, parser: Parser) -> None: ...
    def emit(self, record: logging.LogRecord) -> None: ...

class Species:
    label: str
    thermo: Nasa7 | Nasa9 | None
    transport: TransportData | None
    sites: float | None
    composition: dict[str, float] | None
    note: str | None
    def __init__(self, label: str, sites: float | None = None) -> None: ...
    def __str__(self) -> str: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Species) -> dict[str, str]: ...

class Nasa7:
    Tmin: float
    Tmax: float
    Tmid: float | None
    low_coeffs: Sequence[float]
    high_coeffs: Sequence[float] | None
    note: str
    def __init__(
        self,
        *,
        Tmin: float,
        Tmax: float,
        Tmid: float | None,
        low_coeffs: Sequence[float],
        high_coeffs: Sequence[float],
        note: str = "",
    ) -> None: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Nasa7) -> MappingNode: ...

class Nasa9:
    note: str
    data: list[tuple[list[float], list[float]]]
    Tranges: list[float]
    def __init__(
        self,
        *,
        parser: Parser,
        data: list[tuple[list[float], list[float]]],
        note: str = "",
    ) -> None: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Nasa9) -> MappingNode: ...

class Reaction:
    index: int
    reactants: list[tuple[float, Species]]
    products: list[tuple[float, Species]]
    kinetics: KineticsModel
    reversible: bool
    duplicate: bool
    forward_orders: dict[str, float]
    third_body: str
    comment: str
    def __init__(
        self,
        index: int = -1,
        reactants: list[tuple[float, Species]] | None = None,
        products: list[tuple[float, Species]] | None = None,
        kinetics: KineticsModel | None = None,
        reversible: bool = True,
        duplicate: bool = False,
        forward_orders: dict[str, float] | None = None,
        third_body: str | None = None,
    ) -> None: ...
    def _coeff_string(self, coeffs: list[tuple[float, Species]]) -> str: ...
    def __str__(self) -> str: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Reaction) -> MappingNode: ...

class KineticsModel:
    pressure_dependent: ClassVar[bool | None] = None
    efficiencies: dict[str, float]
    def __init__(self) -> None: ...
    def reaction_string_suffix(self, species: Species) -> Literal[""]: ...
    def reduce(self, output: CommentedMap) -> Never: ...

class Arrhenius:
    A: tuple[float, str]
    b: float
    Ea: tuple[float, str]
    parser: Parser
    def __init__(
        self,
        A: tuple[float, str] | float = 0.0,
        b: float = 0.0,
        Ea: tuple[float, str] | float = 0.0,
        *,
        parser: Parser,
    ) -> None: ...
    def as_yaml(self, extra: Any = ()) -> CommentedMap: ...

class ElementaryRate(KineticsModel):
    pressure_dependent: ClassVar[bool | None] = False
    rate: Arrhenius
    def __init__(self, rate: Arrhenius, **kwargs: Any) -> None: ...
    def reduce(self, output: CommentedMap) -> None: ...  # type: ignore[override]

class SurfaceRate(KineticsModel):
    pressure_dependent: ClassVar[bool | None] = False
    rate: Arrhenius
    coverages: list[tuple[str, float, float, float]]
    is_sticking: bool
    motz_wise: bool
    def __init__(
        self,
        *,
        rate: Arrhenius,
        coverages: list[tuple[str, float, float, float]],
        is_sticking: bool,
        motz_wise: bool,
        **kwargs: Any,
    ) -> None: ...
    def reduce(self, output: CommentedMap) -> None: ...  # type: ignore[override]

class PDepArrhenius(KineticsModel):
    pressure_dependent: ClassVar[bool | None] = True
    pressures: list[float]
    pressure_units: str
    arrhenius: list[Arrhenius]
    def __init__(
        self,
        *,
        pressures: list[float],
        pressure_units: str,
        arrhenius: list[Arrhenius],
        **kwargs: Any,
    ) -> None: ...
    def reduce(self, output: CommentedMap) -> None: ...  # type: ignore[override]

class Chebyshev(KineticsModel):
    pressure_dependent: ClassVar[bool | None] = True
    Tmin: float
    Tmax: float
    Pmin: tuple[float, str]
    Pmax: tuple[float, str]
    coeffs: Array
    quantity_units: str
    def __init__(
        self,
        coeffs: Array,
        *,
        Tmin: float,
        Tmax: float,
        Pmin: tuple[float, str],
        Pmax: tuple[float, str],
        quantity_units: str,
        **kwargs: Any,
    ) -> None: ...
    def reduce(self, output: CommentedMap) -> None: ...  # type: ignore[override]

class ThreeBody(KineticsModel):
    pressure_dependent: ClassVar[bool | None] = True
    high_rate: Arrhenius
    efficiencies: dict[str, float]
    def __init__(
        self,
        high_rate: Arrhenius | None = None,
        efficiencies: dict[str, float] | None = None,
        **kwargs: Any,
    ) -> None: ...
    def reaction_string_suffix(self, species: Species) -> Literal[" + M"]: ...  # type: ignore[override]
    def reduce(self, output: CommentedMap) -> None: ...  # type: ignore[override]

class Falloff(ThreeBody):
    low_rate: Arrhenius
    high_rate: Arrhenius
    efficiencies: dict[str, float]
    F: Troe | Sri
    def __init__(
        self,
        low_rate: Arrhenius | None = None,
        high_rate: Arrhenius | None = None,
        efficiencies: dict[str, float] | None = None,
        F: Troe | Sri | None = None,
        **kwargs: Any,
    ) -> None: ...
    def reaction_string_suffix(self, species: Species) -> str: ...  # type: ignore[override]
    def reduce(self, output: CommentedMap) -> None: ...  # type: ignore[override]

class ChemicallyActivated(ThreeBody):
    low_rate: Arrhenius
    high_rate: Arrhenius
    efficiencies: dict[str, float]
    F: Troe | Sri
    def __init__(
        self,
        low_rate: Arrhenius | None = None,
        high_rate: Arrhenius | None = None,
        efficiencies: dict[str, float] | None = None,
        F: Troe | Sri | None = None,
        **kwargs: Any,
    ) -> None: ...
    def reaction_string_suffix(self, species: Species) -> str: ...  # type: ignore[override]
    def reduce(self, output: CommentedMap) -> None: ...  # type: ignore[override]

class Troe:
    A: float
    T3: float
    T1: float
    T2: float | None
    def __init__(
        self,
        A: float = 0.0,
        T3: float = 0.0,
        T1: float = 0.0,
        T2: float | None = None,
    ) -> None: ...
    def reduce(self, output: CommentedMap) -> None: ...

class Sri:
    A: float
    B: float
    C: float
    D: float | None
    E: float | None
    def __init__(
        self,
        *,
        A: float,
        B: float,
        C: float,
        D: float | None = None,
        E: float | None = None,
    ) -> None: ...
    def reduce(self, output: CommentedMap) -> None: ...

class TransportData:
    geometry_flags: ClassVar[list[str]] = ["atom", "linear", "nonlinear"]
    geometry: Literal["atom", "linear", "nonlinear"]
    well_depth: float
    collision_diameter: float
    dipole_moment: float
    polarizability: float
    z_rot: float
    note: str
    def __init__(
        self,
        parser: Parser,
        label: str,
        geometry: Literal[0, 1, 2, "0", "1", "2"],
        well_depth: str | float,
        collision_diameter: str | float,
        dipole_moment: str | float,
        polarizability: str | float,
        z_rot: str | float,
        note: str = "",
    ) -> None: ...
    @classmethod
    def to_yaml(
        cls, representer: SafeRepresenter, node: TransportData
    ) -> MappingNode: ...

def fortFloat(s: str) -> float: ...
def get_index(seq: str | Sequence[str], value: str) -> int | None: ...
def contains(seq: str | Sequence[str], value: str) -> bool: ...

class Surface:
    name: str
    site_density: float
    species_list: list[Species]
    reactions: list[Reaction]
    def __init__(self, name: str, site_density: float) -> None: ...

class Parser:
    max_loglevel: int
    handler: ParserLogger
    processed_units: bool
    energy_units: str
    output_energy_units: str
    quantity_units: str
    output_quantity_units: str
    motz_wise: bool | None
    single_intermediate_temperature: bool
    skip_undeclared_species: bool
    exit_on_error: bool
    permissive: bool
    verbose: bool
    elements: list[str]
    element_weights: dict[str, float]
    species_list: list[Species]
    species_dict: dict[str, Species]
    surfaces: list[Surface]
    reactions: list[Reaction]
    header_lines: list[str]
    extra: CommentedMap
    files: list[str]
    raw_lines: list[str]
    current_range: list[int]
    def __init__(self) -> None: ...
    def __del__(self) -> None: ...
    def entry(self, where: str = "entry", kind: str = "Error") -> str: ...
    @staticmethod
    def parse_composition(
        elements: str, nElements: int, width: int
    ) -> dict[str, int]: ...
    @staticmethod
    def get_rate_constant_units(
        length_dims: int,
        length_units: str,
        quantity_dims: int,
        quantity_units: str,
        time_dims: int = 1,
        time_units: str = "s",
    ) -> str: ...
    def read_NASA7_entry(
        self, lines: list[str], TintDefault: float, comments: list[str]
    ) -> tuple[str, Nasa7, dict[str, float]]: ...
    def parse_nasa7_section(self, lines: list[tuple[int, str, str, str]]) -> None: ...
    def read_NASA9_entry(
        self, entry: list[str], comments: list[str]
    ) -> tuple[str, Nasa9, dict[str, float]]: ...
    def read_kinetics_entry(
        self, entry: list[str], surface: bool
    ) -> tuple[Reaction, Reaction | None]: ...
    def load_extra_file(self, path: Path | str) -> None: ...
    def load_data_file(
        self,
        path: Path | str,
        load_method: Callable[Concatenate[str, ...], None],
        kind: str,
        **kwargs: Any,
    ) -> None: ...
    def load_chemkin_file(self, path: Path | str, surface: bool = False) -> None: ...
    def parse_elements_section(
        self, lines: list[tuple[int, str, str, str]]
    ) -> None: ...
    def parse_species_section(self, lines: list[tuple[int, str, str, str]]) -> None: ...
    def parse_site_section(self, lines: list[tuple[int, str, str, str]]) -> None: ...
    def parse_nasa9_section(self, lines: list[tuple[int, str, str, str]]) -> None: ...
    species_tokens: set[str]
    other_tokens: dict[str, str]
    Slen: int
    def parse_reactions_section(
        self, lines: list[tuple[int, str, str, str]], surface: bool
    ) -> None: ...
    def load_transport_file(self, path: Path | str) -> None: ...
    def parse_transport_section(
        self, lines: list[tuple[int, str, str, str]]
    ) -> None: ...
    def write_yaml(
        self, name: str = "gas", out_name: str = "mech.yaml"
    ) -> list[str]: ...
    @staticmethod
    def convert_mech(
        input_file: Path | str,
        thermo_file: Path | str | None = None,
        transport_file: Path | str | None = None,
        surface_file: Path | str | None = None,
        phase_name: str = "gas",
        extra_file: str | None = None,
        out_name: str | None = None,
        single_intermediate_temperature: bool = False,
        quiet: bool = False,
        permissive: bool | None = None,
        verbose: bool = False,
        exit_on_error: bool = False,
    ) -> tuple[Parser, list[str]]: ...
    def show_duplicate_reactions(self, error_message: str) -> None: ...

def convert(
    input_file: Path | str,
    thermo_file: Path | str | None = None,
    transport_file: Path | str | None = None,
    surface_file: Path | str | None = None,
    *,
    phase_name: str = "gas",
    extra_file: str | None = None,
    out_name: Path | str | None = None,
    single_intermediate_temperature: bool = False,
    quiet: bool = False,
    permissive: bool | None = None,
    verbose: bool = False,
) -> list[str]: ...
def main(argv: Sequence[str] | None = None) -> None: ...
def create_argparser() -> ArgumentParser: ...
