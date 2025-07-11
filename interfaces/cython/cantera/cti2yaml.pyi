# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from argparse import ArgumentParser
from collections import OrderedDict
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Literal, TypeAlias, TypedDict

from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.representer import SafeRepresenter

yaml_version: tuple[int, int, int]
yaml_min_version: tuple[int, int, int]

class InputError(Exception):
    def __init__(self, msg: str, *args: Any) -> None: ...

BlockMap: CommentedMap

def FlowMap(*args: Any, **kwargs: Any) -> CommentedMap: ...
def FlowList(*args: Any, **kwargs: Any) -> CommentedSeq: ...
def float2string(data: float) -> str: ...
def represent_float(self: SafeRepresenter, data: float) -> str: ...
def applyUnits(value: float | tuple[float, str]) -> float | str: ...

OldKineticsModel: TypeAlias = Literal["GasKinetics", "Interface", "Edge"]
KineticsModel: TypeAlias = Literal["gas", "surface", "edge"]
OldTransportModel: TypeAlias = Literal["Mix", "Multi", "Ion"]
TransportModel: TypeAlias = Literal["mixture-averaged", "multicomponent", "ionized-gas"]
OldConcentrationBasis: TypeAlias = Literal["molar_volume", "molar_volume", "unity"]
ConcentrationBasis: TypeAlias = Literal[
    "species-molar-volume", "solvent-molar-volume", "unity"
]

OneAtm: float
OneBar: float
eV: float
ElectronMass: float

def enable_motz_wise() -> None: ...
def disable_motz_wise() -> None: ...
def validate(species: str = "yes", reactions: str = "yes") -> None: ...
def dataset(nm: str) -> None: ...
def standard_pressure(p0: float) -> None: ...
def units(
    length: str = "",
    quantity: str = "",
    mass: str = "",
    time: str = "",
    act_energy: str = "",
    energy: str = "",
    pressure: str = "",
) -> None: ...
def get_composition(
    atoms: dict[str, float] | str,
) -> dict[str, float] | OrderedDict[str, float]: ...

class element:
    symbol: str
    atomic_weight: float
    atomic_number: int | None
    def __init__(
        self,
        symbol: str = "",
        atomic_mass: float = 0.01,
        atomic_number: int | None = None,
    ) -> None: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: element) -> CommentedMap: ...

class RkPure(TypedDict):
    a: tuple[float, str] | float
    b: float

class RkBinary(TypedDict):
    a: float
    b: float

class species:
    name: str
    atoms: OrderedDict[str, float]
    size: float
    comment: str
    transport: gas_transport | None
    standard_state: constantIncompressible | None
    rk_pure: RkPure
    rk_binary: RkBinary
    density: float | str | None
    def __init__(
        self,
        name: str,
        atoms: str = "",
        note: str = "",
        thermo: thermo | Sequence[thermo] | None = None,
        transport: gas_transport | None = None,
        charge: float | None = None,
        size: float = 1.0,
        standardState: constantIncompressible | None = None,
    ) -> None: ...
    thermo: thermo  # Note: If this is moved above __init__, it overrides the class `thermo` for the static type checker
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: species) -> CommentedMap: ...

class thermo:
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: thermo) -> CommentedMap: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class NASA(thermo):
    model: Literal["NASA7"]
    T_range: Sequence[float]
    pref: float
    coeffs: Sequence[float]
    def __init__(
        self,
        Trange: Sequence[float] = (0.0, 0.0),
        coeffs: Sequence[float] = (),
        p0: float | None = None,
    ) -> None: ...

class NASA9(thermo):
    model: Literal["NASA9"]
    T_range: Sequence[float]
    pref: float
    coeffs: Sequence[float]
    def __init__(
        self,
        Trange: Sequence[float] = (0.0, 0.0),
        coeffs: Sequence[float] = (),
        p0: float | None = None,
    ) -> None: ...

class MultiPolyThermo(thermo):
    pref: float
    Tranges: Sequence[float]
    model: str
    data: Sequence[float]
    def __init__(self, regions: Sequence[thermo]) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class Shomate(thermo):
    model: Literal["Shomate"]
    T_range: Sequence[float]
    pref: float
    coeffs: Sequence[float]
    def __init__(
        self,
        Trange: Sequence[float] = (0.0, 0.0),
        coeffs: Sequence[float] = (),
        p0: float | None = None,
    ) -> None: ...

class const_cp(thermo):
    model: Literal["constant-cp"]
    pref: None
    t0: float
    h0: float
    s0: float
    cp0: float
    tmin: float
    tmax: float
    def __init__(
        self,
        t0: float | None = None,
        cp0: float | None = None,
        h0: float | None = None,
        s0: float | None = None,
        tmax: float | None = None,
        tmin: float | None = None,
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class gas_transport:
    geometry: Literal["atom", "linear", "nonlinear"]
    diameter: float
    well_depth: float
    dipole: float
    polarizability: float
    rot_relax: float
    acentric_factor: float
    disp_coeff: float
    quad_polar: float
    def __init__(
        self,
        geom: Literal["atom", "linear", "nonlinear"],
        diam: float,
        well_depth: float,
        dipole: float = 0.0,
        polar: float = 0.0,
        rot_relax: float = 0.0,
        acentric_factor: float | None = None,
        disp_coeff: float = 0.0,
        quad_polar: float = 0.0,
    ) -> None: ...
    @classmethod
    def to_yaml(
        cls, representer: SafeRepresenter, node: gas_transport
    ) -> CommentedMap: ...

CoverageParameters: TypeAlias = Sequence[str | float]

class Arrhenius:
    A: float
    b: float
    E: float
    coverage: list[CoverageParameters] | None
    def __init__(
        self,
        A: float = 0.0,
        b: float = 0.0,
        E: float = 0.0,
        coverage: CoverageParameters | list[CoverageParameters] = (),
    ) -> None: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: Arrhenius) -> CommentedMap: ...

class stick(Arrhenius):
    motz_wise: bool | None
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...

ReactionOptions: TypeAlias = Literal[
    "duplicate", "negative_A", "negative_orders", "nonreactant_orders"
]

class reaction:
    equation: str
    order: dict[str, float] | OrderedDict[str, float]
    number: int
    id: str
    options: Sequence[ReactionOptions]
    kf: Arrhenius
    type: str = "elementary"
    def __init__(
        self,
        equation: str,
        kf: Sequence[float] | Arrhenius,
        id: str = "",
        order: str = "",
        options: ReactionOptions | Sequence[ReactionOptions] = (),
    ) -> None: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: reaction) -> CommentedMap: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class three_body_reaction(reaction):
    type: Literal["three-body"]
    efficiencies: dict[str, float] | OrderedDict[str, float]
    def __init__(
        self,
        equation: str,
        kf: Sequence[float] | Arrhenius,
        efficiencies: str = "",
        id: str = "",
        options: ReactionOptions | Sequence[ReactionOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

FalloffFunction: TypeAlias = Troe | SRI | Lindemann

class falloff_base(reaction):
    k_low: Arrhenius
    k_high: Arrhenius
    falloff: FalloffFunction | None
    efficiencies: dict[str, float] | OrderedDict[str, float]
    def __init__(
        self,
        equation: str,
        klow: Sequence[float] | Arrhenius,
        khigh: Sequence[float] | Arrhenius,
        efficiencies: str,
        falloff: FalloffFunction | None,
        id: str,
        options: ReactionOptions | Sequence[ReactionOptions],
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class falloff_reaction(falloff_base):
    type: Literal["falloff"]
    def __init__(
        self,
        equation: str,
        kf0: Sequence[float] | Arrhenius,
        kf: Sequence[float] | Arrhenius,
        efficiencies: str = "",
        falloff: FalloffFunction | None = None,
        id: str = "",
        options: ReactionOptions | Sequence[ReactionOptions] = (),
    ) -> None: ...

class chemically_activated_reaction(falloff_base):
    type: Literal["chemically-activated"]
    def __init__(
        self,
        equation: str,
        kLow: Sequence[float] | Arrhenius,
        kHigh: Sequence[float] | Arrhenius,
        efficiencies: str = "",
        falloff: FalloffFunction | None = None,
        id: str = "",
        options: ReactionOptions | Sequence[ReactionOptions] = (),
    ) -> None: ...

class pdep_arrhenius(reaction):
    arrhenius: tuple[Sequence[float]]
    type: Literal["pressure-dependent-Arrhenius"]
    def __init__(
        self, equation: str, *args: Sequence[float], **kwargs: Any
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class chebyshev_reaction(reaction):
    type: Literal["Chebyshev"]
    Pmin: tuple[float, str]
    Pmax: tuple[float, str]
    Tmin: float
    Tmax: float
    coeffs: Sequence[Sequence[float]]
    def __init__(
        self,
        equation: str,
        Tmin: float = 300.0,
        Tmax: float = 2500.0,
        Pmin: tuple[float, str] = (0.001, "atm"),
        Pmax: tuple[float, str] = (100.0, "atm"),
        coeffs: Sequence[Sequence[float]] = (),
        **kwargs: Any,
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class surface_reaction(reaction):
    type: Literal["surface", "edge"] = "surface"
    sticking: bool
    beta: float | None
    rate_coeff_type: Literal["", "exchangecurrentdensity"]
    def __init__(
        self,
        equation: str,
        kf: Sequence[float] | Arrhenius,
        id: str = "",
        order: str = "",
        beta: float | None = None,
        options: ReactionOptions | Sequence[ReactionOptions] = (),
        rate_coeff_type: Literal["", "exchangecurrentdensity"] = "",
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class edge_reaction(surface_reaction):
    type: Literal["edge"]
    def __init__(
        self,
        equation: str,
        kf: Sequence[float] | Arrhenius,
        id: str = "",
        order: str = "",
        beta: float | None = None,
        options: ReactionOptions | Sequence[ReactionOptions] = (),
        rate_coeff_type: Literal["", "exchangecurrentdensity"] = "",
    ) -> None: ...

class state:
    t: float
    rho: float
    p: float
    X: str
    Y: str
    coverages: str
    molalities: str
    def __init__(
        self,
        temperature: float | None = None,
        pressure: float | None = None,
        mole_fractions: str | None = None,
        mass_fractions: str | None = None,
        density: float | None = None,
        coverages: str | None = None,
        solute_molalities: str | None = None,
    ) -> None: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: state) -> CommentedMap: ...

PhaseOptions: TypeAlias = Literal[
    "skip_undeclared_elements", "skip_undeclared_third_bodies"
]

class phase:
    name: str
    elements: str
    species: list[tuple[str, CommentedSeq]]
    reactions: list[list[str]]
    thermo_model: str | None = None
    kinetics: KineticsModel | None = None
    transport: TransportModel | None = None
    comment: str
    options: Sequence[PhaseOptions]
    initial_state: state | None
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        reactions: str | Sequence[str] = "none",
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    @classmethod
    def to_yaml(cls, representer: SafeRepresenter, node: phase) -> CommentedMap: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class ideal_gas(phase):
    kinetics: KineticsModel | None
    transport: TransportModel | None
    thermo_model: Literal["ideal-gas"]
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        reactions: str | Sequence[str] = "none",
        kinetics: OldKineticsModel | Literal["None"] = "GasKinetics",
        transport: OldTransportModel | Literal["None"] | None = None,
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...

class stoichiometric_solid(phase):
    thermo_model: Literal["fixed-stoichiometry"]
    density: float
    transport: TransportModel | None
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        density: float | None = None,
        transport: OldTransportModel | Literal["None"] = "None",
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class stoichiometric_liquid(stoichiometric_solid): ...

class metal(phase):
    thermo_model: Literal["electron-cloud"]
    density: float
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        density: float = -1.0,
        transport: OldTransportModel | Literal["None"] = "None",
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class liquid_vapor(phase):
    pure_fluids: dict[int, str]
    thermo_model: Literal["pure-fluid"]
    substance_flag: int
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        substance_flag: int = 0,
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class pureFluidParameters:
    species: str
    a_coeff: Sequence[float]
    b_coeff: float
    def __init__(
        self,
        species: str | None = None,
        a_coeff: Sequence[float] = (),
        b_coeff: float = 0,
    ) -> None: ...

class crossFluidParameters:
    species1: str
    species2: str
    a_coeff: Sequence[float]
    b_coeff: Sequence[float]
    def __init__(
        self,
        species: str | None = None,
        a_coeff: Sequence[float] = (),
        b_coeff: Sequence[float] = (),
    ) -> None: ...

class RedlichKwongMFTP(phase):
    thermo_model: Literal["Redlich-Kwong"]
    kinetics: KineticsModel | None
    transport: TransportModel | None
    activity_coefficients: Sequence[pureFluidParameters | crossFluidParameters]
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        reactions: str | Sequence[str] = "none",
        kinetics: OldKineticsModel = "GasKinetics",
        initial_state: state | None = None,
        activity_coefficients: Sequence[pureFluidParameters | crossFluidParameters]
        | None = None,
        transport: OldTransportModel | Literal["None"] = "None",
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class constantIncompressible:
    molar_volume: float
    def __init__(self, molarVolume: float = 0.0) -> None: ...

class IdealSolidSolution(phase):
    thermo_model: Literal["ideal-condensed", "binary-solution-tabulated"] = (
        "ideal-condensed"
    )
    standard_concentration: ConcentrationBasis
    transport: TransportModel | None
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        transport: OldTransportModel | Literal["None"] = "None",
        initial_state: state | None = None,
        standard_concentration: OldConcentrationBasis | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class table:
    x: tuple[Sequence[float], str]
    h: tuple[Sequence[float], str]
    s: tuple[Sequence[float], str]
    def __init__(
        self,
        moleFraction: tuple[Sequence[float], str] = ([], ""),
        enthalpy: tuple[Sequence[float], str] = ([], ""),
        entropy: tuple[Sequence[float], str] = ([], ""),
    ) -> None: ...

class BinarySolutionTabulatedThermo(IdealSolidSolution):
    thermo_model: Literal["binary-solution-tabulated"]
    tabulated_species: str
    tabulated_thermo: table
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        transport: OldTransportModel | Literal["None"] = "None",
        initial_state: state | None = None,
        standard_concentration: OldConcentrationBasis | None = None,
        tabulated_species: str | None = None,
        tabulated_thermo: table | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class lattice(phase):
    thermo_model: Literal["lattice"]
    site_density: float
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str | Sequence[str] = "",
        note: str = "",
        reactions: str | Sequence[str] = "none",
        transport: OldTransportModel | Literal["None"] = "None",
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
        site_density: float | None = None,
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class ideal_interface(phase):
    thermo_model: Literal["ideal-surface", "edge"] = "ideal-surface"
    kinetics: KineticsModel
    transport: TransportModel
    site_density: float
    adjacent_phases: list[str]
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str = "",
        note: str = "",
        reactions: str = "none",
        site_density: float = 0.0,
        phases: str
        | Sequence[str] = (),  # Note: Does not actually accept Sequence input
        kinetics: OldKineticsModel | Literal["None"] = "Interface",
        transport: OldTransportModel | Literal["None"] = "None",
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class edge(ideal_interface):
    thermo_model: Literal["edge"]
    def __init__(
        self,
        name: str = "",
        elements: str = "",
        species: str = "",
        note: str = "",
        reactions: str = "none",
        site_density: float = 0.0,
        phases: str
        | Sequence[str] = (),  # Note: Does not actually accept Sequence input
        kinetics: OldKineticsModel | Literal["None"] = "Edge",
        transport: OldTransportModel | Literal["None"] = "None",
        initial_state: state | None = None,
        options: PhaseOptions | Sequence[PhaseOptions] = (),
    ) -> None: ...

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
    def get_yaml(self, out: CommentedMap) -> None: ...

class SRI:
    A: float
    B: float
    C: float
    D: float | None
    E: float | None
    def __init__(
        self,
        A: float = 0.0,
        B: float = 0.0,
        C: float = 0.0,
        D: float | None = None,
        E: float | None = None,
    ) -> None: ...
    def get_yaml(self, out: CommentedMap) -> None: ...

class Lindemann:
    def get_yaml(self, out: CommentedMap) -> None: ...

# Note: Many more encodings available, but you probably don't want most of them.
# See: https://docs.python.org/3/library/codecs.html#standard-encodings
Encoding: TypeAlias = Literal["utf-8", "latin-1", "ascii"]

def convert(
    filename: Path | str | None = None,
    output_name: str | None = None,
    text: str | None = None,
    encoding: Encoding = "latin-1",
) -> tuple[int, int, list[ideal_interface], Path]: ...
def create_argparser() -> ArgumentParser: ...
def main() -> None: ...
