# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from typing import Literal, TypeAlias, TypedDict
from typing_extensions import Never

from ._types import Array, ArrayLike
from .solutionbase import _SolutionBase

# Avoid fixed options unless we can find a way to support custom extensions:
# _TransportModel: TypeAlias = Literal[
#     "none",
#     "unity-Lewis-number",
#     "mixture-averaged",
#     "mixture-averaged-CK",
#     "multicomponent",
#     "multicomponent-CK",
#     "ionized-gas",
#     "water",
#     "high-pressure",
#     "high-pressure=Chung",
# ]
_TransportModel: TypeAlias = str

_GeometryOptions: TypeAlias = Literal["atom", "linear", "nonlinear"]

class _GasTransportInput(TypedDict, total=False):
    model: Literal["gas"]
    geometry: _GeometryOptions
    diameter: float
    well_depth: float
    dipole: float
    polarizability: float
    rotational_relaxation: float
    acentric_factor: float
    dispersion_coefficient: float
    quadrupole_polarizability: float
    note: str

_TransportFittingErrors = TypedDict(
    "_TransportFittingErrors",
    {
        "viscosity-max-abs-error": float,
        "viscosity-max-rel-error": float,
        "conductivity-max-abs-error": float,
        "conductivity-max-rel-error": float,
        "diff-coeff-max-abs-error": float,
        "diff-coeff-max-rel-error": float,
    },
)

class GasTransportData:
    def __init__(
        self,
        geometry: _GeometryOptions | Literal[""] = "",
        diameter: float = -1,
        well_depth: float = -1,
        dipole: float = 0.0,
        polarizability: float = 0.0,
        rotational_relaxation: float = 0.0,
        acentric_factor: float = 0.0,
        dispersion_coefficient: float = 0.0,
        quadrupole_polarizability: float = 0.0,
        *,
        init: bool = True,
    ) -> None: ...
    def set_customary_units(
        self,
        geometry: str,
        diameter: float,
        well_depth: float,
        dipole: float = 0.0,
        polarizability: float = 0.0,
        rotational_relaxation: float = 0.0,
        acentric_factor: float = 0.0,
        dispersion_coefficient: float = 0.0,
        quadrupole_polarizability: float = 0.0,
    ) -> None: ...
    @property
    def input_data(self) -> _GasTransportInput: ...
    def update_user_data(self, data: GasTransportData) -> None: ...
    def clear_user_data(self) -> None: ...
    @property
    def geometry(self) -> _GeometryOptions: ...
    @geometry.setter
    def geometry(self, geometry: _GeometryOptions) -> None: ...
    @property
    def diameter(self) -> float: ...
    @diameter.setter
    def diameter(self, diameter: float) -> None: ...
    @property
    def well_depth(self) -> float: ...
    @well_depth.setter
    def well_depth(self, well_depth: float) -> None: ...
    @property
    def dipole(self) -> float: ...
    @dipole.setter
    def dipole(self, dipole: float) -> None: ...
    @property
    def polarizability(self) -> float: ...
    @polarizability.setter
    def polarizability(self, polarizability: float) -> None: ...
    @property
    def rotational_relaxation(self) -> float: ...
    @rotational_relaxation.setter
    def rotational_relaxation(self, rotational_relaxation: float) -> None: ...
    @property
    def acentric_factor(self) -> float: ...
    @acentric_factor.setter
    def acentric_factor(self, acentric_factor: float) -> None: ...
    @property
    def dispersion_coefficient(self) -> float: ...
    @dispersion_coefficient.setter
    def dispersion_coefficient(self, dispersion_coefficient: float) -> None: ...
    @property
    def quadrupole_polarizability(self) -> float: ...
    @quadrupole_polarizability.setter
    def quadrupole_polarizability(self, quadrupole_polarizability: float) -> None: ...

class Transport(_SolutionBase):
    @property
    def transport_model(self) -> _TransportModel: ...
    @transport_model.setter
    def transport_model(self, model: _TransportModel) -> None: ...
    @property
    def CK_mode(self) -> bool: ...
    @property
    def viscosity(self) -> float: ...
    @property
    def species_viscosities(self) -> Array: ...
    @property
    def electrical_conductivity(self) -> float: ...
    @property
    def thermal_conductivity(self) -> float: ...
    @property
    def mix_diff_coeffs(self) -> Array: ...
    @property
    def mix_diff_coeffs_mass(self) -> Array: ...
    @property
    def mix_diff_coeffs_mole(self) -> Array: ...
    @property
    def thermal_diff_coeffs(self) -> Array: ...
    @property
    def multi_diff_coeffs(self) -> Array: ...
    @property
    def binary_diff_coeffs(self) -> Array: ...
    @property
    def mobilities(self) -> Array: ...
    def get_viscosity_polynomial(self, i: int) -> Array: ...
    def get_thermal_conductivity_polynomial(self, i: int) -> Array: ...
    def get_binary_diff_coeffs_polynomial(self, i: int, j: int) -> Array: ...
    def get_collision_integral_polynomials(
        self, i: int, j: int
    ) -> tuple[Array, Array, Array]: ...
    def set_viscosity_polynomial(self, i: int, values: ArrayLike) -> None: ...
    def set_thermal_conductivity_polynomial(
        self, i: int, values: ArrayLike
    ) -> None: ...
    def set_binary_diff_coeffs_polynomial(
        self, i: int, j: int, values: ArrayLike
    ) -> None: ...
    def set_collision_integral_polynomial(
        self,
        i: int,
        j: int,
        avalues: ArrayLike,
        bvalues: ArrayLike,
        cvalues: ArrayLike,
        actualT: bool = False,
    ) -> None: ...
    @property
    def transport_fitting_errors(self) -> _TransportFittingErrors: ...

class DustyGasTransport(Transport):
    @property
    def porosity(self) -> Never: ...
    @porosity.setter
    def porosity(self, value: float) -> None: ...
    @property
    def tortuosity(self) -> Never: ...
    @tortuosity.setter
    def tortuosity(self, value: float) -> None: ...
    @property
    def mean_pore_radius(self) -> Never: ...
    @mean_pore_radius.setter
    def mean_pore_radius(self, value: float) -> None: ...
    @property
    def mean_particle_diameter(self) -> Never: ...
    @mean_particle_diameter.setter
    def mean_particle_diameter(self, value: float) -> None: ...
    @property
    def permeability(self) -> Never: ...
    @permeability.setter
    def permeability(self, value: float) -> None: ...
    def molar_fluxes(
        self,
        T1: float,
        T2: float,
        rho1: float,
        rho2: float,
        Y1: ArrayLike,
        Y2: ArrayLike,
        delta: float,
    ) -> Array: ...
