# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3
# pyright: reportMissingImports=false, reportAttributeAccessIssue=false
# pyright: reportUndefinedVariable=false, reportUnboundVariable=false
# pyright: reportInvalidTypeArguments=false, reportAssignmentType=false
# pyright: reportIndexIssue=false, reportInvalidTypeForm=false

from pathlib import Path as _Path
from collections.abc import Sequence as _Sequence
from typing import (Any as _Any, Literal as _Literal, TypeAlias as _TypeAlias,
                    TypedDict as _TypedDict, TYPE_CHECKING)
from typing_extensions import Never as _Never

if TYPE_CHECKING:
    from .thermo import Species as _Species, ThermoPhase as _ThermoPhase
    from .reaction import Reaction as _Reaction

import numpy as np

import cython
import cython.cimports.numpy as cnp
from cython.cimports.cantera._utils import (stringify, pystr, anymap_to_py,
                                            py_to_anymap)
from .solutionbase import _SolutionBase
from cython.cimports.cantera.thermo import ThermoPhase

from ._types import Array as _Array, ArrayLike as _ArrayLike


# Avoid fixed options unless we can find a way to support custom extensions
_TransportModel: _TypeAlias = str

_GeometryOptions: _TypeAlias = _Literal["atom", "linear", "nonlinear"]

class _GasTransportInput(_TypedDict, total=False):
    model: _Literal["gas"]
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

_TransportFittingErrors = _TypedDict(
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


# NOTE: These cdef functions cannot be members of Transport because they would
# cause "layout conflicts" when creating derived classes with multiple bases,
# such as class Solution. [Cython 0.16]
@cython.cfunc
def get_transport_1d(tran: Transport, method: transportMethod1d) -> np.ndarray:
    data = np.empty(tran.thermo.nSpecies())
    cdata: cython.double[::1] = data
    method(tran.transport, span[cython.double](
        cython.address(cdata[0]), cython.cast(cython.size_t, data.size)))
    if tran._selected_species.size:
        return data[tran._selected_species]
    else:
        return data


@cython.cfunc
def get_transport_2d(tran: Transport, method: transportMethod2d) -> np.ndarray:
    kk: cython.size_t = tran.thermo.nSpecies()
    # getMultiDiffCoeffs returns an array using Fortran ordering
    data = np.empty((kk, kk), order='F')
    cdata: cython.double[::1, :] = data
    method(tran.transport, kk, span[cython.double](
        cython.address(cdata[0, 0]), cython.cast(cython.size_t, kk * kk)))
    return data


@cython.cfunc
def get_transport_polynomial(tran: Transport, method: transportPolyMethod1i,
                             index: cython.int, n_coeffs: cython.int) -> np.ndarray:
    data = np.empty(n_coeffs)
    cdata: cython.double[::1] = data
    method(tran.transport, index, span[cython.double](
        cython.address(cdata[0]), cython.cast(cython.size_t, data.size)))
    return data


@cython.cfunc
def get_binary_transport_polynomial(tran: Transport, method: transportPolyMethod2i,
                                    indexi: cython.int, indexj: cython.int,
                                    n_coeffs: cython.int) -> np.ndarray:
    data = np.empty(n_coeffs)
    cdata: cython.double[::1] = data
    method(tran.transport, indexi, indexj, span[cython.double](
        cython.address(cdata[0]), cython.cast(cython.size_t, data.size)))
    return data


@cython.cclass
class GasTransportData:
    """
    Transport data for a single gas-phase species which can be used in
    mixture-averaged or multicomponent transport models.

    The arguments passed to the constructor are equivalent to the properties of
    the object, with values in MKS units. To set properties in non-MKS units,
    use the `set_customary_units` method.
    """
    def __cinit__(self, geometry='', diameter=-1, well_depth=-1,
                  dipole=0.0, polarizability=0.0, rotational_relaxation=0.0,
                  acentric_factor=0.0, dispersion_coefficient=0.0,
                  quadrupole_polarizability=0.0, *, init=True):
        if init:
            self._data = make_shared[CxxGasTransportData](stringify(geometry),
                cython.cast(cython.double, diameter),
                cython.cast(cython.double, well_depth),
                cython.cast(cython.double, dipole),
                cython.cast(cython.double, polarizability),
                cython.cast(cython.double, rotational_relaxation),
                cython.cast(cython.double, acentric_factor),
                cython.cast(cython.double, dispersion_coefficient),
                cython.cast(cython.double, quadrupole_polarizability))
            self.data = cython.cast(cython.pointer(CxxGasTransportData),
                                    self._data.get(), typecheck=True)

    def __init__(
        self,
        geometry: _GeometryOptions | _Literal[""] = "",
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
    ) -> None:
        pass

    @cython.cfunc
    def _assign(self, other: shared_ptr[CxxTransportData]):
        self._data = other
        self.data = cython.cast(cython.pointer(CxxGasTransportData),
                                self._data.get(), typecheck=True)

    def set_customary_units(
        self,
        geometry: _GeometryOptions,
        diameter: float,
        well_depth: float,
        dipole: float = 0.0,
        polarizability: float = 0.0,
        rotational_relaxation: float = 0.0,
        acentric_factor: float = 0.0,
        dispersion_coefficient: float = 0.0,
        quadrupole_polarizability: float = 0.0
    ) -> None:
        """
        Set the parameters using "customary" units: diameter in Angstroms, well
        depth in Kelvin, dipole in Debye, and polarizability in Angstroms³.
        These are the units used in in CK-style input files.
        """
        self.data.setCustomaryUnits(stringify(geometry), diameter, well_depth,
            dipole, polarizability, rotational_relaxation, acentric_factor,
            dispersion_coefficient, quadrupole_polarizability)

    @property
    def input_data(self) -> _GasTransportInput:
        """
        Get input data defining this GasTransportData object, along with any
        user-specified data provided with its input (YAML) definition.
        """
        return anymap_to_py(self.data.parameters(True))

    def update_user_data(self, data: _GasTransportInput) -> None:
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.data.input.update(py_to_anymap(data), False)

    def clear_user_data(self) -> None:
        """
        Clear all saved input data, so that the data given by `input_data` or
        `Solution.write_yaml` will only include values generated by Cantera based on
        the current object state.
        """
        self.data.input.clear()

    @property
    def geometry(self) -> _GeometryOptions:
        """
        Get/Set the string specifying the molecular geometry. One of `atom`,
        `linear`, or `nonlinear`.
        """
        return pystr(self.data.geometry)

    @geometry.setter
    def geometry(self, geometry: _GeometryOptions) -> None:
        self.data.geometry = stringify(geometry)

    @property
    def diameter(self) -> float:
        """ Get/Set the Lennard-Jones collision diameter [m] """
        return self.data.diameter

    @diameter.setter
    def diameter(self, diameter: float) -> None:
        self.data.diameter = diameter

    @property
    def well_depth(self) -> float:
        """ Get/Set the Lennard-Jones well depth [J] """
        return self.data.well_depth

    @well_depth.setter
    def well_depth(self, well_depth: float) -> None:
        self.data.well_depth = well_depth

    @property
    def dipole(self) -> float:
        """ Get/Set the permanent dipole moment of the molecule [Coulomb-m]. """
        return self.data.dipole

    @dipole.setter
    def dipole(self, dipole: float) -> None:
        self.data.dipole = dipole

    @property
    def polarizability(self) -> float:
        """Get/Set the polarizability of the molecule [m³]."""
        return self.data.polarizability

    @polarizability.setter
    def polarizability(self, polarizability: float) -> None:
        self.data.polarizability = polarizability

    @property
    def rotational_relaxation(self) -> float:
        """
        Get/Set the rotational relaxation number (the number of collisions it
        takes to equilibrate the rotational degrees of freedom with the
        temperature).
        """
        return self.data.rotational_relaxation

    @rotational_relaxation.setter
    def rotational_relaxation(self, rotational_relaxation: float) -> None:
        self.data.rotational_relaxation = rotational_relaxation

    @property
    def acentric_factor(self) -> float:
        """Get/Set Pitzer's acentric factor [dimensionless]."""
        return self.data.acentric_factor

    @acentric_factor.setter
    def acentric_factor(self, acentric_factor: float) -> None:
        self.data.acentric_factor = acentric_factor

    @property
    def dispersion_coefficient(self) -> float:
        """Get/Set dispersion coefficient [m⁵]."""
        return self.data.dispersion_coefficient

    @dispersion_coefficient.setter
    def dispersion_coefficient(self, dispersion_coefficient: float) -> None:
        self.data.dispersion_coefficient = dispersion_coefficient

    @property
    def quadrupole_polarizability(self) -> float:
        """Get/Set quadrupole polarizability [m⁵]."""
        return self.data.quadrupole_polarizability

    @quadrupole_polarizability.setter
    def quadrupole_polarizability(self, quadrupole_polarizability: float) -> None:
        self.data.quadrupole_polarizability = quadrupole_polarizability


@cython.cclass
class Transport(_SolutionBase):
    """
    This class is used to compute transport properties for a phase of matter.

    Not all transport properties are implemented in all transport models.
    """
    def __init__(
        self,
        infile: _Path | str = "",
        name: str = "",
        adjacent: _Sequence[_ThermoPhase] = (),
        *,
        origin: _SolutionBase | None = None,
        yaml: str | None = None,
        thermo: str | None = None,
        species: _Sequence[_Species] | None = (),
        kinetics: str | None = None,
        reactions: _Sequence[_Reaction] | None = (),
        init: bool = True,
        **kwargs: _Any,
    ) -> None:
        super().__init__(infile, name, adjacent, origin=origin, yaml=yaml,
                         thermo=thermo, species=species, kinetics=kinetics,
                         reactions=reactions, init=init, **kwargs)
        if self._references is None:
            raise ValueError(
                "Cannot instantiate stand-alone 'Transport' object as it requires an "
                "associated thermo phase.\nAll 'Transport' methods should be accessed "
                "from a 'Solution' object.")

    @property
    def transport_model(self) -> _TransportModel:
        """
        Get/Set the string specifying the transport model, such as `multicomponent`,
        `mixture-averaged`, or `unity-Lewis-number`.

        Setting a new transport model deletes the underlying C++ Transport
        object and replaces it with a new one implementing the specified model.
        """
        return pystr(self.transport.transportModel())

    @transport_model.setter
    def transport_model(self, model: _TransportModel) -> None:
        self.base.setTransportModel(stringify(model))

    @property
    def CK_mode(self) -> bool:
        """Boolean to indicate if the chemkin interpretation is used."""
        return self.transport.CKMode()

    @property
    def viscosity(self) -> float:
        """Viscosity [Pa-s]."""
        return self.transport.viscosity()

    @property
    def species_viscosities(self) -> _Array:
        """Pure species viscosities [Pa-s]"""
        return get_transport_1d(self, tran_getSpeciesViscosities)

    @property
    def electrical_conductivity(self) -> float:
        """Electrical conductivity. [S/m]."""
        return self.transport.electricalConductivity()

    @property
    def thermal_conductivity(self) -> float:
        """
        Thermal conductivity. [W/m/K]
        """
        return self.transport.thermalConductivity()

    @property
    def mix_diff_coeffs(self) -> _Array:
        """
        Mixture-averaged diffusion coefficients [m²/s] relating the
        mass-averaged diffusive fluxes (with respect to the mass averaged
        velocity) to gradients in the species mole fractions.
        """
        return get_transport_1d(self, tran_getMixDiffCoeffs)

    @property
    def mix_diff_coeffs_mass(self) -> _Array:
        """
        Mixture-averaged diffusion coefficients [m²/s] relating the
        diffusive mass fluxes to gradients in the species mass fractions.
        """
        return get_transport_1d(self, tran_getMixDiffCoeffsMass)

    @property
    def mix_diff_coeffs_mole(self) -> _Array:
        """
        Mixture-averaged diffusion coefficients [m²/s] relating the
        molar diffusive fluxes to gradients in the species mole fractions.
        """
        return get_transport_1d(self, tran_getMixDiffCoeffsMole)

    @property
    def thermal_diff_coeffs(self) -> _Array:
        """
        Return a one-dimensional array of the species thermal diffusion
        coefficients [kg/m/s].
        """
        return get_transport_1d(self, tran_getThermalDiffCoeffs)

    @property
    def multi_diff_coeffs(self) -> _Array:
        """Multicomponent diffusion coefficients, D[i,j], the diffusion
        coefficient for species i due to concentration gradients in
        species j [m**2/s]."""
        return get_transport_2d(self, tran_getMultiDiffCoeffs)

    @property
    def binary_diff_coeffs(self) -> _Array:
        """Binary diffusion coefficients [m²/s]."""
        return get_transport_2d(self, tran_getBinaryDiffCoeffs)

    @property
    def mobilities(self) -> _Array:
        """
        Electrical mobilities of charged species [m²/s/V].
        """
        return get_transport_1d(self, tran_getMobilities)

    def get_viscosity_polynomial(self, i: int) -> _Array:
        """
        Get the coefficients of the polynomial fit for the viscosity of species ``i``
        as a function of temperature. See :ct:`GasTransport::fitProperties` for the
        functional form.
        """
        n_values = 4 if self.transport.CKMode() else 5
        return get_transport_polynomial(self, tran_getViscosityPolynomial, i, n_values)

    def get_thermal_conductivity_polynomial(self, i: int) -> _Array:
        """
        Get the coefficients of the polynomial fit for the thermal conductivity of
        species ``i`` as a function of temperature. See
        :ct:`GasTransport::fitProperties` for the functional form.
        """
        n_values = 4 if self.transport.CKMode() else 5
        return get_transport_polynomial(self, tran_getConductivityPolynomial, i,
                                        n_values)

    def get_binary_diff_coeffs_polynomial(self, i: int, j: int) -> _Array:
        """
        Get the coefficients of the temperature-dependent polynomial used as part of the
        fit for the binary diffusion coefficient for species ``i`` and ``j``. See
        :ct:`GasTransport::fitDiffCoeffs` for the functional form.
        """
        n_values = 4 if self.transport.CKMode() else 5
        return get_binary_transport_polynomial(self, tran_getBinDiffusivityPolynomial,
                                               i, j, n_values)

    def get_collision_integral_polynomials(
        self, i: int, j: int
    ) -> tuple[_Array, _Array, _Array]:
        """
        Get the coefficients of the polynomial fit of the collision integral for
        species ``i`` and ``j`` as a function of the reduced temperature. See
        :ct:`GasTransport::fitCollisionIntegrals`.
        """
        n_values = 7 if self.transport.CKMode() else 9
        adata = np.empty(n_values)
        bdata = np.empty(n_values)
        cdata_arr = np.empty(n_values)
        acdata: cython.double[::1] = adata
        bcdata: cython.double[::1] = bdata
        ccdata: cython.double[::1] = cdata_arr
        self.transport.getCollisionIntegralPolynomial(i, j,
            span[cython.double](cython.address(acdata[0]),
                                cython.cast(cython.size_t, adata.size)),
            span[cython.double](cython.address(bcdata[0]),
                                cython.cast(cython.size_t, bdata.size)),
            span[cython.double](cython.address(ccdata[0]),
                                cython.cast(cython.size_t, cdata_arr.size)))
        return adata, bdata, cdata_arr

    def set_viscosity_polynomial(self, i: int, values: _ArrayLike) -> None:
        """
        Set the coefficients of the polynomial fit for the viscosity of species ``i``
        as a function of temperature. See :ct:`GasTransport::fitProperties` for the
        functional form.
        """
        n_values = 4 if self.transport.CKMode() else 5
        data = np.ascontiguousarray(values, dtype=np.double)
        if len(data) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(data)}.""")
        cdata: cython.double[::1] = data
        tran_setViscosityPolynomial(self.transport, i,
                                    span[cython.double](cython.address(cdata[0]),
                                                        cython.cast(cython.size_t, data.size)))

    def set_thermal_conductivity_polynomial(self, i: int, values: _ArrayLike) -> None:
        """
        Set the coefficients of the polynomial fit for the thermal conductivity of
        species ``i`` as a function of temperature. See
        :ct:`GasTransport::fitProperties` for the functional form.
        """
        n_values = 4 if self.transport.CKMode() else 5
        data = np.ascontiguousarray(values, dtype=np.double)
        if len(data) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(data)}.""")
        cdata: cython.double[::1] = data
        tran_setConductivityPolynomial(self.transport, i,
                                       span[cython.double](cython.address(cdata[0]),
                                                           cython.cast(cython.size_t, data.size)))

    def set_binary_diff_coeffs_polynomial(
        self, i: int, j: int, values: _ArrayLike
    ) -> None:
        """
        Set the coefficients of the polynomial fit of the binary diffusion coefficient
        for species ``i`` and ``j`` as a function of temperature. See
        :ct:`GasTransport::fitDiffCoeffs` for the functional form.
        """
        n_values = 4 if self.transport.CKMode() else 5
        data = np.ascontiguousarray(values, dtype=np.double)
        if len(data) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(data)}.""")
        cdata: cython.double[::1] = data
        tran_setBinDiffusivityPolynomial(self.transport, i, j,
                                         span[cython.double](cython.address(cdata[0]),
                                                             cython.cast(cython.size_t, data.size)))

    def set_collision_integral_polynomial(
        self, i: int, j: int, avalues: _ArrayLike, bvalues: _ArrayLike,
        cvalues: _ArrayLike, actualT: bool = False
    ) -> None:
        r"""
        Set the coefficients of the polynomial fit of the collision integral for
        species ``i`` and ``j`` as a function of temperature. See
        :ct:`GasTransport::fitCollisionIntegrals`. The flag ``actualT`` determines
        whether the fit is done in terms of the reduced temperature
        (:math:`T^*_{ij} = T k_B / \epsilon_{ij}`; default) or the actual temperature.
        """
        n_values = 7 if self.transport.CKMode() else 9
        adata = np.ascontiguousarray(avalues, dtype=np.double)
        bdata = np.ascontiguousarray(bvalues, dtype=np.double)
        cdata_arr = np.ascontiguousarray(cvalues, dtype=np.double)
        if len(adata) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(adata)}.""")
        if len(bdata) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(bdata)}.""")
        if len(cdata_arr) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(cdata_arr)}.""")
        acdata: cython.double[::1] = adata
        bcdata: cython.double[::1] = bdata
        ccdata: cython.double[::1] = cdata_arr
        self.transport.setCollisionIntegralPolynomial(i, j,
            span[cython.double](cython.address(acdata[0]),
                                cython.cast(cython.size_t, adata.size)),
            span[cython.double](cython.address(bcdata[0]),
                                cython.cast(cython.size_t, bdata.size)),
            span[cython.double](cython.address(ccdata[0]),
                                cython.cast(cython.size_t, cdata_arr.size)), actualT)

    @property
    def transport_fitting_errors(self) -> _TransportFittingErrors:
        """
        Get error metrics about any functional fits calculated for pure species
        transport properties. See :ct:`GasTransport::fitDiffCoeffs` and
        :ct:`GasTransport::fitProperties`.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.

        .. versionadded:: 3.1
        """
        stats: CxxAnyMap = self.transport.fittingErrors()
        return anymap_to_py(stats)


@cython.cclass
class DustyGasTransport(Transport):
    """
    Implements the "dusty gas" model for transport in porous media.

    As implemented here, only species transport (`~Transport.multi_diff_coeffs`)
    is handled. The viscosity, thermal conductivity, and thermal diffusion
    coefficients are not implemented.
    """
    def __init__(
        self,
        infile: _Path | str = "",
        name: str = "",
        adjacent: _Sequence[_ThermoPhase] = (),
        *,
        origin: _SolutionBase | None = None,
        yaml: str | None = None,
        thermo: str | None = None,
        species: _Sequence[_Species] | None = (),
        kinetics: str | None = None,
        reactions: _Sequence[_Reaction] | None = (),
        init: bool = True,
        **kwargs: _Any,
    ) -> None:
        self.base.setTransportModel(stringify("DustyGas"))
        self.transport = self.base.transport().get()
        super().__init__(infile, name, adjacent, origin=origin, yaml=yaml,
                         thermo=thermo, species=species, kinetics=kinetics,
                         reactions=reactions, init=init, **kwargs)

    @property
    def porosity(self) -> _Never:
        """Porosity of the porous medium [dimensionless]."""
        raise AttributeError("unreadable attribute 'porosity'")

    @porosity.setter
    def porosity(self, value: float) -> None:
        cython.cast(cython.pointer(CxxDustyGasTransport), self.transport).setPorosity(value)

    @property
    def tortuosity(self) -> _Never:
        """Tortuosity of the porous medium [dimensionless]."""
        raise AttributeError("unreadable attribute 'tortuosity'")

    @tortuosity.setter
    def tortuosity(self, value: float) -> None:
        cython.cast(cython.pointer(CxxDustyGasTransport), self.transport).setTortuosity(value)

    @property
    def mean_pore_radius(self) -> _Never:
        """Mean pore radius of the porous medium [m]."""
        raise AttributeError("unreadable attribute 'mean_pore_radius'")

    @mean_pore_radius.setter
    def mean_pore_radius(self, value: float) -> None:
        cython.cast(cython.pointer(CxxDustyGasTransport), self.transport).setMeanPoreRadius(value)

    @property
    def mean_particle_diameter(self) -> _Never:
        """Mean particle diameter of the porous medium [m]."""
        raise AttributeError("unreadable attribute 'mean_particle_diameter'")

    @mean_particle_diameter.setter
    def mean_particle_diameter(self, value: float) -> None:
        cython.cast(cython.pointer(CxxDustyGasTransport), self.transport).setMeanParticleDiameter(value)

    @property
    def permeability(self) -> _Never:
        """Permeability of the porous medium [m²]."""
        raise AttributeError("unreadable attribute 'permeability'")

    @permeability.setter
    def permeability(self, value: float) -> None:
        cython.cast(cython.pointer(CxxDustyGasTransport), self.transport).setPermeability(value)

    @property
    def thermal_conductivity(self) -> float:
        """
        Thermal conductivity. [W/m/K]
        Returns the thermal conductivity of the ideal gas object using the
        multicomponent model. The value is not specific to the dusty gas model.
        """
        return cython.cast(cython.pointer(CxxDustyGasTransport),
                           self.transport).gasTransport().thermalConductivity()

    def molar_fluxes(
        self, T1: float, T2: float, rho1: float, rho2: float,
        Y1: _ArrayLike, Y2: _ArrayLike, delta: float
    ) -> _Array:
        """
        Get the molar fluxes [kmol/m²/s], given the thermodynamic state at
        two nearby points.

        :param T1:
            Temperature [K] at the first point
        :param T2:
            Temperature [K] at the second point
        :param rho1:
            Density [kg/m³] at the first point
        :param rho2:
            Density [kg/m³] at the second point
        :param Y1:
            Array of mass fractions at the first point. Length `n_species`.
        :param Y2:
            Array of mass fractions at the second point. Length `n_species`.
        :param delta:
            Distance [m] between the two points.
        """

        state1 = np.empty(self.n_species + 2)
        state2 = np.empty(self.n_species + 2)
        fluxes = np.empty(self.n_species)

        state1[0] = T1
        state1[1] = rho1
        state1[2:] = Y1
        state2[0] = T2
        state2[1] = rho2
        state2[2:] = Y2

        cstate1: cython.double[::1] = state1
        cstate2: cython.double[::1] = state2
        cfluxes: cython.double[::1] = fluxes

        cython.cast(cython.pointer(CxxDustyGasTransport), self.transport).getMolarFluxes(
            span[const_double](cython.address(cstate1[0]),
                               cython.cast(cython.size_t, state1.size)),
            span[const_double](cython.address(cstate2[0]),
                               cython.cast(cython.size_t, state2.size)),
            delta,
            span[cython.double](cython.address(cfluxes[0]),
                                cython.cast(cython.size_t, fluxes.size)))
        return fluxes
