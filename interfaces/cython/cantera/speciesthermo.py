# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3

from typing import Any as _Any, ClassVar as _ClassVar, TypedDict as _TypedDict
from typing_extensions import Required as _Required

import numpy as np

import cython
from cython.cimports.cantera._utils import anymap_to_py, py_to_anymap
from .constants import gas_constant
from ._types import Array as _Array, ArrayLike as _ArrayLike

# These match the definitions in speciesThermoTypes.h
SPECIES_THERMO_CONSTANT_CP = cython.declare(cython.int, 1)
SPECIES_THERMO_NASA2 = cython.declare(cython.int, 4)
SPECIES_THERMO_SHOMATE2 = cython.declare(cython.int, 8)
SPECIES_THERMO_NASA9MULTITEMP = cython.declare(cython.int, 513)
SPECIES_THERMO_MU0_INTERP = cython.declare(cython.int, 64)

_SpeciesThermoInput = _TypedDict(
    "_SpeciesThermoInput",
    {
        "model": _Required[str],
        "temperature-ranges": list[float],
        "data": list[list[float]],
        "note": str,
    },
    total=False,
)


@cython.cclass
class SpeciesThermo:
    """
    Base class for representing the reference-state thermodynamic properties of
    a pure species. These properties are a function of temperature. Derived
    classes implement a parameterization of this temperature dependence. This is
    a wrapper for the C++ class :ct:`SpeciesThermoInterpType`.

    :param T_low:
        The minimum temperature [K] at which the parameterization is valid
    :param T_high:
        The maximum temperature [K] at which the parameterization is valid
    :param P_ref:
        The reference pressure [Pa] for the parameterization
    :param coeffs:
        An array of coefficients for the parameterization. The length of this
        array and the meaning of each element depends on the specific
        parameterization.
    """
    derived_type: _ClassVar[int]

    def __cinit__(self, T_low=None, T_high=None, P_ref=None, coeffs=None, *args,
                  init=True, **kwargs):
        if not init:
            return

        if not self._check_n_coeffs(len(coeffs)):
            raise ValueError("Coefficient array has incorrect length")
        data = np.ascontiguousarray(coeffs, dtype=np.double)
        cdata: cython.double[::1] = data
        view: span[cython.double] = span[cython.double](
            cython.address(cdata[0]), cython.cast(cython.size_t, data.size))
        self._spthermo.reset(CxxNewSpeciesThermo(self.derived_type, T_low,
                                                 T_high, P_ref, view))
        self.spthermo = self._spthermo.get()

    def __init__(self, T_low: float | None = None, T_high: float | None = None,
                 P_ref: float | None = None, coeffs: _ArrayLike | None = None,
                 *args: _Any, init: bool = True, **kwargs: _Any) -> None:
        # The C++ object is constructed in __cinit__; this typed __init__ exists so
        # that mypy/pyright (which do not recognize Cython's __cinit__) publish the
        # constructor signature. Constructor parameters are documented in the class
        # docstring above.
        pass

    @cython.cfunc
    def _assign(self, other: shared_ptr[CxxSpeciesThermo]):
        self._spthermo = other
        self.spthermo = self._spthermo.get()

    @property
    def min_temp(self) -> float:
        """ Minimum temperature [K] at which the parameterization is valid."""
        return self.spthermo.minTemp()

    @property
    def max_temp(self) -> float:
        """ Maximum temperature [K] at which the parameterization is valid."""
        return self.spthermo.maxTemp()

    @property
    def reference_pressure(self) -> float:
        """ Reference pressure [Pa] for the parameterization."""
        return self.spthermo.refPressure()

    @property
    def n_coeffs(self) -> int:
        """ Number of parameters for the parameterization."""
        return self.spthermo.nCoeffs()

    @property
    def coeffs(self) -> _Array:
        """
        Array of coefficients for the parameterization. The length of this
        array and the meaning of each element depends on the specific
        parameterization.
        """
        index: cython.size_t = 0
        thermo_type: cython.int = 0
        T_low: cython.double = 0
        T_high: cython.double = 0
        P_ref: cython.double = 0
        data = np.empty(self.n_coeffs)
        cdata: cython.double[::1] = data
        view: span[cython.double] = span[cython.double](
            cython.address(cdata[0]), cython.cast(cython.size_t, self.n_coeffs))
        self.spthermo.reportParameters(index, thermo_type, T_low,
                                       T_high, P_ref, view)
        return data

    def _check_n_coeffs(self, n: int) -> bool:
        """
        Check whether number of coefficients is compatible with a given
        parameterization prior to instantiation of the underlying C++ object.
        """
        raise NotImplementedError('Needs to be overloaded')

    @property
    def input_data(self) -> _SpeciesThermoInput:
        """
        Get input data defining this SpeciesThermo object, along with any user-specified
        data provided with its input (YAML) definition.
        """
        return anymap_to_py(self.spthermo.parameters(True))

    def update_user_data(self, data: dict[str, _Any]) -> None:
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.spthermo.input().update(py_to_anymap(data), False)

    def clear_user_data(self) -> None:
        """
        Clear all saved input data, so that the data given by `input_data` or
        `Solution.write_yaml` will only include values generated by Cantera based on
        the current object state.
        """
        self.spthermo.input().clear()

    def cp(self, T: float) -> float:
        """
        Molar heat capacity at constant pressure [J/kmol/K] at temperature *T*.
        """
        cp_r: cython.double = 0
        h_rt: cython.double = 0
        s_r: cython.double = 0
        self.spthermo.updatePropertiesTemp(T, cp_r, h_rt, s_r)
        return cp_r * gas_constant

    def h(self, T: float) -> float:
        """ Molar enthalpy [J/kmol] at temperature *T* """
        cp_r: cython.double = 0
        h_rt: cython.double = 0
        s_r: cython.double = 0
        self.spthermo.updatePropertiesTemp(T, cp_r, h_rt, s_r)
        return h_rt * gas_constant * T

    def s(self, T: float) -> float:
        """ Molar entropy [J/kmol/K] at temperature *T* """
        cp_r: cython.double = 0
        h_rt: cython.double = 0
        s_r: cython.double = 0
        self.spthermo.updatePropertiesTemp(T, cp_r, h_rt, s_r)
        return s_r * gas_constant


@cython.cclass
class ConstantCp(SpeciesThermo):
    r"""
    Thermodynamic properties for a species that has a constant specific heat
    capacity. This is a wrapper for the C++ class :ct:`ConstCpPoly`.

    :param coeffs:
        An array of 4 elements:

        - `coeffs[0]` = :math:`T_0` [K]
        - `coeffs[1]` = :math:`H^o(T_0, p_{ref})` [J/kmol]
        - `coeffs[2]` = :math:`S^o(T_0, p_{ref})` [J/kmol-K]
        - `coeffs[3]` = :math:`c_p^o(T_0, p_{ref})` [J/kmol-K]
    """
    derived_type = SPECIES_THERMO_CONSTANT_CP

    def _check_n_coeffs(self, n: int) -> bool:
        return n == 4


@cython.cclass
class Mu0Poly(SpeciesThermo):
    """
    Thermodynamic properties for a species which is parameterized using an
    interpolation of the Gibbs free energy based on a piecewise constant heat
    capacity approximation. This is a wrapper for the C++ class :ct:`Mu0Poly`.

    :param coeffs:
        An array of `2 + 2*npoints` elements, in the following order:

        - `coeffs[0]`: number of points (integer)
        - `coeffs[1]`: h^o(298.15 K) [J/kmol]
        - `coeffs[2]`: T_1 [Kelvin]
        - `coeffs[3]`: \mu^o(T_1) [J/kmol]
        - `coeffs[4]`: T_2 [Kelvin]
        - `coeffs[5]`: \mu^o(T_2) [J/kmol]
        - ...
    """
    derived_type = SPECIES_THERMO_MU0_INTERP

    def _check_n_coeffs(self, n: int) -> bool:
        return n > 3 and n % 2 == 0


@cython.cclass
class NasaPoly2(SpeciesThermo):
    """
    Thermodynamic properties for a species which is parameterized using the
    7-coefficient NASA polynomial form in two temperature ranges. This is a
    wrapper for the C++ class :ct:`NasaPoly2`.

    :param coeffs:
        An array of 15 elements, in the following order:

        - `coeffs[0]`: The mid-point temperature [K] between the two
          parameterizations
        - `coeffs[1:8]`: The 7 coefficients of the high-temperature
          parameterization
        - `coeffs[8:15]`: The 7 coefficients of the low-temperature
          parameterization

        This is the coefficient order used in the standard fixed-format NASA
        input files.
    """
    derived_type = SPECIES_THERMO_NASA2

    def _check_n_coeffs(self, n: int) -> bool:
        return n == 15


@cython.cclass
class Nasa9PolyMultiTempRegion(SpeciesThermo):
    """
    Thermodynamic properties for a species which is parameterized using the
    9-coefficient NASA polynomial form encompassing multiple temperature ranges.
    This is a wrapper for the C++ class :ct:`Nasa9PolyMultiTempRegion`.

    :param coeffs:
        An array of `1 + 11*nzones` elements, in the following order:

        - `coeffs[0]`: Number of zones (`nzones`)
        - `coeffs[1 + 11*zone]`: minimum temperature within zone
        - `coeffs[2 + 11*zone]`: maximum temperature within zone
        - `coeffs[3:11 + 11*zone]`: 9 coefficients of the parameterization

        where `zone` runs from zero to `nzones-1`.
    """
    derived_type = SPECIES_THERMO_NASA9MULTITEMP

    def _check_n_coeffs(self, n: int) -> bool:
        return n > 11 and ((n - 1) % 11) == 0


@cython.cclass
class ShomatePoly2(SpeciesThermo):
    """
    Thermodynamic properties for a species which is parameterized using the
    Shomate equation in two temperature ranges. This is a wrapper for the C++
    class :ct:`ShomatePoly2`.

    :param coeffs:
        An array of 15 elements, in the following order:

        - `coeffs[0]`: The mid-point temperature [K] between the two
          parameterizations
        - `coeffs[1:8]`: The 7 coefficients of the low-temperature
          parameterization
        - `coeffs[8:15]`: The 7 coefficients of the high-temperature
          parameterization

        These coefficients should be provided in their customary units (that is,
        such that :math:`c_p^o` is in J/gmol-K and :math:`H^o` is in kJ/gmol,
        as in the NIST Chemistry WebBook).
    """
    derived_type = SPECIES_THERMO_SHOMATE2

    def _check_n_coeffs(self, n: int) -> bool:
        return n == 15


@cython.cfunc
def wrapSpeciesThermo(spthermo: shared_ptr[CxxSpeciesThermo]):
    """
    Wrap a C++ SpeciesThermoInterpType object with a Python object of the
    correct derived type.
    """
    thermo_type: cython.int = spthermo.get().reportType()

    if thermo_type == SPECIES_THERMO_NASA2:
        st = NasaPoly2(init=False)
    elif thermo_type == SPECIES_THERMO_NASA9MULTITEMP:
        st = Nasa9PolyMultiTempRegion(init=False)
    elif thermo_type == SPECIES_THERMO_CONSTANT_CP:
        st = ConstantCp(init=False)
    elif thermo_type == SPECIES_THERMO_MU0_INTERP:
        st = Mu0Poly(init=False)
    elif thermo_type == SPECIES_THERMO_SHOMATE2:
        st = ShomatePoly2(init=False)
    else:
        return None

    st._assign(spthermo)
    return st
