# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3
# pyright: reportMissingImports=false, reportAttributeAccessIssue=false
# pyright: reportUndefinedVariable=false, reportUnboundVariable=false
# pyright: reportInvalidTypeArguments=false, reportAssignmentType=false
# pyright: reportIndexIssue=false, reportInvalidTypeForm=false

from ctypes import c_int as _c_int
import warnings
from collections.abc import Sequence as _Sequence
from pathlib import Path as _Path
from typing import Any as _Any, TypeAlias as _TypeAlias, TypedDict as _TypedDict, \
    TYPE_CHECKING
import numpy as np

import cython
import cython.cimports.numpy as cnp  # Required: triggers import_array() for PyArray_SIZE
from cython.cimports.cantera.reaction import CustomRate, ExtensibleRate, Reaction
from cython.cimports.cantera._utils import stringify, pystr, anymap_to_py, py_to_anymap
from .solutionbase import _SolutionBase
from ._types import Array as _Array
from . import _utils

if TYPE_CHECKING:
    from .reaction import CustomRate as _CustomRate, Reaction as _Reaction
    from .thermo import Species as _Species, ThermoPhase as _ThermoPhase
    from .units import UnitSystem as _UnitSystem, _UnitDict as _UnitDict

# Avoid fixed options unless we can find a way to support custom extensions:
# _KineticsType: TypeAlias = Literal["none", "bulk", "edge", "surface"]
_KineticsType: _TypeAlias = str

_DerivativeSettings = _TypedDict(
    "_DerivativeSettings",
    {
        "skip-third-bodies": bool,
        "skip-falloff": bool,
        "rtol-delta": float,
        "skip-coverage-dependence": bool,
        "skip-electrochemistry": bool,
        "skip-nonideal": bool,
        "skip-flow-devices": bool,
        "skip-walls": bool,
        "skip-connector-composition-dependence": bool,
        "skip-connector-pressure-composition-dependence": bool,
    },
    total=False,
)

# NOTE: These cdef functions cannot be members of Kinetics because they would
# cause "layout conflicts" when creating derived classes with multiple bases,
# such as class Solution. [Cython 0.16]
@cython.cfunc
def get_species_array(kin: Kinetics, method: kineticsMethod1d) -> np.ndarray:
    data = np.empty(kin.n_total_species)
    if kin.n_total_species == 0:
        return data
    cdata: cython.double[::1] = data
    view: span[cython.double] = span[cython.double](
        cython.address(cdata[0]), cython.cast(cython.size_t, data.size))
    method(kin.kinetics, view)
    # @todo: Fix _selected_species to work with interface kinetics
    if kin._selected_species.size:
        return data[kin._selected_species]
    else:
        return data


@cython.cfunc
def get_reaction_array(kin: Kinetics, method: kineticsMethod1d) -> np.ndarray:
    data = np.empty(kin.n_reactions)
    if kin.n_reactions == 0:
        return data
    cdata: cython.double[::1] = data
    view: span[cython.double] = span[cython.double](
        cython.address(cdata[0]), cython.cast(cython.size_t, data.size))
    method(kin.kinetics, view)
    return data


@cython.cfunc
def get_dense(smat: CxxSparseMatrix) -> np.ndarray:
    length: cython.size_t = smat.nonZeros()
    if length == 0:
        return np.zeros((smat.rows(), 0))

    # index/value triplets
    rows = np.empty(length, dtype=_c_int)
    cols = np.empty(length, dtype=_c_int)
    data = np.empty(length)

    crows: cython.int[::1] = rows
    ccols: cython.int[::1] = cols
    cdata: cython.double[::1] = data

    size = CxxSparseTriplets(smat, cython.address(crows[0]), cython.address(ccols[0]),
                             cython.address(cdata[0]), length)
    out = np.zeros((smat.rows(), smat.cols()))
    for i in range(length):
        out[rows[i], cols[i]] = data[i]
    return out


@cython.cfunc
def get_sparse(smat: CxxSparseMatrix):
    # pointers to values and inner indices of CSC storage
    length: cython.size_t = smat.nonZeros()
    value = np.empty(length)
    inner = np.empty(length, dtype=_c_int)

    # pointers outer indices of CSC storage
    ncols: cython.size_t = smat.outerSize()
    outer = np.empty(ncols + 1, dtype=_c_int)

    cvalue: cython.double[::1] = value
    cinner: cython.int[::1] = inner
    couter: cython.int[::1] = outer

    CxxSparseCscData(smat, cython.address(cvalue[0]), cython.address(cinner[0]),
                     cython.address(couter[0]))
    return value, inner, outer


@cython.cfunc
def get_from_sparse(smat: CxxSparseMatrix, rows: cython.int, cols: cython.int):
    if _utils._USE_SPARSE:
        tup = get_sparse(smat)
        return _utils._scipy_sparse.csc_matrix(tup, shape=(rows, cols))
    else:
        return get_dense(smat)


@cython.cclass
class Kinetics(_SolutionBase):
    """
    Instances of class `Kinetics` are responsible for evaluating reaction rates
    of progress, species production rates, and other quantities pertaining to
    a reaction mechanism.
    """
    #: List holding references to CustomRate objects (prevents garbage collection)
    _custom_rates: list[_CustomRate] = []

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
                "Cannot instantiate stand-alone 'Kinetics' object as it requires an "
                "associated thermo phase.\nAll 'Kinetics' methods should be accessed "
                "from a 'Solution' object.")

    @property
    def kinetics_model(self) -> str:
        """
        Return type of kinetics.
        """
        return pystr(self.kinetics.kineticsType())

    @property
    def n_total_species(self) -> int:
        """
        Total number of species in all phases participating in the kinetics
        mechanism.
        """
        return self.kinetics.nTotalSpecies()

    @property
    def n_reactions(self) -> int:
        """Number of reactions in the reaction mechanism."""
        return self.kinetics.nReactions()

    @property
    def n_phases(self) -> int:
        """Number of phases in the reaction mechanism."""
        return self.kinetics.nPhases()

    def _check_phase_index(self, n):
        if not 0 <= n < self.n_phases:
            raise ValueError("Phase index ({0}) out of range".format(n))

    def _check_reaction_index(self, n):
        if not 0 <= n < self.n_reactions:
            raise ValueError("Reaction index ({0}) out of range".format(n))

    def _check_kinetics_species_index(self, n):
        if not 0 <= n < self.n_total_species:
            raise ValueError("Kinetics Species index ({0}) out of range".format(n))

    def kinetics_species_index(self, species: str | int, phase: int = 0) -> int:
        """
        The index of species ``species`` of phase ``phase`` within arrays returned
        by methods of class `Kinetics`. If ``species`` is a string, the ``phase``
        argument is unused.
        """
        k: cython.int
        if isinstance(species, (str, bytes)):
            return self.kinetics.kineticsSpeciesIndex(stringify(species))
        else:
            k = species
            self._check_kinetics_species_index(k)
            self._check_phase_index(k)
            return self.kinetics.kineticsSpeciesIndex(k, phase)

    def kinetics_species_name(self, k: int) -> str:
        """
        Name of the species with index ``k`` in the arrays returned by methods
        of class `Kinetics`.
        """
        return pystr(self.kinetics.kineticsSpeciesName(k))

    @property
    def kinetics_species_names(self) -> list[str]:
        """
        A list of all species names, corresponding to the arrays returned by
        methods of class `Kinetics`.
        """
        return [self.kinetics_species_name(k)
                for k in range(self.n_total_species)]

    def reaction(self, i_reaction: int) -> _Reaction:
        """
        Return a `Reaction` object representing the reaction with index
        ``i_reaction``. Changes to this object do not affect the `Kinetics` or
        `Solution` object until the `modify_reaction` function is called.
        """
        self._check_reaction_index(i_reaction)
        return Reaction.wrap(self.kinetics.reaction(i_reaction))

    def reactions(self) -> list[_Reaction]:
        """
        Return a list of all `Reaction` objects. Changes to these objects do not
        affect the `Kinetics` or `Solution` object until the `modify_reaction`
        function is called.
        """
        return [self.reaction(i) for i in range(self.n_reactions)]

    def modify_reaction(self, irxn: int, rxn: _Reaction) -> None:
        """
        Modify the `Reaction` with index ``irxn`` to have the same rate parameters as
        ``rxn``. ``rxn`` must have the same reactants and products and use the same rate
        parameterization (for example, `ArrheniusRate`, `FalloffRate`, `PlogRate`, etc.)
        as the existing reaction. This method does not modify the third-body
        efficiencies, reaction orders, or reversibility of the reaction.
        """
        self.kinetics.modifyReaction(irxn, cython.cast(Reaction, rxn)._reaction)

    def add_reaction(self, rxn: _Reaction) -> None:
        """ Add a new reaction to this phase. """
        self.kinetics.addReaction(cython.cast(Reaction, rxn)._reaction)
        if isinstance(rxn.rate, (CustomRate, ExtensibleRate)):
            # prevent premature garbage collection
            self._custom_rates.append(rxn.rate)

    def multiplier(self, i_reaction: int) -> float:
        """
        A scaling factor applied to the rate coefficient for reaction
        ``i_reaction``. Can be used to carry out sensitivity analysis or to
        selectively disable a particular reaction. See `set_multiplier`.
        """
        self._check_reaction_index(i_reaction)
        return self.kinetics.multiplier(i_reaction)

    def set_multiplier(self, value: float, i_reaction: int = -1) -> None:
        """
        Set the multiplier for for reaction ``i_reaction`` to ``value``.
        If ``i_reaction`` is not specified, then the multiplier for all reactions
        is set to ``value``. See `multiplier`.
        """
        if i_reaction == -1:
            for i_reaction in range(self.n_reactions):
                self.kinetics.setMultiplier(i_reaction, value)
        else:
            self._check_reaction_index(i_reaction)
            self.kinetics.setMultiplier(i_reaction, value)

    def reaction_equations(self, indices: _Sequence[int] | None = None) -> list[str]:
        """
        Returns a list containing the reaction equation for all reactions in the
        mechanism if ``indices`` is unspecified, or the equations for each
        reaction in the sequence ``indices``. For example::

            >>> gas.reaction_equations()
            ['2 O + M <=> O2 + M', 'O + H + M <=> OH + M', 'O + H2 <=> H + OH', ...]
            >>> gas.reaction_equations([2,3])
            ['O + H + M <=> OH + M', 'O + H2 <=> H + OH']

        See also `reaction_equation`.
        """
        if indices is None:
            return list(r.equation for r in self.reactions())
        else:
            return [self.reaction(i).equation for i in indices]

    def reactant_stoich_coeff(self, k_spec: str | int, i_reaction: int) -> float:
        """
        The stoichiometric coefficient of species ``k_spec`` as a reactant in
        reaction ``i_reaction``.
        """
        k: cython.int
        if isinstance(k_spec, (str, bytes)):
            k = self.kinetics_species_index(k_spec)
        else:
            k = k_spec
            self._check_kinetics_species_index(k_spec)

        self._check_reaction_index(i_reaction)
        return self.kinetics.reactantStoichCoeff(k, i_reaction)

    def product_stoich_coeff(self, k_spec: str | int, i_reaction: int) -> float:
        """
        The stoichiometric coefficient of species ``k_spec`` as a product in
        reaction ``i_reaction``.
        """
        k: cython.int
        if isinstance(k_spec, (str, bytes)):
            k = self.kinetics_species_index(k_spec)
        else:
            k = k_spec
            self._check_kinetics_species_index(k_spec)

        self._check_reaction_index(i_reaction)
        return self.kinetics.productStoichCoeff(k, i_reaction)

    @property
    def reactant_stoich_coeffs(self) -> _Array:
        """
        The array of reactant stoichiometric coefficients. Element ``[k,i]`` of
        this array is the reactant stoichiometric coefficient of species ``k`` in
        reaction ``i``.

        For sparse output, set ``ct.use_sparse(True)``.

        .. versionchanged:: 3.0

            Method was changed to a property in Cantera 3.0.
        """
        return get_from_sparse(self.kinetics.reactantStoichCoeffs(),
                               self.n_total_species, self.n_reactions)

    @property
    def product_stoich_coeffs(self) -> _Array:
        """
        The array of product stoichiometric coefficients. Element ``[k,i]`` of
        this array is the product stoichiometric coefficient of species ``k`` in
        reaction ``i``.

        For sparse output, set ``ct.use_sparse(True)``.

        .. versionchanged:: 3.0

            Method was changed to a property in Cantera 3.0.
        """
        return get_from_sparse(self.kinetics.productStoichCoeffs(),
                               self.n_total_species, self.n_reactions)

    @property
    def product_stoich_coeffs_reversible(self) -> _Array:
        """
        The array of product stoichiometric coefficients of reversible reactions.
        Element ``[k,i]`` of this array is the product stoichiometric coefficient
        of species ``k`` in reaction ``i``.

        For sparse output, set ``ct.use_sparse(True)``.
        """
        return get_from_sparse(self.kinetics.revProductStoichCoeffs(),
                               self.n_total_species, self.n_reactions)

    @property
    def forward_rates_of_progress(self) -> _Array:
        """
        Forward rates of progress for the reactions. [kmol/m³/s] for bulk
        phases or [kmol/m²/s] for surface phases.
        """
        return get_reaction_array(self, kin_getFwdRatesOfProgress)

    @property
    def reverse_rates_of_progress(self) -> _Array:
        """
        Reverse rates of progress for the reactions. [kmol/m³/s] for bulk
        phases or [kmol/m²/s] for surface phases.
        """
        return get_reaction_array(self, kin_getRevRatesOfProgress)

    @property
    def net_rates_of_progress(self) -> _Array:
        """
        Net rates of progress for the reactions. [kmol/m³/s] for bulk phases
        or [kmol/m²/s] for surface phases.
        """
        return get_reaction_array(self, kin_getNetRatesOfProgress)

    @property
    def equilibrium_constants(self) -> _Array:
        """Equilibrium constants in concentration units for all reactions."""
        return get_reaction_array(self, kin_getEquilibriumConstants)

    @property
    def forward_rate_constants(self) -> _Array:
        """
        Forward rate constants for all reactions.

        The computed values include all temperature-dependent and pressure-dependent
        contributions. By default, third-body concentrations are only considered if
        they are part of the reaction rate definition; for a legacy implementation that
        includes third-body concentrations, see `use_legacy_rate_constants`.
        """
        return get_reaction_array(self, kin_getFwdRateConstants)

    @property
    def reverse_rate_constants(self) -> _Array:
        """
        Reverse rate constants for all reactions.

        The computed values include all temperature-dependent and pressure-dependent
        contributions. By default, third-body concentrations are only considered if
        they are part of the reaction rate definition; for a legacy implementation that
        includes third-body concentrations, see `use_legacy_rate_constants`.
        """
        return get_reaction_array(self, kin_getRevRateConstants)

    @property
    def creation_rates(self) -> _Array:
        """
        Creation rates for each species. [kmol/m³/s] for bulk phases or
        [kmol/m²/s] for surface phases.
        """
        return get_species_array(self, kin_getCreationRates)

    @property
    def destruction_rates(self) -> _Array:
        """
        Destruction rates for each species. [kmol/m³/s] for bulk phases or
        [kmol/m²/s] for surface phases.
        """
        return get_species_array(self, kin_getDestructionRates)

    @property
    def net_production_rates(self) -> _Array:
        """
        Net production rates for each species. [kmol/m³/s] for bulk phases or
        [kmol/m²/s] for surface phases.
        """
        return get_species_array(self, kin_getNetProductionRates)

    @property
    def derivative_settings(self) -> _DerivativeSettings:
        """
        Property setting behavior of derivative evaluation.

        Derivative settings are updated using a dictionary::

            >>> gas.derivative_settings = {"skip-falloff": True}

        Passing an empty dictionary will reset all values to their defaults.

        For :ct:`BulkKinetics`, the following keyword/value pairs are supported:

        -  ``skip-third-bodies`` (boolean) ... if `False` (default), third body
           concentrations are considered for the evaluation of derivatives

        -  ``skip-falloff`` (boolean) ... if `True` (default), third-body effects
           on reaction rates are not considered.

        - `skip-nonideal` (boolean): if `false` (default), derivatives are only
          supported for ideal ThermoPhase models; if `true`, derivatives for
          non-ideal phases are evaluated using idealized approximations that
          neglect non-ideal contributions.

        -  ``rtol-delta`` (double) ... relative tolerance used to perturb properties
           when calculating numerical derivatives. The default value is 1e-8.

        For :ct:`InterfaceKinetics`, the following keyword/value pairs are supported:

        - `skip-coverage-dependence` (boolean): if `false` (default), rate constant
          coverage dependence is not considered when evaluating derivatives.

        - `skip-electrochemistry` (boolean): if `false` (default), electrical charge
          is not considered in evaluating the derivatives and these reactions are
          treated as normal surface reactions.

        - `rtol-delta` (double): relative tolerance used to perturb properties
          when calculating numerical derivatives. The default value is 1e-8.
        """
        settings: CxxAnyMap
        self.kinetics.getDerivativeSettings(settings)
        return anymap_to_py(settings)

    @derivative_settings.setter
    def derivative_settings(self, settings: _DerivativeSettings) -> None:
        self.kinetics.setDerivativeSettings(py_to_anymap(settings))

    @property
    def forward_rate_constants_ddT(self) -> _Array:
        """
        Calculate derivatives for forward rate constants with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getFwdRateConstants_ddT)

    @property
    def forward_rate_constants_ddP(self) -> _Array:
        """
        Calculate derivatives for forward rate constants with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getFwdRateConstants_ddP)

    @property
    def forward_rate_constants_ddC(self) -> _Array:
        """
        Calculate derivatives for forward rate constants with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_reaction_array(self, kin_getFwdRateConstants_ddC)

    @property
    def forward_rates_of_progress_ddT(self) -> _Array:
        """
        Calculate derivatives for forward rates-of-progress with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getFwdRatesOfProgress_ddT)

    @property
    def forward_rates_of_progress_ddP(self) -> _Array:
        """
        Calculate derivatives for forward rates-of-progress with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getFwdRatesOfProgress_ddP)

    @property
    def forward_rates_of_progress_ddC(self) -> _Array:
        """
        Calculate derivatives for forward rates-of-progress with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_reaction_array(self, kin_getFwdRatesOfProgress_ddC)

    @property
    def forward_rates_of_progress_ddX(self) -> _Array:
        """
        Calculate derivatives for forward rates-of-progress with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_from_sparse(self.kinetics.fwdRatesOfProgress_ddX(),
                               self.n_reactions, self.n_total_species)

    @property
    def forward_rates_of_progress_ddCi(self) -> _Array:
        """
        Calculate derivatives for forward rates-of-progress with respect to species
        concentrations at constant temperature, pressure and remaining species
        concentrations.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`c_i`, all other :math:`c_j`
        are held constant.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.

        .. versionadded:: 3.0
        """
        return get_from_sparse(self.kinetics.fwdRatesOfProgress_ddCi(),
                               self.n_reactions, self.n_total_species)

    @property
    def reverse_rates_of_progress_ddT(self) -> _Array:
        """
        Calculate derivatives for reverse rates-of-progress with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getRevRatesOfProgress_ddT)

    @property
    def reverse_rates_of_progress_ddP(self) -> _Array:
        """
        Calculate derivatives for reverse rates-of-progress with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getRevRatesOfProgress_ddP)

    @property
    def reverse_rates_of_progress_ddC(self) -> _Array:
        """
        Calculate derivatives for reverse rates-of-progress with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_reaction_array(self, kin_getRevRatesOfProgress_ddC)

    @property
    def reverse_rates_of_progress_ddX(self) -> _Array:
        """
        Calculate derivatives for reverse rates-of-progress with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_from_sparse(self.kinetics.revRatesOfProgress_ddX(),
                               self.n_reactions, self.n_total_species)

    @property
    def reverse_rates_of_progress_ddCi(self) -> _Array:
        """
        Calculate derivatives for reverse rates-of-progress with respect to species
        concentrations at constant temperature, pressure and remaining species
        concentrations.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`c_i`, all other :math:`c_j`
        are held constant.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.

        .. versionadded:: 3.0
        """
        return get_from_sparse(self.kinetics.revRatesOfProgress_ddCi(),
                               self.n_reactions, self.n_total_species)

    @property
    def net_rates_of_progress_ddT(self) -> _Array:
        """
        Calculate derivatives for net rates-of-progress with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getNetRatesOfProgress_ddT)

    @property
    def net_rates_of_progress_ddP(self) -> _Array:
        """
        Calculate derivatives for net rates-of-progress with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        return get_reaction_array(self, kin_getNetRatesOfProgress_ddP)

    @property
    def net_rates_of_progress_ddC(self) -> _Array:
        """
        Calculate derivatives for net rates-of-progress with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_reaction_array(self, kin_getNetRatesOfProgress_ddC)

    @property
    def net_rates_of_progress_ddX(self) -> _Array:
        """
        Calculate derivatives for net rates-of-progress with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_from_sparse(self.kinetics.netRatesOfProgress_ddX(),
                               self.n_reactions, self.n_total_species)

    @property
    def net_rates_of_progress_ddCi(self) -> _Array:
        """
        Calculate derivatives for net rates-of-progress with respect to species
        concentrations at constant temperature, pressure and remaining species
        concentrations. For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`c_i`, all other :math:`c_j`
        are held constant.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.

        .. versionadded:: 3.0
        """
        return get_from_sparse(self.kinetics.netRatesOfProgress_ddCi(),
                               self.n_reactions, self.n_total_species)

    @property
    def creation_rates_ddT(self) -> _Array:
        """
        Calculate derivatives of species creation rates with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        return get_species_array(self, kin_getCreationRates_ddT)

    @property
    def creation_rates_ddP(self) -> _Array:
        """
        Calculate derivatives of species creation rates with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        return get_species_array(self, kin_getCreationRates_ddP)

    @property
    def creation_rates_ddC(self) -> _Array:
        """
        Calculate derivatives of species creation rates with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_species_array(self, kin_getCreationRates_ddC)

    @property
    def creation_rates_ddX(self) -> _Array:
        """
        Calculate derivatives for species creation rates with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_from_sparse(self.kinetics.creationRates_ddX(),
                               self.n_total_species, self.n_total_species)

    @property
    def creation_rates_ddCi(self) -> _Array:
        """
        Calculate derivatives for species creation rates with respect to species
        concentration at constant temperature, pressure, and concentration of all other
        species. For sparse output, set ``ct.use_sparse(True)``.

        The method returns a matrix with `n_total_species` rows and `n_total_species`
        columns.

        For a derivative with respect to :math:`c_i`, all other :math:`c_i` are
        held constant.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.

        .. versionadded:: 3.0
        """
        return get_from_sparse(self.kinetics.creationRates_ddCi(),
                               self.n_total_species, self.n_total_species)

    @property
    def destruction_rates_ddT(self) -> _Array:
        """
        Calculate derivatives of species destruction rates with respect to temperature
        at constant pressure, molar concentration and mole fractions.
        """
        return get_species_array(self, kin_getDestructionRates_ddT)

    @property
    def destruction_rates_ddP(self) -> _Array:
        """
        Calculate derivatives of species destruction rates with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        return get_species_array(self, kin_getDestructionRates_ddP)

    @property
    def destruction_rates_ddC(self) -> _Array:
        """
        Calculate derivatives of species destruction rates with respect to molar
        concentration at constant temperature, pressure and mole fractions.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_species_array(self, kin_getDestructionRates_ddC)

    @property
    def destruction_rates_ddX(self) -> _Array:
        """
        Calculate derivatives for species destruction rates with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_from_sparse(self.kinetics.destructionRates_ddX(),
                               self.n_total_species, self.n_total_species)

    @property
    def destruction_rates_ddCi(self) -> _Array:
        """
        Calculate derivatives for species destruction rates with respect to species
        concentration at constant temperature, pressure, and concentration of all other
        species. For sparse output, set ``ct.use_sparse(True)``.

        The method returns a matrix with `n_total_species` rows and `n_total_species` columns.
        For a derivative with respect to :math:`c_i`, all other :math:`c_i` are
        held constant.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.

        .. versionadded:: 3.0
        """
        return get_from_sparse(self.kinetics.destructionRates_ddCi(),
                               self.n_total_species, self.n_total_species)

    @property
    def net_production_rates_ddT(self) -> _Array:
        """
        Calculate derivatives of species net production rates with respect to
        temperature at constant pressure, molar concentration and mole fractions.
        """
        return get_species_array(self, kin_getNetProductionRates_ddT)

    @property
    def net_production_rates_ddP(self) -> _Array:
        """
        Calculate derivatives of species net production rates with respect to pressure
        at constant temperature, molar concentration and mole fractions.
        """
        return get_species_array(self, kin_getNetProductionRates_ddP)

    @property
    def net_production_rates_ddC(self) -> _Array:
        """
        Calculate derivatives of species net production rates with respect to molar
        density at constant temperature, pressure and mole fractions.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_species_array(self, kin_getNetProductionRates_ddC)

    @property
    def net_production_rates_ddX(self) -> _Array:
        """
        Calculate derivatives for species net production rates with respect to species
        concentrations at constant temperature, pressure and molar concentration.
        For sparse output, set ``ct.use_sparse(True)``.

        Note that for derivatives with respect to :math:`X_i`, all other :math:`X_j`
        are held constant, rather than enforcing :math:`\sum X_j = 1`.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_from_sparse(self.kinetics.netProductionRates_ddX(),
                               self.n_total_species, self.n_total_species)

    @property
    def net_production_rates_ddCi(self) -> _Array:
        """
        Calculate derivatives for species net production rates with respect to species
        concentration at constant temperature, pressure, and concentration of all other
        species. For sparse output, set ``ct.use_sparse(True)``.

        The method returns a matrix with `n_total_species` rows and `n_total_species` columns.
        For a derivative with respect to :math:`c_i`, all other :math:`c_i` are
        held constant.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.

        .. versionadded:: 3.0
        """
        return get_from_sparse(self.kinetics.netProductionRates_ddCi(),
                               self.n_total_species, self.n_total_species)

    @property
    def delta_enthalpy(self) -> _Array:
        """Change in enthalpy for each reaction [J/kmol]."""
        return get_reaction_array(self, kin_getDeltaEnthalpy)

    @property
    def delta_gibbs(self) -> _Array:
        """Change in Gibbs free energy for each reaction [J/kmol]."""
        return get_reaction_array(self, kin_getDeltaGibbs)

    @property
    def delta_entropy(self) -> _Array:
        """Change in entropy for each reaction [J/kmol/K]."""
        return get_reaction_array(self, kin_getDeltaEntropy)

    @property
    def delta_standard_enthalpy(self) -> _Array:
        """
        Change in standard-state enthalpy (independent of composition) for
        each reaction [J/kmol].
        """
        return get_reaction_array(self, kin_getDeltaSSEnthalpy)

    @property
    def third_body_concentrations(self) -> _Array:
        """
        Effective third-body concentrations used by individual reactions; values
        are only defined for reactions involving third-bodies and are set to
        not-a-number otherwise.
        """
        return get_reaction_array(self, kin_getThirdBodyConcentrations)

    @property
    def delta_standard_gibbs(self) -> _Array:
        """
        Change in standard-state Gibbs free energy (independent of composition)
        for each reaction [J/kmol].
        """
        return get_reaction_array(self, kin_getDeltaSSGibbs)

    @property
    def delta_standard_entropy(self) -> _Array:
        """
        Change in standard-state entropy (independent of composition) for
        each reaction [J/kmol/K].
        """
        return get_reaction_array(self, kin_getDeltaSSEntropy)

    @property
    def heat_release_rate(self) -> float:
        """
        Get the total volumetric heat release rate [W/m³].
        """
        return - np.sum(self.partial_molar_enthalpies *
                        self.net_production_rates, 0)

    @property
    def heat_production_rates(self) -> _Array:
        """
        Get the volumetric heat production rates [W/m³] on a per-reaction
        basis. The sum over all reactions results in the total volumetric heat
        release rate.
        Example: C. K. Law: Combustion Physics (2006), Fig. 7.8.6

        >>> gas.heat_production_rates[1]  # heat production rate of the 2nd reaction
        """
        return - self.net_rates_of_progress * self.delta_enthalpy


@cython.cclass
class InterfaceKinetics(Kinetics):
    """
    A kinetics manager for heterogeneous reaction mechanisms. The
    reactions are assumed to occur at an interface between bulk phases.
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
        if not init:
            return
        if pystr(self.kinetics.kineticsType()) not in ("surface", "edge"):
            raise TypeError("Underlying Kinetics class is not of the correct type.")
        self._setup_phase_indices()

    def _setup_phase_indices(self):
        self._phase_indices = {}
        for name, phase in list(self.adjacent.items()) + [(self.name, self)]:
            i = self.kinetics.phaseIndex(stringify(name), True)
            self._phase_indices[phase] = i
            self._phase_indices[name] = i
            self._phase_indices[i] = i

    def advance_coverages(self, dt: float, rtol: float = 1e-7, atol: float = 1e-14,
                          max_step_size: float = 0.0, max_steps: int = 20000,
                          max_error_test_failures: int = 7) -> None:
        """
        This method carries out a time-accurate advancement of the surface
        coverages for a specified amount of time.
        """
        cython.cast(cython.pointer(CxxInterfaceKinetics), self.kinetics).advanceCoverages(
            dt, rtol, atol, max_step_size, max_steps, max_error_test_failures)

    def advance_coverages_to_steady_state(self, loglevel: int = 0) -> None:
        """
        This method advances the surface coverages to steady state.
        """
        cython.cast(cython.pointer(CxxInterfaceKinetics),
                    self.kinetics).solvePseudoSteadyStateProblem(loglevel)

    def phase_index(self, phase: _ThermoPhase | str | int) -> int:
        """
        Get the index of the phase ``phase``, where ``phase`` may specified using
        the phase object, the name, or the index itself.
        """
        return self._phase_indices[phase]

    def _phase_slice(self, phase: _ThermoPhase | str | int) -> slice:
        p = self.phase_index(phase)
        k1 = self.kinetics_species_index(0, p)

        if p == self.n_phases - 1:
            k2 = self.n_total_species
        else:
            k2 = self.kinetics_species_index(0, p+1)

        return slice(k1, k2)

    def get_creation_rates(self, phase: _ThermoPhase | str | int) -> _Array:
        """
        Creation rates for each species in phase ``phase``. Use the
        `creation_rates` property to get the creation rates for species in all
        phases.
        """
        return self.creation_rates[self._phase_slice(phase)]

    def get_destruction_rates(self, phase: _ThermoPhase | str | int) -> _Array:
        """
        Destruction rates for each species in phase ``phase``. Use the
        `destruction_rates` property to get the destruction rates for species
        in all phases.
        """
        return self.destruction_rates[self._phase_slice(phase)]

    def get_net_production_rates(self, phase: _ThermoPhase | str | int) -> _Array:
        """
        Net production rates for each species in phase ``phase``. Use the
        `net_production_rates` property to get the net production rates for
        species in all phases.
        """
        return self.net_production_rates[self._phase_slice(phase)]

    def interface_current(self, phase: _ThermoPhase | str | int) -> float:
        """
        The interface current is useful when charge transfer reactions occur at
        an interface. It is defined here as the net positive charge entering the
        phase ``phase`` (Units: A/m² for a surface, A/m for an edge reaction).
        """
        i_phase = self.phase_index(phase)
        return cython.cast(cython.pointer(CxxInterfaceKinetics),
                           self.kinetics).interfaceCurrent(i_phase)

    def write_yaml(  # type: ignore[override]
        self,
        filename: _Path | str,
        phases: _Sequence[_ThermoPhase] | _ThermoPhase | None = None,
        units: _UnitSystem | _UnitDict | None = None,
        precision: int | None = None,
        skip_user_defined: bool | None = None,
    ) -> None:
        """
        See `_SolutionBase.write_yaml <cantera._cantera._SolutionBase.write_yaml>`.
        """
        if phases is not None:
            phases = list(phases)
        else:
            phases = []

        for phase in self._phase_indices:
            if isinstance(phase, _SolutionBase) and phase is not self:
                phases.append(phase)

        super().write_yaml(filename, phases, units, precision, skip_user_defined)
