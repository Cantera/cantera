# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3
# pyright: reportMissingImports=false, reportAttributeAccessIssue=false
# pyright: reportUndefinedVariable=false, reportUnboundVariable=false
# pyright: reportInvalidTypeArguments=false, reportAssignmentType=false
# pyright: reportIndexIssue=false, reportInvalidTypeForm=false

from collections.abc import Sequence as _Sequence
from typing import Literal as _Literal, overload as _overload

import numpy as np

import cython
from cython.cimports.cantera._utils import stringify, pystr
from cython.cimports.cantera.solutionbase import _SolutionBase

from .thermo import ThermoPhase as _ThermoPhase
from ._types import (
    Array as _Array,
    ArrayLike as _ArrayLike,
    EquilibriumSolver as _EquilibriumSolver,
    LogLevel as _LogLevel,
    PropertyPair as _PropertyPair,
)


@cython.cclass
class Mixture:
    """
    Class Mixture represents mixtures of one or more phases of matter.  To
    construct a mixture, supply a list of phases to the constructor, each
    paired with the number of moles for that phase::

        >>> gas = cantera.Solution("gas.yaml")
        >>> gas.species_names
        ['H2', 'H', 'O2', 'O', 'OH']
        >>> graphite = cantera.Solution("graphite.yaml")
        >>> graphite.species_names
        ['C(g)']
        >>> mix = cantera.Mixture([(gas, 1.0), (graphite, 0.1)])
        >>> mix.species_names
        ['H2', 'H', 'O2', 'O', 'OH', 'C(g)']

    Note that the objects representing each phase compute only the intensive
    state of the phase -- they do not store any information on the amount of
    this phase. Mixture objects, on the other hand, represent the full
    extensive state.

    .. caution::
       The Mixture class exists mainly for the purpose of providing input for multiphase
       equilibrium calculations. Its functionality for modifying the state of the
       mixture and computing properties is quite limited. Mixture objects cannot be used
       in conjunction with reactor networks.

       Furthermore, the multiphase equilibrium solvers currently have a number of
       problems that lead to solver failures or incorrect results for some inputs. See
       the `list of issues on GitHub <https://github.com/Cantera/cantera/issues?q=is%3Aopen+is%3Aissue+label%3AEquilibrium>`_
       for more information.

    Mixture objects are 'lightweight' in the sense that they do not store
    parameters needed to compute thermodynamic or kinetic properties of the
    phases. These are contained in the ('heavyweight') phase objects. Multiple
    mixture objects may be constructed using the same set of phase objects.
    Each one stores its own state information locally, and synchronizes the
    phases objects whenever it requires phase properties.
    """
    def __cinit__(self, phases):
        self._mix = make_shared[CxxMultiPhase]()
        self.mix = self._mix.get()
        self._phases = []

        phase: _SolutionBase
        if isinstance(phases[0], _SolutionBase):
            # Assign default composition to the list of phases
            phases = [(p, 1 if i == 0 else 0) for i, p in enumerate(phases)]

        for phase, moles in phases:
            self.mix.addPhase(phase.base.thermo(), moles)
            self._phases.append(phase)

        self.mix.init()
        if self._phases:
            self.P = self._phases[0].P
            self.T = self._phases[0].T

    def __init__(
        self,
        phases: _Sequence[tuple[_ThermoPhase, float]] | _Sequence[_ThermoPhase],
    ) -> None:
        # The C++ object is constructed in __cinit__; this typed __init__ exists so that
        # mypy/pyright (which do not recognize Cython's __cinit__) publish the constructor
        # signature.
        pass

    def report(self, threshold: float = 1e-14) -> str:
        """
        Generate a report describing the thermodynamic state of this mixture. To
        print the report to the screen, simply call the mixture object. The
        following two statements are equivalent::

        >>> mix()
        >>> print(mix.report())
        """
        self.mix.updatePhases()
        s = []
        for i, phase in enumerate(self._phases):
            s.append('************ Phase {0} ************'.format(phase.name))
            s.append('Moles:  {0}'.format(self.phase_moles(i)))
            s.append(phase.report(threshold=threshold))

        return '\n'.join(s)

    def __call__(self) -> None:
        print(self.report())

    @property
    def n_elements(self) -> int:
        """Total number of elements present in the mixture."""
        return self.mix.nElements()

    @cython.ccall
    def element_index(self, element: str | int) -> cython.int:
        """Index of element with name 'element'.

            >>> mix.element_index('H')
            2
        """
        if isinstance(element, (str, bytes)):
            return self.mix.elementIndex(stringify(element), True)

        if isinstance(element, (int, float)):
            return self.mix.checkElementIndex(cython.cast(cython.int, element))

        raise TypeError("'element' must be a string or a number. "
                        f"Got {element!r}.")

    @property
    def n_species(self) -> int:
        """Number of species."""
        return self.mix.nSpecies()

    def species_name(self, k: cython.int) -> str:
        """Name of the species with index ``k``. Note that index numbers
        are assigned in order as phases are added."""
        return pystr(self.mix.speciesName(k))

    @property
    def species_names(self) -> list[str]:
        """ Get the names of the species from all phases in the mixture """
        return [self.species_name(k) for k in range(self.n_species)]

    def species_index(self, phase: _ThermoPhase | str | int, species: str | int) -> int:
        """
        :param phase:
            Phase object, index or name
        :param species:
            Species name or index

        Returns the global index of species ``species`` in phase ``phase``.
        """
        p = self.phase_index(phase)

        if isinstance(species, (str, bytes, int, float)):
            k = self.phase(p).species_index(species)
        else:
            raise TypeError("'species' must be a string or number")

        return self.mix.speciesIndex(k, p)

    def n_atoms(self, k: int, m: str | int) -> float:
        """
        Number of atoms of element ``m`` in the species with global index ``k``.
        The element may be referenced either by name or by index.

        >>> n = mix.n_atoms(3, 'H')
        4.0
        """
        if not 0 <= k < self.n_species:
            raise IndexError('Species index ({0}) out of range (0 < {1})'.format(k, self.n_species))
        return self.mix.nAtoms(k, self.element_index(m))

    @property
    def n_phases(self) -> int:
        """Number of phases"""
        return len(self._phases)

    def phase(self, n: int) -> _ThermoPhase:
        """ Return the ThermoPhase object for phase number ``n`` in the mixture. """
        return self._phases[n]

    def phase_index(self, p: _ThermoPhase | str | int) -> int:
        """Index of the phase named ``p``."""
        if isinstance(p, _ThermoPhase):
            p = p.name

        if isinstance(p, (int, float)):
            if p == int(p) and 0 <= p < self.n_phases:
                return int(p)
            else:
                raise IndexError("Phase index '{0}' out of range.".format(p))
        elif isinstance(p, (str, bytes)):
            for i, phase in enumerate(self._phases):
                if phase.name == p:
                    return i
        raise KeyError("No such phase: '{0}'".format(p))

    @property
    def phase_names(self) -> list[str]:
        """Names of all phases in the order added."""
        return [phase.name for phase in self._phases]

    @property
    def T(self) -> float:
        """
        Get or set the Temperature [K] of all phases in the mixture. When set,
        the pressure of the mixture is held fixed.
        """
        return self.mix.temperature()

    @T.setter
    def T(self, T: float) -> None:
        self.mix.setTemperature(T)

    @property
    def min_temp(self) -> float:
        """
        The minimum temperature for which all species in multi-species
        solutions have valid thermo data. Stoichiometric phases are not
        considered in determining min_temp.
        """
        return self.mix.minTemp()

    @property
    def max_temp(self) -> float:
        """
        The maximum temperature for which all species in multi-species
        solutions have valid thermo data. Stoichiometric phases are not
        considered in determining max_temp.
        """
        return self.mix.maxTemp()

    @property
    def P(self) -> float:
        """Get or set the Pressure [Pa] of all phases in the mixture. When set,
         the temperature of the mixture is held fixed."""
        return self.mix.pressure()

    @P.setter
    def P(self, P: float) -> None:
        self.mix.setPressure(P)

    @property
    def charge(self) -> float:
        """The total charge in Coulombs, summed over all phases."""
        return self.mix.charge()

    def phase_charge(self, p: _ThermoPhase | str | int) -> float:
        """The charge of phase ``p`` in Coulombs."""
        return self.mix.phaseCharge(self.phase_index(p))

    @_overload
    def phase_moles(self, p: _ThermoPhase | str | int) -> float: ...
    @_overload
    def phase_moles(self, p: None = None) -> list[float]: ...
    def phase_moles(self, p=None):
        """
        Moles in phase ``p``, if ``p`` is specified, otherwise the number of
        moles in all phases.
        """
        if p is None:
            return [self.mix.phaseMoles(n) for n in range(self.n_phases)]
        else:
            return self.mix.phaseMoles(self.phase_index(p))

    def set_phase_moles(self, p: _ThermoPhase | str | int, moles: float) -> None:
        """
        Set the number of moles of phase ``p`` to ``moles``.
        """
        self.mix.setPhaseMoles(self.phase_index(p), moles)

    @property
    def species_moles(self) -> _Array:
        """
        Get or set the number of moles of each species. May be set either as a
        string or as an array. If an array is used, it must be dimensioned at
        least as large as the total number of species in the mixture. Note that
        the species may belong to any phase, and unspecified species are set to
        zero.

        >>> mix.species_moles = 'C(s):1.0, CH4:2.0, O2:0.2'
        """
        data = np.empty(self.n_species)
        for k in range(self.n_species):
            data[k] = self.mix.speciesMoles(k)
        return data

    @species_moles.setter
    def species_moles(self, moles: str | _ArrayLike) -> None:
        if isinstance(moles, (str, bytes)):
            self.mix.setMolesByName(stringify(moles))
            return

        data = np.ascontiguousarray(moles, dtype=np.double)
        cdata: cython.double[::1] = data
        self.mix.setMoles(span[cython.double](cython.address(cdata[0]),
                                              cython.cast(cython.size_t, data.size)))

    def element_moles(self, e: str | int) -> float:
        """
        Total number of moles of element ``e``, summed over all species.
        The element may be referenced either by index number or by name.
        """
        return self.mix.elementMoles(self.element_index(e))

    @property
    def chemical_potentials(self) -> _Array:
        """The chemical potentials of all species [J/kmol]."""
        data = np.empty(self.n_species)
        cdata: cython.double[::1] = data
        view: span[cython.double] = span[cython.double](
            cython.address(cdata[0]), cython.cast(cython.size_t, self.n_species))
        self.mix.getChemPotentials(view)
        return data

    def equilibrate(self, XY: _PropertyPair, solver: _EquilibriumSolver = "auto",
                    rtol: float = 1e-9, max_steps: int = 50000, max_iter: int = 100,
                    estimate_equil: _Literal[-1, 0, 1] = 0,
                    log_level: _LogLevel = 0) -> None:
        """
        Set to a state of chemical equilibrium holding property pair ``XY``
        constant. This method uses a version of the VCS algorithm to find the
        composition that minimizes the total Gibbs free energy of the mixture,
        subject to element conservation constraints. For a description of the
        theory, see Smith and Missen, "Chemical Reaction Equilibrium."

        :param XY:
            A two-letter string, which must be one of the set::

                ['TP', 'HP', 'SP']
        :param solver: Set to either 'auto', 'vcs', or 'gibbs' to choose
            implementation of the solver to use. 'vcs' uses the solver
            implemented in the C++ class 'VCSnonideal', 'gibbs' uses the one
            implemented in class 'MultiPhaseEquil'. 'auto' will try the 'vcs'
            solver first and then the 'gibbs' solver if that fails.
        :param rtol:
            Error tolerance. Iteration will continue until (Delta mu)/RT is
            less than this value for each reaction. Note that this default is
            very conservative, and good equilibrium solutions may be obtained
            with larger error tolerances.
        :param max_steps:
            Maximum number of steps to take while solving the equilibrium
            problem for specified *T* and *P*.
        :param max_iter:
            Maximum number of temperature and/or pressure iterations.
            This is only relevant if a property pair other than (T,P) is
            specified.
        :param estimate_equil:
            Flag indicating whether the solver should estimate its own initial
            condition. If 0, the initial mole fraction vector in the phase
            objects are used as the initial condition. If 1, the initial mole
            fraction vector is used if the element abundances are satisfied.
            if -1, the initial mole fraction vector is thrown out, and an
            estimate is formulated.
        :param log_level:
            Determines the amount of output displayed during the solution
            process. 0 indicates no output, while larger numbers produce
            successively more verbose information.
        """
        self.mix.equilibrate(stringify(XY.upper()), stringify(solver), rtol,
                             max_steps, max_iter, estimate_equil, log_level)
