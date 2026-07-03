# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3
# pyright: reportMissingImports=false, reportAttributeAccessIssue=false
# pyright: reportUndefinedVariable=false, reportUnboundVariable=false
# pyright: reportInvalidTypeArguments=false, reportAssignmentType=false
# pyright: reportIndexIssue=false, reportInvalidTypeForm=false

import warnings
from shutil import get_terminal_size as _get_terminal_size
from pathlib import Path as _Path
from collections.abc import Callable as _Callable, Sequence as _Sequence
from typing import (
    Any as _Any,
    ClassVar as _ClassVar,
    Literal as _Literal,
    TypeAlias as _TypeAlias,
    TypedDict as _TypedDict,
    overload as _overload,
    TYPE_CHECKING,
)
from typing_extensions import Never as _Never

import numpy as np

import cython
from cython.cimports.cython import view  # for view.array (sized-pointer→memoryview idiom)
import cython.cimports.numpy as cnp  # Required: triggers import_array() for numpy C-API
from cython.cimports.cantera._utils import stringify, pystr, anymap_to_py

from ._utils import CanteraError
from .interrupts import no_op
from .solutionbase import _SolutionBase
from .transport import _TransportModel

from ._types import (
    Array as _Array, ArrayLike as _ArrayLike, Basis as _Basis,
    CompositionLike as _CompositionLike, CompressionLevel as _CompressionLevel,
    LogLevel as _LogLevel, RefineCriteria as _RefineCriteria)

if TYPE_CHECKING:
    from .func1 import Func1 as _Func1

# Parametrized generics (tuple[...]) are coerced by Cython's annotation_typing just
# like bare builtins, rejecting a list where a tuple was published; route through a
# TypeAlias (not coerced) to keep the runtime accepting any sequence
_BoundsPair: _TypeAlias = tuple[float, float]

# `anymap_to_py` may return an `AnyMap` (a `dict` subclass) rather than a plain
# `dict`; an inline `dict[str, str]` return annotation is coerced by Cython 3 and
# rejects that subclass instance at runtime. Route through a TypeAlias (not coerced).
_RestoreMetadata: _TypeAlias = dict[str, str]

_ToleranceSettings = _TypedDict(
    "_ToleranceSettings",
    {
        "transient-abstol": float,
        "steady-abstol": float,
        "transient-reltol": float,
        "steady-reltol": float,
    },
)

class _PhaseSettings(_TypedDict):
    name: str
    source: str

class _FixedPointSettings(_TypedDict):
    location: float
    temperature: float

_Domain1DSettings = _TypedDict(
    "_Domain1DSettings",
    {
        "type": str,
        "points": int,
        "tolerances": _ToleranceSettings,
        "transport-model": _TransportModel,
        "phase": _PhaseSettings,
        "radiation-enabled": bool,
        "energy-enabled": bool,
        "Soret-enabled": bool,
        "flux-gradient-basis": _Literal[0, 1],
        "refine-criteria": _RefineCriteria,
        "fixed-point": _FixedPointSettings,
    },
)


@cython.cclass
class Domain1D:
    _domain_type: _ClassVar[str] = "none"
    have_user_tolerances: bool

    def __cinit__(self, phase: _SolutionBase, *args, **kwargs):
        self.domain = cython.NULL

    def __init__(self, phase: _SolutionBase, *args: _Any, **kwargs: _Any) -> None:
        if self.domain is cython.NULL:
            raise TypeError("Can't instantiate abstract class Domain1D.")

        self.gas = phase
        self.set_default_tolerances()

    @property
    def phase(self) -> _SolutionBase:
        """
        Phase describing the domain (that is, a gas phase or surface phase).
        """
        return self.gas

    @property
    def index(self) -> int:
        """
        Index of this domain in a stack. Returns -1 if this domain is not part
        of a stack.
        """
        return self.domain.domainIndex()

    @property
    def domain_type(self) -> str:
        """
        String indicating the domain implemented.
        """
        return pystr(self.domain.domainType())

    @property
    def n_components(self) -> int:
        """Number of solution components at each grid point."""
        return self.domain.nComponents()

    @property
    def n_points(self) -> int:
        """Number of grid points belonging to this domain."""
        return self.domain.nPoints()

    def _to_array(
        self, dest: SolutionArrayBase, normalize: cbool
    ) -> SolutionArrayBase:
        """
        Retrieve domain data as a `SolutionArray` object. Service method used by
        `FlameBase.to_array`, which adds information not available in Cython.

        .. versionadded:: 3.0
        """
        dest._base = self.domain.toArray(normalize)
        dest.base = dest._base.get()
        return dest

    def _from_array(self, arr: SolutionArrayBase) -> None:
        """
        Restore domain data from a `SolutionArray` object. Service method used by
        `FlameBase.from_array`.

        .. versionadded:: 3.0
        """
        self.domain.fromArray(arr._base)

    def component_name(self, n: cython.int) -> str:
        """Name of the nth component."""
        return pystr(self.domain.componentName(n))

    @property
    def component_names(self) -> list[str]:
        """List of the names of all components of this domain."""
        return [self.component_name(n) for n in range(self.n_components)]

    def component_index(self, name: str) -> int:
        """Index of the component with name 'name'"""
        return self.domain.componentIndex(stringify(name))

    def global_component_index(self, component: str, point: int) -> int:
        """
        The index of the component named ``component`` at grid point ``point``
        within the global solution vector of a containing `Sim1D`. Mirrors
        `ReactorNet.global_component_index`.

        .. versionadded:: 4.0
        """
        return self.domain.globalComponentIndex(stringify(component), point)

    def _has_component(self, name: str) -> bool:
        """Check whether `Domain1D` has component"""
        return self.domain.hasComponent(stringify(name))

    @_overload
    def info(
        self,
        keys: "_Sequence[str] | None" = None,
        rows: int = 10,
        width: int | None = None,
        display: _Literal[False] = False,
    ) -> str: ...
    @_overload
    def info(
        self,
        keys: "_Sequence[str] | None" = None,
        rows: int = 10,
        width: int | None = None,
        display: bool = True,
    ) -> None: ...
    def info(self, keys=None, rows=10, width=None, display=True):
        """
        Display or return a concise summary of a `Domain1D`.

        :param keys: List of components to be displayed; if `None`, all components are
            considered.
        :param rows: Maximum number of rendered rows.
        :param width: Maximum width of rendered output.
        :param display: If `True`, display result (default), otherwise, return a string.

        .. versionadded:: 3.2

        .. todo::

            Consolidate with `Sim1D.show`
        """
        cxx_keys: vector[string]
        if keys is not None:
            for key in keys:
                cxx_keys.push_back(stringify(key))
        if width is None:
            try:
                width = _get_terminal_size().columns
            except OSError:
                width = 100

        ret = pystr(self.domain.info(cxx_keys, rows, width))
        if not display:
            return ret
        print(ret)

    def update_state(self, loc: cython.int) -> None:
        """
        Set the state of the `Solution` object used for calculations to the temperature
        and composition at the point with index ``point``.
        """
        self.domain.updateState(loc)

    @property
    def grid(self) -> _Array:
        """The grid for this domain."""
        grid_span: span[const_double] = self.domain.grid()
        garr: view.array = view.array(shape=(grid_span.size(),),
                                      itemsize=cython.sizeof(cython.double), format="d",
                                      allocate_buffer=False)
        garr.data = cython.cast(cython.p_char, grid_span.data())
        return np.array(garr, copy=True)

    @grid.setter
    def grid(self, grid: _ArrayLike) -> None:
        grid_vec: vector[cython.double]
        for g in grid:
            grid_vec.push_back(g)
        self.domain.setupGrid(span[const_double](grid_vec.data(), grid_vec.size()))

    def value(self, component: str) -> float:
        """
        Component value at a boundary.

        :param component:
            component name

        >>> t = b.value("T")

        .. versionadded:: 3.2
        """
        return self.domain.value(stringify(component))

    def set_value(self, component: str, value: float) -> None:
        """
        Set the value of one component at a boundary.

        :param component:
            component name
        :param value:
            numerical value

        >>> b.set("T", 500.)

        .. versionadded:: 3.2
        """
        return self.domain.setValue(stringify(component), value)

    def values(self, component: str) -> _Array:
        """
        Retrieve spatial profile of a component.

        :param component:
            component name

        >>> T = d.values("T")

        .. versionadded:: 3.2
        """
        values: vector[cython.double] = self.domain.values(stringify(component))
        return np.asarray(values)

    def set_values(self, component: str, values: _Array) -> None:
        """
        Specify spatial profile of a component.

        :param component:
            component name
        :param values:
            array containing values

        >>> d.set_values("T", T)

        .. versionadded:: 3.2
        """
        values_arr = np.ascontiguousarray(values, dtype=np.double)
        cvalues: cython.double[::1] = values_arr
        self.domain.setValues(stringify(component),
                              span[const_double](cython.address(cvalues[0]),
                                                cython.cast(cython.size_t, values_arr.size)))

    def residuals(self, component: str) -> _Array:
        """
        Retrieve internal work array value at one point. After calling `Sim1D.eval`,
        this array contains the values of the residual function.

        :param component:
            component name

        >>> T = d.residuals("T")

        .. versionadded:: 3.2
        """
        values: vector[cython.double] = self.domain.residuals(stringify(component))
        return np.asarray(values)

    def set_profile(
        self, component: str, positions: _Sequence[float], values: _Sequence[float]
    ) -> None:
        """
        Set an initial estimate for a profile of one component in one domain.

        :param component:
            component name
        :param positions:
            sequence of relative positions, from 0 on the left to 1 on the right
        :param values:
            sequence of values at the relative positions specified in ``positions``

        >>> d.set_profile('T', [0.0, 0.2, 1.0], [400.0, 800.0, 1500.0])

        .. versionadded:: 3.2
        """
        pos_vec: vector[cython.double]
        val_vec: vector[cython.double]
        for p in positions:
            pos_vec.push_back(p)
        for v in values:
            val_vec.push_back(v)

        pos_span: span[const_double] = span[const_double](pos_vec.data(), pos_vec.size())
        val_span: span[const_double] = span[const_double](val_vec.data(), val_vec.size())
        self.domain.setProfile(stringify(component), pos_span, val_span)

    def set_flat_profile(self, component: str, value: float) -> None:
        """
        Set a flat profile for a component.

        :param component:
            component name
        :param v:
            value

        >>> d.set_flat_profile('u', -3.0)

        .. versionadded:: 3.2
        """
        self.domain.setFlatProfile(stringify(component), value)

    def set_bounds(
        self,
        *,
        default: "_BoundsPair | None" = None,
        Y: "_BoundsPair | None" = None,
        **kwargs: "_BoundsPair",
    ) -> None:
        """
        Set the lower and upper bounds on the solution.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (lower bound, upper bound) tuples as
        the values.  The keyword ``default`` may be used to specify default
        bounds for all unspecified components. The keyword ``Y`` can be used to
        stand for all species mass fractions in flow domains.

        >>> d.set_bounds(default=(0, 1), Y=(-1.0e-5, 2.0))
        """
        if default is not None:
            for n in range(self.n_components):
                self.domain.setBounds(n, default[0], default[1])

        if Y is not None:
            k0 = self.component_index(self.gas.species_name(0))
            for n in range(k0, k0 + self.gas.n_species):
                self.domain.setBounds(n, Y[0], Y[1])

        for name,(lower,upper) in kwargs.items():
            self.domain.setBounds(self.component_index(name), lower, upper)

    def set_steady_tolerances(
        self,
        *,
        default: "_BoundsPair | None" = None,
        Y: "_BoundsPair | None" = None,
        abs: "_BoundsPair | None" = None,
        rel: "_BoundsPair | None" = None,
        **kwargs: "_BoundsPair",
    ) -> None:
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with component names as
        keywords and (relative tolerance, absolute tolerance) tuples as the values.
        The keyword ``default`` may be used to specify default bounds for all
        unspecified components. The keyword ``Y`` can be used to stand for all
        species mass fractions in flow domains. Alternatively, the keywords
        ``abs`` and ``rel`` can be used to specify arrays for the absolute and
        relative tolerances for each solution component.
        """
        self.have_user_tolerances = True
        if default is not None:
            self.domain.setSteadyTolerances(default[0], default[1])

        if abs is not None or rel is not None:
            if rel is None:
                rel = self.steady_reltol()
            if abs is None:
                abs = self.steady_abstol()
            assert len(abs) == len(rel) == self.n_components
            for n, (r, a) in enumerate(zip(rel, abs)):
                self.domain.setSteadyTolerances(r, a, n)

        if Y is not None:
            k0 = self.component_index(self.gas.species_name(0))
            for n in range(k0, k0 + self.gas.n_species):
                self.domain.setSteadyTolerances(Y[0], Y[1], n)

        for name, (rtol, atol) in kwargs.items():
            self.domain.setSteadyTolerances(rtol, atol,
                                            self.component_index(name))

    def set_transient_tolerances(
        self,
        *,
        default: "_BoundsPair | None" = None,
        Y: "_BoundsPair | None" = None,
        abs: "_BoundsPair | None" = None,
        rel: "_BoundsPair | None" = None,
        **kwargs: "_BoundsPair",
    ) -> None:
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with component names as
        keywords and (relative tolerance, absolute tolerance) tuples as the values.
        The keyword ``default`` may be used to specify default bounds for all
        unspecified components. The keyword ``Y`` can be used to stand for all
        species mass fractions in flow domains. Alternatively, the keywords
        ``abs`` and ``rel`` can be used to specify arrays for the absolute and
        relative tolerances for each solution component.
        """
        self.have_user_tolerances = True
        if default is not None:
            self.domain.setTransientTolerances(default[0], default[1])

        if abs is not None or rel is not None:
            if rel is None:
                rel = self.transient_reltol()
            if abs is None:
                abs = self.transient_abstol()
            assert len(abs) == len(rel) == self.n_components
            for n, (r, a) in enumerate(zip(rel, abs)):
                self.domain.setTransientTolerances(r, a, n)

        if Y is not None:
            k0 = self.component_index(self.gas.species_name(0))
            for n in range(k0, k0 + self.gas.n_species):
                self.domain.setTransientTolerances(Y[0], Y[1], n)

        for name, (rtol, atol) in kwargs.items():
            self.domain.setTransientTolerances(rtol, atol,
                                               self.component_index(name))

    def set_default_tolerances(self) -> None:
        """
        Set all tolerances to their default values
        """
        self.set_steady_tolerances(default=(1e-4, 1e-9))
        self.set_transient_tolerances(default=(1e-4, 1e-11))
        self.have_user_tolerances = False

    def bounds(self, component: str) -> tuple[float, float]:
        """
        Return the (lower, upper) bounds for a solution component.

        >>> d.bounds('T')
        (200.0, 5000.0)
        """
        n = self.component_index(component)
        return self.domain.lowerBound(n), self.domain.upperBound(n)

    def tolerances(self, component: str) -> tuple[float, float]:
        """
        Return the (relative, absolute) error tolerances for a solution
        component.

        >>> rtol, atol = d.tolerances('u')
        """
        k = self.component_index(component)
        return self.domain.rtol(k), self.domain.atol(k)

    @_overload
    def steady_reltol(self, component: str) -> float: ...
    @_overload
    def steady_reltol(self, component: None = None) -> _Array: ...
    def steady_reltol(self, component=None):
        """
        Return the relative error tolerance for the steady state problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.steady_rtol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.steady_rtol(self.component_index(component))

    @_overload
    def steady_abstol(self, component: str) -> float: ...
    @_overload
    def steady_abstol(self, component: None = None) -> _Array: ...
    def steady_abstol(self, component=None):
        """
        Return the absolute error tolerance for the steady state problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.steady_atol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.steady_atol(self.component_index(component))

    @_overload
    def transient_reltol(self, component: str) -> float: ...
    @_overload
    def transient_reltol(self, component: None = None) -> _Array: ...
    def transient_reltol(self, component=None):
        """
        Return the relative error tolerance for the transient problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.transient_rtol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.transient_rtol(self.component_index(component))

    @_overload
    def transient_abstol(self, component: str) -> float: ...
    @_overload
    def transient_abstol(self, component: None = None) -> _Array: ...
    def transient_abstol(self, component=None):
        """
        Return the absolute error tolerance for the transient problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.transient_atol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.transient_atol(self.component_index(component))

    @property
    def jacobian_mode(self) -> str:
        """
        Method used to evaluate this domain's Jacobian columns: ``'auto'`` (default),
        ``'analytic'``, or ``'finite-difference'``.

        In ``'auto'`` mode the domain computes the species Jacobian columns analytically
        where supported and silently falls back to finite differences otherwise.
        ``'analytic'`` behaves the same but raises an exception if analytic evaluation
        was requested yet cannot be used -- because the kinetics object lacks the
        required composition derivatives, or multicomponent transport is active.
        ``'finite-difference'`` always uses finite differences.

        .. versionadded:: 4.0
        """
        return pystr(self.domain.jacobianMode())

    @jacobian_mode.setter
    def jacobian_mode(self, mode: str) -> None:
        self.domain.setJacobianMode(stringify(mode))

    @property
    def name(self) -> str:
        """ The name / id of this domain """
        return pystr(self.domain.id())

    @name.setter
    def name(self, name: str) -> None:
        self.domain.setID(stringify(name))

    def __reduce__(self) -> _Never:
        raise NotImplementedError('Domain1D object is not picklable')

    def __copy__(self) -> _Never:
        raise NotImplementedError('Domain1D object is not copyable')

    @property
    def settings(self) -> _Domain1DSettings:
        """
        Return comprehensive dictionary describing type, name, and simulation settings
        that are specific to domain types.

        .. versionchanged:: 3.0

            Added missing domain-specific simulation settings and updated structure.
        """
        arr: shared_ptr[CxxSolutionArray]
        arr = self.domain.toArray(False)
        return anymap_to_py(arr.get().meta())


@cython.cclass
class Boundary1D(Domain1D):
    """
    Base class for boundary domains.

    :param phase:
        The (gas) phase corresponding to the adjacent flow domain
    """
    def __cinit__(self, phase: _SolutionBase, *args, name="", **kwargs):
        if self._domain_type in {"None"}:
            self.boundary = cython.NULL
        else:
            self._domain = CxxNewDomain1D(
                stringify(self._domain_type), phase._base, stringify(name))
            self.domain = self._domain.get()
            self.boundary = cython.cast(cython.pointer(CxxBoundary1D), self.domain)

    def __init__(self, phase: _SolutionBase, name: str | None = None) -> None:
        if self.boundary is cython.NULL:
            raise TypeError("Can't instantiate abstract class Boundary1D.")
        Domain1D.__init__(self, phase, name=name)

    @property
    def T(self) -> float:
        """ The temperature [K] at this boundary. """
        return self.boundary.temperature()

    @T.setter
    def T(self, T: float) -> None:
        self.boundary.setTemperature(T)

    @property
    def mdot(self) -> float:
        """The mass flow rate per unit area [kg/s/m²]"""
        return self.boundary.mdot()

    @mdot.setter
    def mdot(self, mdot: float) -> None:
        self.boundary.setMdot(mdot)

    @property
    def X(self) -> _Array:
        """
        Species mole fractions at this boundary. May be set as either a string
        or as an array.
        """
        self.gas.TPY = self.gas.T, self.gas.P, self.Y
        return self.gas.X

    @X.setter
    def X(self, X: _CompositionLike) -> None:
        self.gas.TPX = None, None, X
        data = np.ascontiguousarray(self.gas.X, dtype=np.double)
        cdata: cython.double[::1] = data
        self.boundary.setMoleFractions(span[const_double](cython.address(cdata[0]),
                                                          cython.cast(cython.size_t, data.size)))

    @property
    def Y(self) -> _Array:
        """
        Species mass fractions at this boundary. May be set as either a string
        or as an array.
        """
        nsp: cython.int = self.boundary.nSpecies()
        Y = np.empty(nsp)
        k: cython.int
        for k in range(nsp):
            Y[k] = self.boundary.massFraction(k)
        return Y

    @Y.setter
    def Y(self, Y: _CompositionLike) -> None:
        self.gas.TPY = self.gas.T, self.gas.P, Y
        self.X = self.gas.X

    @property
    def spread_rate(self) -> float:
        """
        Get/set the tangential velocity gradient [1/s] at this boundary.
        """
        return self.boundary.spreadRate()

    @spread_rate.setter
    def spread_rate(self, s: float) -> None:
        self.boundary.setSpreadRate(s)


@cython.cclass
class Inlet1D(Boundary1D):
    """
    A one-dimensional inlet. Note that an inlet can only be a terminal
    domain - it must be either the leftmost or rightmost domain in a stack.
    """
    _domain_type = "inlet"


@cython.cclass
class Outlet1D(Boundary1D):
    """
    A one-dimensional outlet. An outlet imposes a zero-gradient boundary
    condition on the flow.
    """
    _domain_type = "outlet"


@cython.cclass
class OutletReservoir1D(Boundary1D):
    """
    A one-dimensional outlet into a reservoir.
    """
    _domain_type = "outlet-reservoir"


@cython.cclass
class SymmetryPlane1D(Boundary1D):
    """A symmetry plane."""
    _domain_type = "symmetry-plane"


@cython.cclass
class Surface1D(Boundary1D):
    """A solid surface."""
    _domain_type = "surface"


@cython.cclass
class ReactingSurface1D(Boundary1D):
    """A reacting solid surface.

    :param phase:
        The (surface) phase corresponding to the boundary

    .. versionchanged:: 3.0

        Starting in Cantera 3.0, parameter `phase` should reference surface instead of
        gas phase.
    """
    _domain_type = "reacting-surface"

    def __init__(self, phase: _SolutionBase, name: str | None = None) -> None:
        gas = None
        for val in phase._adjacent.values():
            if val.phase_of_matter == "gas":
                gas = val
                break
        if gas is None:
            raise CanteraError("ReactingSurface1D needs an adjacent gas phase")
        super().__init__(gas, name=name)

        self.surface = phase

    @property
    def phase(self) -> _SolutionBase:
        """
        Get the `Interface` object representing species and reactions on the surface
        """
        return self.surface

    @property
    def coverage_enabled(self) -> bool:
        """Controls whether or not to solve the surface coverage equations."""
        return cython.cast(cython.pointer(CxxReactingSurf1D), self.domain).coverageEnabled()

    @coverage_enabled.setter
    def coverage_enabled(self, value: bool) -> None:
        cython.cast(cython.pointer(CxxReactingSurf1D), self.domain).enableCoverageEquations(
            cython.cast(cbool, value))


@cython.cclass
class FlowBase(Domain1D):
    """ Base class for 1D flow domains """
    def __cinit__(self, phase: _SolutionBase, *args, name="", **kwargs):
        self._domain = CxxNewDomain1D(
            stringify(self._domain_type), phase._base, stringify(name))
        self.domain = self._domain.get()
        self.flow = cython.cast(cython.pointer(CxxFlow1D), self.domain)

    def __init__(self, *args: _Any, **kwargs: _Any) -> None:
        super().__init__(*args, **kwargs)
        self.P = self.gas.P
        self.flow.solveEnergyEqn()

    def __getattr__(self, name: str, /) -> _Array:
        # used to access fields by alias rather than canonical names;
        # also provides access to species that are valid Python identifiers
        # (replicates SolutionArray behavior)
        component_name = name
        if (not self._has_component(name) and
            self._has_component(name.replace("_", "-"))):
            component_name = name.replace("_", "-")

        if self._has_component(component_name):
            return self.values(component_name)

        raise AttributeError(
            f"{self.__class__.__name__!r} object has no attribute {name!r}")

    @property
    def P(self) -> float:
        """Pressure [Pa]"""
        return self.flow.pressure()

    @P.setter
    def P(self, P: float) -> None:
        self.flow.setPressure(P)

    @property
    def T(self) -> _Array:
        """
        Array containing the temperature [K] at each grid point.

        .. versionadded:: 3.2
        """
        return self.values("T")

    @property
    def velocity(self) -> _Array:
        """
        Array containing the velocity [m/s] normal to the flame at each point.

        .. versionadded:: 3.2
        """
        return self.values("velocity")

    @property
    def spread_rate(self) -> _Array:
        """
        Array containing the tangential velocity gradient [1/s] (that is, radial
        velocity divided by radius) at each point. Note: This value is named
        ``spreadRate`` in the C++ code and is only defined for axisymmetric flows.

        .. versionadded:: 3.2
        """
        return self.values("spread_rate")

    @property
    def radial_pressure_gradient(self) -> _Array:
        """
        Array containing the radial pressure gradient (1/r)(dP/dr) [N/m⁴] at
        each point. Note: This value is named ``Lambda`` in the C++ code and is only
        defined for axisymmetric flows.

        .. versionadded:: 3.2
        """
        return self.values("Lambda")

    @property
    def electric_field(self) -> _Array:
        """
        Array containing the electric field strength at each point.
        Note: This value is named ``eField`` in the C++ code and is only defined if
        the transport model is ``ionized-gas``.

        .. versionadded:: 3.2
        """
        return self.values("eField")

    @property
    def oxidizer_velocity(self) -> _Array:
        """
        Array containing the oxidizer velocity (right boundary velocity) [m/s] at
        each point.
        Note: This value is named ``Uo`` in the C++ code and is only defined when using
        two-point control.

        .. versionchanged:: 3.2
        """
        return self.values("Uo")

    @property
    def transport_model(self) -> str:
        """
        Get/set the transport model used for calculating transport properties.

        .. versionadded:: 3.0
        """
        return pystr(self.flow.transportModel())

    @transport_model.setter
    def transport_model(self, model: str) -> None:
        self.flow.setTransportModel(stringify(model))
        # ensure that transport remains accessible
        self.gas.transport = self.gas.base.transport().get()

    def set_default_tolerances(self) -> None:
        """
        Set all tolerances to their default values
        """
        super().set_default_tolerances()
        if self.transport_model != "ionized-gas":
            return

        chargetol = {}
        for S in self.gas.species():
            if S.composition == {'E': 1.0}:
                chargetol[S.name] = (1e-5, 1e-20)
            elif S.charge != 0:
                chargetol[S.name] = (1e-5, 1e-16)
        self.set_steady_tolerances(**chargetol)
        self.set_transient_tolerances(**chargetol)
        self.have_user_tolerances = False

    @property
    def soret_enabled(self) -> bool:
        """
        Determines whether or not to include diffusive mass fluxes due to the
        Soret effect. Enabling this option only works for multicomponent and
        mixture-averaged diffusion models.
        """
        return self.flow.withSoret()

    @soret_enabled.setter
    def soret_enabled(self, enable: bool) -> None:
        self.flow.enableSoret(cython.cast(cbool, enable))

    @property
    def flux_gradient_basis(self) -> _Basis:
        """
        Get/Set whether or not species diffusive fluxes are computed with
        respect to their mass fraction gradients ('mass')
        or mole fraction gradients ('molar', default) when
        using the mixture-averaged transport model.
        """
        if self.flow.fluxGradientBasis() == ThermoBasis.molar:
            return 'molar'
        else:
            return 'mass'

    @flux_gradient_basis.setter
    def flux_gradient_basis(self, basis: _Basis) -> None:
        if basis == 'molar':
            self.flow.setFluxGradientBasis(ThermoBasis.molar)
        elif basis == 'mass':
            self.flow.setFluxGradientBasis(ThermoBasis.mass)
        else:
            raise ValueError("Valid choices are 'mass' or 'molar'."
                             " Got {!r}.".format(basis))

    @property
    def energy_enabled(self) -> bool:
        """ Determines whether or not to solve the energy equation."""
        return self.flow.doEnergy(0)

    @energy_enabled.setter
    def energy_enabled(self, enable: bool) -> None:
        if enable:
            self.flow.solveEnergyEqn()
        else:
            self.flow.fixTemperature()

    def set_fixed_temp_profile(self, pos: _ArrayLike, T: _ArrayLike) -> None:
        """Set the fixed temperature profile. This profile is used
        whenever the energy equation is disabled.

        :param pos:
            array of relative positions from 0 to 1
        :param temp:
            array of temperature values

        >>> d.set_fixed_temp_profile(array([0.0, 0.5, 1.0]),
        ...                          array([500.0, 1500.0, 2000.0])
        """
        x: vector[cython.double]
        y: vector[cython.double]
        for p in pos:
            x.push_back(p)
        for t in T:
            y.push_back(t)
        xspan: span[const_double] = span[const_double](x.data(), x.size())
        yspan: span[const_double] = span[const_double](y.data(), y.size())
        self.flow.setFixedTempProfile(xspan, yspan)

    def get_settings3(self) -> _Domain1DSettings:
        """
        Temporary method returning new behavior of settings getter.

        .. versionadded:: 3.0
        """
        return self.settings

    @property
    def boundary_emissivities(self) -> tuple[float, float]:
        """ Set/get boundary emissivities. """
        return self.flow.leftEmissivity(), self.flow.rightEmissivity()

    @boundary_emissivities.setter
    def boundary_emissivities(self, epsilon: "_BoundsPair") -> None:
        if len(epsilon) != 2:
            raise ValueError('Setting the boundary emissivities requires a '
                             'tuple of length 2.')
        self.flow.setBoundaryEmissivities(epsilon[0], epsilon[1])

    @property
    def radiation_enabled(self) -> bool:
        """ Determines whether or not to include radiative heat transfer """
        return self.flow.radiationEnabled()

    @radiation_enabled.setter
    def radiation_enabled(self, do_radiation: bool) -> None:
        self.flow.enableRadiation(cython.cast(cbool, do_radiation))

    @property
    def radiative_heat_loss(self) -> _Array:
        """
        Return radiative heat loss (only non-zero if radiation is enabled).
        """
        j: cython.int
        data = np.empty(self.n_points)
        for j in range(self.n_points):
            data[j] = self.flow.radiativeHeatLoss(j)
        return data

    def set_free_flow(self) -> None:
        """
        Set flow configuration for freely-propagating flames, using an internal
        point with a fixed temperature as the condition to determine the inlet
        mass flux.
        """
        self.flow.setFreeFlow()

    def set_axisymmetric_flow(self) -> None:
        """
        Set flow configuration for axisymmetric counterflow or burner-stabilized
        flames, using specified inlet mass fluxes.
        """
        self.flow.setAxisymmetricFlow()

    @property
    def type(self) -> str:
        """
        Return the type of flow domain being represented.

        Examples:
        - ``free-flow``/``free-ion-flow``,
        - ``axisymmetric-flow``/``axisymmetric-ion-flow``,
        - ``unstrained-flow``/``unstrained-ion-flow``
        """
        return pystr(self.flow.domainType())

    @property
    def electric_field_enabled(self) -> bool:
        """
        Determines whether or not to solve the electric field equation (only relevant
        if transport model is ``ionized-gas``).
        """
        return self.flow.doElectricField()

    @electric_field_enabled.setter
    def electric_field_enabled(self, enable: bool) -> None:
        if enable:
            self.flow.solveElectricField()
        else:
            self.flow.fixElectricField()

    @property
    def left_control_point_temperature(self) -> float:
        """ Get/Set the left control point temperature [K] """
        return self.flow.leftControlPointTemperature()

    @left_control_point_temperature.setter
    def left_control_point_temperature(self, T: float) -> None:
        self.flow.setLeftControlPointTemperature(T)

    @property
    def right_control_point_temperature(self) -> float:
        """ Get/Set the right control point temperature [K] """
        return self.flow.rightControlPointTemperature()

    @right_control_point_temperature.setter
    def right_control_point_temperature(self, T: float) -> None:
        self.flow.setRightControlPointTemperature(T)

    @property
    def left_control_point_coordinate(self) -> float:
        """ Get the left control point coordinate [m] """
        return self.flow.leftControlPointCoordinate()

    @property
    def right_control_point_coordinate(self) -> float:
        """ Get the right control point coordinate [m] """
        return self.flow.rightControlPointCoordinate()

    @property
    def two_point_control_enabled(self) -> bool:
        """ Get/Set the state of the two-point flame control """
        return self.flow.twoPointControlEnabled()

    @two_point_control_enabled.setter
    def two_point_control_enabled(self, enable: bool) -> None:
        self.flow.enableTwoPointControl(cython.cast(cbool, enable))


@cython.cclass
class FreeFlow(FlowBase):
    r"""A free flow domain. The equations solved are standard equations for adiabatic
    one-dimensional flow. The solution variables are:

    *velocity*
        axial velocity
    *T*
        temperature
    *Y_k*
        species mass fractions
    """
    _domain_type = "free-flow"


@cython.cclass
class UnstrainedFlow(FlowBase):
    r"""An unstrained flow domain. The equations solved are standard equations for
    adiabatic one-dimensional flow. The solution variables are:

    *velocity*
        axial velocity
    *T*
        temperature
    *Y_k*
        species mass fractions
    """
    _domain_type = "unstrained-flow"


@cython.cclass
class AxisymmetricFlow(FlowBase):
    r"""
    An axisymmetric flow domain. The equations solved are the similarity equations for
    the flow in a finite-height gap of infinite radial extent. The solution variables
    are:

    *velocity*
        axial velocity
    *spread_rate*
        radial velocity divided by radius
    *T*
        temperature
    *lambda*
        :math:`(1/r)(dP/dr)`
    *Y_k*
        species mass fractions

    It may be shown that if the boundary conditions on these variables are independent
    of radius, then a similarity solution to the exact governing equations exists in
    which these variables are all independent of radius. This solution holds only in
    the low-Mach-number limit, in which case :math:`(dP/dz) = 0`, and :math:`\Lambda` is
    a constant. (Lambda is treated as a spatially-varying solution variable for
    numerical reasons, but in the final solution it is always independent of :math:`z`.)
    As implemented here, the governing equations assume an ideal gas mixture. Arbitrary
    chemistry is allowed, as well as arbitrary variation of the transport properties.
    """
    _domain_type = "axisymmetric-flow"


@cython.cclass
class Sim1D:
    """
    Class Sim1D is a container for one-dimensional domains. It also holds the
    multi-domain solution vector, and controls the process of finding the
    solution.

    Domains are ordered left-to-right, with domain number 0 at the left.
    """

    domains: "tuple[Domain1D, ...]"

    def __init__(
        self, domains: _Sequence[Domain1D], *args: _Any, **kwargs: _Any
    ) -> None:
        cxx_domains: vector[shared_ptr[CxxDomain1D]]
        d: Domain1D
        for d in domains:
            cxx_domains.push_back(d._domain)

        self._sim = CxxNewSim1D(cxx_domains)
        self.sim = self._sim.get()
        self.domains = tuple(domains)
        self.set_interrupt(no_op)
        self._initialized = False
        self._initial_guess_args = ()
        self._initial_guess_kwargs = {}

    def set_interrupt(self, f: "_Callable[[float], float] | None") -> None:
        """
        Set an interrupt function to be called each time that :ct:`OneDim::eval` is
        called. The signature of ``f`` is ``float f(float)``. The default
        interrupt function is used to trap `KeyboardInterrupt` exceptions so
        that ``ctrl-c`` can be used to break out of the C++ solver loop.
        """
        if f is None:
            self.sim.setInterrupt(cython.NULL)
            self._interrupt = None
            return

        if not isinstance(f, Func1):
            f = Func1(f)
        self._interrupt = f
        self.sim.setInterrupt(self._interrupt.func)

    def set_time_step_callback(self, f: "_Callable[[float], float] | None") -> None:
        """
        Set a callback function to be called after each successful timestep.
        The signature of ``f`` is ``float f(float)``. The argument passed to ``f`` is
        the size of the timestep. The output is ignored.
        """
        if f is None:
            self.sim.setTimeStepCallback(cython.NULL)
            self._time_step_callback = None
            return

        if not isinstance(f, Func1):
            f = Func1(f)
        self._time_step_callback = f
        self.sim.setTimeStepCallback(self._time_step_callback.func)

    def set_steady_callback(self, f: "_Callable[[float], float] | None") -> None:
        """
        Set a callback function to be called after each successful steady-state
        solve, before regridding. The signature of ``f`` is ``float f(float)``. The
        argument passed to ``f`` is 0.0 and the output is ignored.
        """
        if f is None:
            self.sim.setSteadyCallback(cython.NULL)
            self._steady_callback = None
            return

        if not isinstance(f, Func1):
            f = Func1(f)
        self._steady_callback = f
        self.sim.setSteadyCallback(self._steady_callback.func)

    def domain_index(self, dom: Domain1D | str | int) -> int:
        """
        Get the index of a domain, specified either by name or as a Domain1D object.
        """
        if isinstance(dom, Domain1D):
            return self.domains.index(dom)
        if isinstance(dom, (str, bytes)):
            return self.sim.domainIndex(stringify(dom))
        assert 0 <= dom < len(self.domains)
        return dom

    def _get_indices(self, dom, comp):
        idom = self.domain_index(dom)
        dom = self.domains[idom]
        if isinstance(comp, (str, bytes)):
            kcomp = dom.component_index(comp)
        else:
            kcomp = comp

        assert 0 <= kcomp < dom.n_components

        return idom, kcomp

    def eval(self, rdt: float = 0.0) -> None:
        """
        Evaluate the governing equations using the current solution estimate,
        storing the residual in the array which is accessible with the
        `Domain1D.residuals` function.

        :param rdt:
           Reciprocal of the time-step
        """
        self.sim.eval(rdt)

    def eval_jacobian(self) -> None:
        """
        Evaluate the steady-state Jacobian at the current solution estimate.
        The result can be examined through the `linear_solver` object, e.g. the
        ``jacobian`` property of `EigenSparseDirectJacobian`.

        .. versionadded:: 4.0
        """
        self.sim.evalSSJacobian()

    def show(self) -> None:
        """
        Print the current solution.

        .. todo::

            Consolidate with `Domain1D.info`
        """
        if not self._initialized:
            self.set_initial_guess()
        self.sim.show()

    def set_time_step(self, stepsize: float, n_steps: _Sequence[int]) -> None:
        """Set the sequence of time steps to try when Newton fails.

        :param stepsize:
            initial time step size [s]
        :param n_steps:
            sequence of integer step numbers

        >>> s.set_time_step(1.0e-5, [1, 2, 5, 10])
        """
        data: vector[cython.int]
        for n in n_steps:
            data.push_back(n)
        self.sim.setTimeStep(stepsize, span[int](data))

    @property
    def max_time_step_count(self) -> int:
        """
        Get/Set the maximum number of time steps allowed before reaching the
        steady-state solution
        """
        return self.sim.maxTimeStepCount()

    @max_time_step_count.setter
    def max_time_step_count(self, nmax: int) -> None:
        self.sim.setMaxTimeStepCount(nmax)

    def set_jacobian_perturbation(
        self, relative: float, absolute: float, threshold: float
    ) -> None:
        """
        Configure perturbations used to evaluate finite difference Jacobian

        :param relative:
            Relative perturbation (multiplied by the absolute value of each component).
            Default ``1.0e-5``.
        :param absolute:
            Absolute perturbation (independent of component value). Default ``1.0e-10``.
        :param threshold:
            Threshold below which to exclude elements from the Jacobian. Default ``0.0``.
        """
        self.sim.setJacobianPerturbation(relative, absolute, threshold)

    @property
    def linear_solver(self) -> SystemJacobian:
        """
        Get/Set the the linear solver used to hold the Jacobian matrix and solve linear
        systems as part of each Newton iteration. The default is a banded, direct
        solver. See :ref:`sec-python-jacobians` for available solvers.
        """
        return SystemJacobian.wrap(self.sim.linearSolver())

    @linear_solver.setter
    def linear_solver(self, precon: SystemJacobian) -> None:
        self.sim.setLinearSolver(precon._base)

    def set_initial_guess(self, *args: _Any, **kwargs: _Any) -> None:
        """
        Store arguments for initial guess and prepare storage for solution.
        """
        self._initial_guess_args = args
        self._initial_guess_kwargs = kwargs
        self._get_initial_solution()
        self._initialized = True

    def _get_initial_solution(self):
        """
        Load the initial solution from each domain into the global solution
        vector.
        """
        self.sim.resize()
        self.sim.getInitialSoln()

    def extinct(self) -> bool:
        """
        Method overloaded for some flame types to indicate if the flame has been
        extinguished. Base class method always returns 'False'
        """
        return False

    def solve(
        self, loglevel: _LogLevel = 1, refine_grid: bool = True, auto: bool = False
    ) -> None:
        """
        Solve the problem.

        :param loglevel:
            integer flag controlling the amount of diagnostic output. Zero
            suppresses all output, and 5 produces very verbose output.
        :param refine_grid:
            if True, enable grid refinement.
        :param auto: if True, sequentially execute the different solution stages
            and attempt to automatically recover from errors. Attempts to first
            solve on the initial grid with energy enabled. If that does not
            succeed, a fixed-temperature solution will be tried followed by
            enabling the energy equation, and then with grid refinement enabled.
            If non-default tolerances have been specified or multicomponent
            transport is enabled, an additional solution using these options
            will be calculated.
        """

        if not auto:
            if not self._initialized:
                self.set_initial_guess()
            self.sim.solve(loglevel, cython.cast(cbool, refine_grid))
            return

        def set_transport(multi):
            for dom in self.domains:
                if isinstance(dom, FlowBase):
                    dom.transport_model = multi

        # Do initial solution steps with default tolerances
        have_user_tolerances = any(dom.have_user_tolerances for dom in self.domains)
        if have_user_tolerances:
            # Save the user-specified tolerances
            atol_ss_final = [dom.steady_abstol() for dom in self.domains]
            rtol_ss_final = [dom.steady_reltol() for dom in self.domains]
            atol_ts_final = [dom.transient_abstol() for dom in self.domains]
            rtol_ts_final = [dom.transient_reltol() for dom in self.domains]

        def restore_tolerances():
            if have_user_tolerances:
                for i in range(len(self.domains)):
                    self.domains[i].set_steady_tolerances(abs=atol_ss_final[i],
                                                        rel=rtol_ss_final[i])
                    self.domains[i].set_transient_tolerances(abs=atol_ts_final[i],
                                                            rel=rtol_ts_final[i])

        for dom in self.domains:
            dom.set_default_tolerances()

        # Do initial steps without Soret diffusion
        soret_doms = [dom for dom in self.domains if getattr(dom, 'soret_enabled', False)]

        def set_soret(soret):
            for dom in soret_doms:
                dom.soret_enabled = soret

        set_soret(False)

        # Do initial solution steps without multicomponent transport
        transport = self.transport_model
        solve_multi = self.transport_model == 'multicomponent'
        if solve_multi:
            set_transport('mixture-averaged')

        def log(msg, *args):
            if loglevel:
                print('\n{:*^78s}'.format(' ' + msg.format(*args) + ' '))

        flow_domains = [D for D in self.domains if isinstance(D, FlowBase)]
        zmin = [D.grid[0] for D in flow_domains]
        zmax = [D.grid[-1] for D in flow_domains]

        # 'data' entry is used for restart
        data = self._initial_guess_kwargs.get('data')
        if data is not None:
           nPoints = [len(flow_domains[0].grid)]
        else:
           nPoints = [len(flow_domains[0].grid), 12, 24, 48]

        for N in nPoints:
            for i,D in enumerate(flow_domains):
                if N > self.get_max_grid_points(D):
                    raise CanteraError('Maximum number of grid points exceeded')

                if N != len(D.grid):
                    D.grid = np.linspace(zmin[i], zmax[i], N)

            if data is None:
                self.set_initial_guess(*self._initial_guess_args,
                                       **self._initial_guess_kwargs)

            # Try solving with energy enabled, which usually works
            log('Solving on {} point grid with energy equation enabled', N)
            self.energy_enabled = True
            try:
                self.sim.solve(loglevel, cython.cast(cbool, False))
                solved = True
            except CanteraError as e:
                log(str(e))
                solved = False
            except Exception as e:
                # restore settings before re-raising exception
                set_transport(transport)
                set_soret(True)
                restore_tolerances()
                raise e

            # If initial solve using energy equation fails, fall back on the
            # traditional fixed temperature solve followed by solving the energy
            # equation
            if not solved or self.extinct():
                if self.extinct():
                    self.set_initial_guess(*self._initial_guess_args,
                                           **self._initial_guess_kwargs)
                log('Initial solve failed; Retrying with energy equation disabled')
                self.energy_enabled = False
                try:
                    self.sim.solve(loglevel, cython.cast(cbool, False))
                    solved = True
                except CanteraError as e:
                    log(str(e))
                    solved = False
                except Exception as e:
                    # restore settings before re-raising exception
                    set_transport(transport)
                    set_soret(True)
                    restore_tolerances()
                    raise e

                if solved:
                    log('Solving on {} point grid with energy equation re-enabled', N)
                    self.energy_enabled = True
                    try:
                        self.sim.solve(loglevel, cython.cast(cbool, False))
                        solved = True
                    except CanteraError as e:
                        log(str(e))
                        solved = False
                    except Exception as e:
                        # restore settings before re-raising exception
                        set_transport(transport)
                        set_soret(True)
                        restore_tolerances()
                        raise e

            if solved and not self.extinct() and refine_grid:
                # Found a non-extinct solution on the fixed grid
                log('Solving with grid refinement enabled')
                try:
                    self.sim.solve(loglevel, cython.cast(cbool, True))
                    solved = True
                except CanteraError as e:
                    log(str(e))
                    solved = False
                except Exception as e:
                    # restore settings before re-raising exception
                    set_transport(transport)
                    set_soret(True)
                    restore_tolerances()
                    raise e

                if solved and not self.extinct():
                    # Found a non-extinct solution on the refined grid
                    break

            if self.extinct():
                log('Flame is extinct on {} point grid', N)

            if not refine_grid:
                break

        if not solved:
            raise CanteraError('Could not find a solution for the 1D problem')

        if solve_multi:
            log('Solving with multicomponent transport')
            set_transport('multicomponent')

        if soret_doms:
            log('Solving with Soret diffusion')
            set_soret(True)

        if have_user_tolerances:
            log('Solving with user-specified tolerances')
            restore_tolerances()

        # Final call with expensive options enabled
        if have_user_tolerances or solve_multi or soret_doms:
            self.sim.solve(loglevel, cython.cast(cbool, refine_grid))

    def refine(self, loglevel: _LogLevel = 1) -> None:
        """
        Refine the grid, adding points where solution is not adequately
        resolved.
        """
        self.sim.refine(loglevel)

    def set_refine_criteria(
        self,
        domain: Domain1D | str | int,
        ratio: float = 10.0,
        slope: float = 0.8,
        curve: float = 0.8,
        prune: float = 0.05,
    ) -> None:
        """
        Set the criteria used to refine one domain.

        :param domain:
            domain object, index, or name
        :param ratio:
            additional points will be added if the ratio of the spacing on
            either side of a grid point exceeds this value
        :param slope:
            maximum difference in value between two adjacent points, scaled by
            the maximum difference in the profile (0.0 < slope < 1.0). Adds
            points in regions of high slope.
        :param curve:
            maximum difference in slope between two adjacent intervals, scaled
            by the maximum difference in the profile (0.0 < curve < 1.0). Adds
            points in regions of high curvature.
        :param prune:
            if the slope or curve criteria are satisfied to the level of
            'prune', the grid point is assumed not to be needed and is removed.
            Set prune significantly smaller than 'slope' and 'curve'. Set to
            zero to disable pruning the grid.

        >>> s.set_refine_criteria(d, ratio=5.0, slope=0.2, curve=0.3, prune=0.03)
        """
        idom = self.domain_index(domain)
        self.sim.setRefineCriteria(idom, ratio, slope, curve, prune)

    def get_refine_criteria(self, domain: Domain1D | str | int) -> _RefineCriteria:
        """
        Get a dictionary of the criteria used to refine one domain. The items in
        the dictionary are the ``ratio``, ``slope``, ``curve``, and ``prune``,
        as defined in `~Sim1D.set_refine_criteria`.

        :param domain:
            domain object, index, or name

        >>> s.set_refine_criteria(d, ratio=5.0, slope=0.2, curve=0.3, prune=0.03)
        >>> s.get_refine_criteria(d)
        {'ratio': 5.0, 'slope': 0.2, 'curve': 0.3, 'prune': 0.03}
        """
        idom = self.domain_index(domain)
        c = self.sim.getRefineCriteria(idom)
        return {'ratio': c[0], 'slope': c[1], 'curve': c[2], 'prune': c[3]}

    def set_grid_min(
        self, dz: float, domain: Domain1D | str | int | None = None
    ) -> None:
        """
        Set the minimum grid spacing on ``domain``. If ``domain`` is `None`, then
        set the grid spacing for all domains.
        """
        if domain is None:
            idom = -1
        else:
            idom = self.domain_index(domain)
        self.sim.setGridMin(idom, dz)

    def set_max_jac_age(self, ss_age: int, ts_age: int) -> None:
        """
        Set the maximum number of times the Jacobian will be used before it
        must be re-evaluated.

        :param ss_age:
            age criterion during steady-state mode
        :param ts_age:
            age criterion during time-stepping mode
        """
        self.sim.setJacAge(ss_age, ts_age)

    def set_time_step_factor(self, tfactor: float) -> None:
        """
        Set the factor by which the time step will be decreased after an
        unsuccessful step.

        :param tfactor:
            Multiplicative reduction factor applied after failed steps.
        """
        self.sim.setTimeStepFactor(tfactor)

    @property
    def time_step_growth_factor(self) -> float:
        """
        Get/Set the factor by which the time step will be increased after a
        successful step when the Jacobian is reused.

        This value is used directly by the ``"fixed-growth"`` strategy, and as
        the accepted growth factor or cap value for the other growth
        strategies.

        :param tfactor:
            Finite growth factor >= 1.0. The default value is 1.5.
        """
        return self.sim.timeStepGrowthFactor()

    @time_step_growth_factor.setter
    def time_step_growth_factor(self, tfactor: float) -> None:
        self.sim.setTimeStepGrowthFactor(tfactor)

    @property
    def time_step_growth_strategy(self) -> str:
        """
        Get/Set the strategy used to grow the time step after a successful
        step that reuses the Jacobian.

        Available options are:

        ``"fixed-growth"``:
            Always apply ``time_step_growth_factor``.
        ``"steady-norm"``:
            Apply growth only when the steady-state residual norm decreases.
        ``"transient-residual"``:
            Apply growth only when the transient residual norm decreases.
        ``"residual-ratio"``:
            Scale the growth factor based on transient residual improvement,
            capped by ``time_step_growth_factor``.
        ``"newton-iterations"``:
            Apply growth only if the most recent Newton solve used at most
            three iterations.
        """
        return pystr(self.sim.timeStepGrowthStrategy())

    @time_step_growth_strategy.setter
    def time_step_growth_strategy(self, strategy: str) -> None:
        self.sim.setTimeStepGrowthStrategy(stringify(strategy))

    @property
    def time_step_regrid(self) -> int:
        """
        Get/Set the maximum number of regrid attempts after a time step
        failure.

        This fallback is used by :meth:`Sim1D.solve` when ``refine_grid=True``.
        Set to zero to disable regrid-on-failure retries. The default value is
        3.

        :param max_tries:
            Maximum retry attempts. Values less than zero are invalid.
        """
        return self.sim.timeStepRegridMax()

    @time_step_regrid.setter
    def time_step_regrid(self, max_tries: int) -> None:
        self.sim.setTimeStepRegridMax(max_tries)

    def set_min_time_step(self, tsmin: float) -> None:
        """ Set the minimum time step. """
        self.sim.setMinTimeStep(tsmin)

    def set_max_time_step(self, tsmax: float) -> None:
        """ Set the maximum time step. """
        self.sim.setMaxTimeStep(tsmax)

    @property
    def fixed_temperature(self) -> float:
        """
        Set the temperature used to fix the spatial location of a freely
        propagating flame.
        """
        return self.sim.fixedTemperature()

    @fixed_temperature.setter
    def fixed_temperature(self, T: float) -> None:
        self.sim.setFixedTemperature(T)

    @property
    def fixed_temperature_location(self) -> float:
        """
        Return the location of the point where temperature is fixed for a freely
        propagating flame.
        """
        return self.sim.fixedTemperatureLocation()

    def set_left_control_point(self, T: float) -> None:
        """
        Set the left control point using the specified temperature. This user-provided
        temperature will be used to locate the closest grid point to that temperature,
        which will serve to locate the left control point's coordinate. Starting from
        the left boundary, the first grid point that is equal to or exceeds the
        specified temperature will be used to locate the left control point's
        coordinate.
        """
        self.sim.setLeftControlPoint(T)

    def set_right_control_point(self, T: float) -> None:
        """
        Set the right control point using a specified temperature. This user-provided
        temperature will be used to locate the closest grid point to that temperature,
        which will serve to locate the right control point's coordinate.Starting from
        the right boundary, the first grid point that is equal to or exceeds the
        specified temperature will be used to locate the right control point's
        coordinate.
        """
        self.sim.setRightControlPoint(T)

    def save(
        self,
        filename: _Path | str = 'soln.yaml',
        name: str = 'solution',
        description: str | None = None,
        loglevel: _LogLevel | None = None,
        *,
        overwrite: bool = False,
        compression: _CompressionLevel = 0,
        basis: _Basis | _Literal["X", "Y"] | None = None,
    ) -> None:
        """
        Save current simulation data to a data file (CSV, YAML or HDF).

        In order to save the content of the current object, individual domains are
        converted to `SolutionArray` objects and saved using the `~SolutionArray.save`
        method. For HDF and YAML output, all domains are written to a single container
        file with shared header information. Simulation settings of individual domains
        are preserved as meta data of the corresponding `SolutionArray` objects.
        For CSV files, only state and auxiliary data of the main 1D domain are saved.

        The complete state of the current object can be restored from HDF and YAML
        container files using the `restore` method, while individual domains can be
        loaded using `SolutionArray.restore` for further analysis. While CSV files do
        not contain complete information, they can be used for setting initial states
        of individual simulation objects (example: `~FreeFlame.set_initial_guess`).

        :param filename:
            Name of output file (CSV, YAML or HDF)
        :param name:
            Identifier of storage location within the container file; this node/group
            contains header information and multiple subgroups holding domain-specific
            `SolutionArray` data (YAML/HDF only).
        :param description:
            Custom comment describing the dataset to be stored (YAML/HDF only).
        :param overwrite:
            Force overwrite if file and/or data entry exists; optional (default=`False`)
        :param compression:
            Compression level (0-9); optional (default=0; HDF only)
        :param basis:
            Output mass (``Y``/``mass``) or mole (``X``/``mole``) fractions;
            if not specified (`None`), the native basis of the underlying `ThermoPhase`
            manager is used.

        >>> s.save(filename='save.yaml', name='energy_off',
        ...        description='solution with energy eqn. disabled')

        .. versionchanged:: 3.0
            Argument loglevel is no longer supported
        """
        if loglevel is not None:
            warnings.warn("Sim1D.save: Argument 'loglevel' is deprecated and will be "
                "ignored.", DeprecationWarning)
        self.sim.save(stringify(str(filename)), stringify(name),
                      stringify(description), overwrite, compression, stringify(basis))

    def restore(
        self,
        filename: _Path | str = 'soln.yaml',
        name: str = 'solution',
        loglevel: _LogLevel | None = None,
    ) -> _RestoreMetadata:
        """Retrieve data and settings from a previously saved simulation.

        This method restores a simulation object from YAML or HDF data previously saved
        using the `save` method.

        :param filename:
            Name of container file (YAML or HDF)
        :param name:
            Identifier of location within the container file; this node/group contains
            header information and subgroups with domain-specific `SolutionArray` data
        :param loglevel:
            Amount of logging information to display while restoring,
            from 0 (disabled) to 2 (most verbose).
        :return:
            Dictionary containing header information

        >>> s.restore(filename='save.yaml', name='energy_off')

        .. versionchanged:: 3.0
            Implemented return value for meta data; loglevel is no longer supported
        """
        if loglevel is not None:
            warnings.warn("Sim1D.restore: Argument 'loglevel' is deprecated and will be"
                 " ignored.", DeprecationWarning)
        header: CxxAnyMap
        header = self.sim.restore(stringify(str(filename)), stringify(name))
        self._initialized = True
        return anymap_to_py(header)

    def restore_time_stepping_solution(self) -> None:
        """
        Set the current solution vector to the last successful time-stepping
        solution. This can be used to examine the solver progress after a failed
        integration.
        """
        self.sim.restoreTimeSteppingSolution()

    def restore_steady_solution(self) -> None:
        """
        Set the current solution vector to the last successful steady-state
        solution. This can be used to examine the solver progress after a
        failure during grid refinement.
        """
        self.sim.restoreSteadySolution()

    def show_stats(self, print_time: bool = True) -> None:
        """
        Show the statistics for the last solution.

        If invoked with no arguments or with a non-zero argument, the timing
        statistics will be printed. Otherwise, the timing will not be printed.
        """
        self.sim.writeStats(print_time)

    def clear_stats(self) -> None:
        """
        Clear solver statistics.
        """
        self.sim.clearStats()

    def solve_adjoint(
        self,
        perturb: "_Callable[[Sim1D, int, float], None]",
        n_params: int,
        dgdx: _Array,
        g: "_Callable[[Sim1D], float] | None" = None,
        dp: float = 1e-5,
    ) -> _Array:
        """
        Find the sensitivities of an objective function using an adjoint method.

        For an objective function :math:`g(x, p)` where :math:`x` is the state
        vector of the system and :math:`p` is a vector of parameters, this
        computes the vector of sensitivities :math:`dg/dp`. This assumes that
        the system of equations has already been solved to find :math:`x`.

        :param perturb:
            A function with the signature ``perturb(sim, i, dp)`` which
            perturbs parameter ``i`` by a relative factor of ``dp``. To
            perturb a reaction rate constant, this function could be defined
            as::

                def perturb(sim, i, dp):
                    sim.gas.set_multiplier(1+dp, i)

            Calling ``perturb(sim, i, 0)`` should restore that parameter to its
            default value.
        :param n_params:
            The length of the vector of sensitivity parameters
        :param dgdx:
            The vector of partial derivatives of the function :math:`g(x, p)`
            with respect to the system state :math:`x`.
        :param g:
            A function with the signature ``value = g(sim)`` which computes the
            value of :math:`g(x,p)` at the current system state. This is used to
            compute :math:`\partial g/\partial p`. If this is identically zero
            (that is, :math:`g` is independent of :math:`p`) then this argument may
            be omitted.
        :param dp:
            A relative value by which to perturb each parameter
        """
        n_vars = self.sim.size()
        L = np.empty(n_vars)
        gg = np.ascontiguousarray(dgdx, dtype=np.double)
        cL: cython.double[::1] = L
        cgg: cython.double[::1] = gg

        self.sim.solveAdjoint(span[const_double](cython.address(cgg[0]),
                                                 cython.cast(cython.size_t, gg.size)),
                              span[double](cython.address(cL[0]),
                                           cython.cast(cython.size_t, L.size)))

        dgdp = np.empty(n_params)
        dfdp = np.empty((n_vars, n_params))
        fplus = np.empty(n_vars)
        fminus = np.empty(n_vars)
        cfplus: cython.double[::1] = fplus
        cfminus: cython.double[::1] = fminus
        gplus = gminus = 0

        for i in range(n_params):
            perturb(self, i, dp)
            if g:
                gplus = g(self)
            self.sim.getResidual(0, span[double](cython.address(cfplus[0]),
                                                  cython.cast(cython.size_t, fplus.size)))

            perturb(self, i, -dp)
            if g:
                gminus = g(self)
            self.sim.getResidual(0, span[double](cython.address(cfminus[0]),
                                                  cython.cast(cython.size_t, fminus.size)))

            perturb(self, i, 0)
            dgdp[i] = (gplus - gminus)/(2*dp)
            dfdp[:,i] = (fplus - fminus) / (2*dp)

        return dgdp - np.dot(L, dfdp)

    @property
    def solver_stats(self) -> dict[str, float]:
        """
        Solver statistics from the most recent solve, as a dict of per-grid arrays.
        Keys: ``grid_points``, ``steps``, ``residual_evals``, ``residual_time``,
        ``jacobian_evals``, ``jacobian_time``, ``factorizations``, ``factor_time``,
        ``linear_solves``, ``solve_time``, and ``total_time``. Time values are
        wall-clock seconds.

        .. versionadded:: 4.0
        """
        return anymap_to_py(self.sim.solverStats())

    @property
    def grid_size_stats(self) -> list[int]:
        """
        Return total grid size in each call to solve().

        .. deprecated:: 4.0
            Use `solver_stats` instead.
        """
        warnings.warn("Sim1D.grid_size_stats is deprecated and will be "
            "removed after Cantera 4.0. Use Sim1D.solver_stats instead.",
            DeprecationWarning)
        return anymap_to_py(self.sim.solverStats())["grid_points"]

    @property
    def jacobian_time_stats(self) -> list[int]:
        """
        Return CPU time spent evaluating Jacobians in each call to solve().

        .. deprecated:: 4.0
            Use `solver_stats` instead.
        """
        warnings.warn("Sim1D.jacobian_time_stats is deprecated and will be "
            "removed after Cantera 4.0. Use Sim1D.solver_stats instead.",
            DeprecationWarning)
        return anymap_to_py(self.sim.solverStats())["jacobian_time"]

    @property
    def jacobian_count_stats(self) -> list[int]:
        """
        Return number of Jacobian evaluations made in each call to solve().

        .. deprecated:: 4.0
            Use `solver_stats` instead.
        """
        warnings.warn("Sim1D.jacobian_count_stats is deprecated and will be "
            "removed after Cantera 4.0. Use Sim1D.solver_stats instead.",
            DeprecationWarning)
        return anymap_to_py(self.sim.solverStats())["jacobian_evals"]

    @property
    def eval_time_stats(self) -> list[float]:
        """
        Return CPU time spent on non-Jacobian function evaluations in each call
        to solve().

        .. deprecated:: 4.0
            Use `solver_stats` instead.
        """
        warnings.warn("Sim1D.eval_time_stats is deprecated and will be "
            "removed after Cantera 4.0. Use Sim1D.solver_stats instead.",
            DeprecationWarning)
        return anymap_to_py(self.sim.solverStats())["residual_time"]

    @property
    def eval_count_stats(self) -> list[int]:
        """
        Return number of non-Jacobian function evaluations made in each call to
        solve().

        .. deprecated:: 4.0
            Use `solver_stats` instead.
        """
        warnings.warn("Sim1D.eval_count_stats is deprecated and will be "
            "removed after Cantera 4.0. Use Sim1D.solver_stats instead.",
            DeprecationWarning)
        return anymap_to_py(self.sim.solverStats())["residual_evals"]

    @property
    def time_step_stats(self) -> list[int]:
        """
        Return number of time steps taken in each call to solve().

        .. deprecated:: 4.0
            Use `solver_stats` instead.
        """
        warnings.warn("Sim1D.time_step_stats is deprecated and will be "
            "removed after Cantera 4.0. Use Sim1D.solver_stats instead.",
            DeprecationWarning)
        return anymap_to_py(self.sim.solverStats())["steps"]

    def set_max_grid_points(self, domain: Domain1D | str | int, npmax: int) -> None:
        """ Set the maximum number of grid points in the specified domain. """
        idom = self.domain_index(domain)
        self.sim.setMaxGridPoints(idom, npmax)

    def get_max_grid_points(self, domain: Domain1D | str | int) -> int:
        """ Get the maximum number of grid points in the specified domain. """
        idom = self.domain_index(domain)
        return self.sim.maxGridPoints(idom)
