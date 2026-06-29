# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# distutils: language = c++
# cython: language_level=3

import warnings
import numbers as _numbers
import numpy as np

from collections.abc import Callable as _Callable, Iterable as _Iterable, \
    Sequence as _Sequence
from typing import (
    Any as _Any,
    ClassVar as _ClassVar,
    Literal as _Literal,
    TypeAlias as _TypeAlias,
    overload as _overload,
    TYPE_CHECKING,
)
from typing_extensions import Never as _Never

import cython
import cython.cimports.numpy as cnp  # Required: triggers import_array() for numpy C-API
from cython.cimports.cython import view

from cython.cimports.cantera._utils import (
    pystr, stringify, comp_map, py_to_anymap, anymap_to_py)
from cython.cimports.cantera.delegator import CxxDelegatorPtr, assign_delegates
from cython.cimports.cantera.solutionbase import _SolutionBase, _wrap_Solution

from ._types import Array as _Array, ArrayLike as _ArrayLike, LogLevel as _LogLevel
from .func1 import _Func1Like

if TYPE_CHECKING:
    from graphviz import Digraph as _Digraph
    from .composite import Solution as _Solution
    from .kinetics import _DerivativeSettings

from .drawnetwork import *

# `anymap_to_py` may return an `AnyMap` (a `dict` subclass) rather than a plain
# `dict`; an inline `dict[str, int]` return annotation is coerced by Cython 3 and
# rejects that subclass instance at runtime. Route through a TypeAlias (not coerced),
# matching the `_onedim._RestoreMetadata` precedent.
_SolverStats: _TypeAlias = dict[str, int]


@cython.cclass
class ReactorBase:
    """
    Common base class for reactors and reservoirs.
    """
    reactor_type: _ClassVar[str] = "none"
    node_attr: dict[str, str] | None

    def __cinit__(self, phase: _SolutionBase, *args, name="(none)", clone=True,
                  **kwargs):
        if not isinstance(self, ReactorSurface):
            self._rbase = newReactorBase(stringify(self.reactor_type),
                                        phase._base, clone, stringify(name))
            self.rbase = self._rbase.get()

    def __init__(self, phase: _SolutionBase | None = None, *args: _Any,
                 clone: bool | None = None, name: str = "(none)",
                 volume: float | None = None,
                 node_attr: dict[str, str] | None = None) -> None:
        self._inlets = []
        self._outlets = []
        self._walls = []
        self._surfaces = []
        self._phase = _wrap_Solution(self.rbase.phase())

        if volume is not None:
            self.volume = volume

        self.node_attr = node_attr or {}

    @property
    def type(self) -> str:
        """The type of the reactor."""
        return pystr(self.rbase.type())

    @property
    def name(self) -> str:
        """The name of the reactor."""
        return pystr(self.rbase.name())

    @name.setter
    def name(self, name: str) -> None:
        self.rbase.setName(stringify(name))

    def component_index(self, name: str) -> int:
        """
        Returns the index of the component named ``name`` in the system. This determines
        the index of the component in the vector of sensitivity coefficients. ``name``
        is either a species name or the name of a reactor state variable, for example
        ``'int_energy'`` or ``'temperature'``, depending on the reactor's equations.
        """
        k = self.rbase.componentIndex(stringify(name))
        if k == -1:
            raise IndexError('No such component: {!r}'.format(name))
        return k

    def component_name(self, i: cython.int) -> str:
        """
        Returns the name of the component with index ``i`` within the array of
        variables returned by `get_state`. This is the inverse of
        `component_index`.
        """
        return pystr(self.rbase.componentName(i))

    @property
    def n_vars(self) -> int:
        """
        The number of state variables in the reactor.
        Equal to:

        `Reactor` and `IdealGasReactor`: `n_species` + 3 (mass, volume,
        internal energy or temperature).

        `ConstPressureReactor` and `IdealGasConstPressureReactor`:
        `n_species` + 2 (mass, enthalpy or temperature).
        """
        return self.rbase.neq()

    def get_state(self) -> _Array:
        """
        Get the state vector of the reactor.

        The order of the variables (that is, rows) is:

        `Reactor` or `IdealGasReactor`:

        - 0  - mass
        - 1  - volume
        - 2  - internal energy or temperature
        - 3+ - mass fractions of the species

        `ConstPressureReactor` or `IdealGasConstPressureReactor`:

        - 0  - mass
        - 1  - enthalpy or temperature
        - 2+ - mass fractions of the species

        You can use the function `component_index` to determine the location of
        a specific component from its name, or `component_name` to determine the
        name from the index.
        """
        if not self.n_vars:
            raise CanteraError('Reactor empty or network not initialized.')
        y = np.zeros(self.n_vars)
        cy: cython.double[::1] = y
        self.rbase.getState(span[double](cython.address(cy[0]),
                                        cython.cast(cython.size_t, y.size)))
        return y

    def get_state_dae(self) -> tuple[_Array, _Array]:
        """
        Get the state vector and its time derivative for reactors formulated as DAEs
        (namely, `FlowReactor`).

        .. versionadded:: 4.0
        """
        y = np.zeros(self.n_vars)
        yp = np.zeros(self.n_vars)
        cy: cython.double[::1] = y
        cyp: cython.double[::1] = yp
        self.rbase.getStateDae(
            span[double](cython.address(cy[0]), cython.cast(cython.size_t, y.size)),
            span[double](cython.address(cyp[0]), cython.cast(cython.size_t, yp.size)))
        return y, yp

    @property
    def atol(self) -> _Array | None:
        """
        Absolute tolerances for this reactor's local state variables.

        If set, this array overrides the network-level scalar absolute tolerance
        for this reactor. Entries are ordered according to `component_name`.
        Setting this property to ``None`` clears the local override.

        .. versionadded:: 4.0
        """
        atol = np.empty(self.n_vars)
        catol: cython.double[::1] = atol
        if self.rbase.getAbsoluteTolerances(
                span[double](cython.address(catol[0]),
                             cython.cast(cython.size_t, atol.size))):
            return atol
        return None

    @atol.setter
    def atol(self, values: _ArrayLike | None) -> None:
        if values is None:
            self.rbase.clearAbsoluteTolerances()
            return
        if len(values) != self.n_vars:
            raise ValueError('array must be of length n_vars')
        data = np.ascontiguousarray(values, dtype=np.double)
        cdata: cython.double[::1] = data
        self.rbase.setAbsoluteTolerances(
            span[const_double](cython.address(cdata[0]),
                               cython.cast(cython.size_t, data.size)))

    def syncState(self) -> None:
        """
        Set the state of the Reactor to match that of the associated
        `ThermoPhase` object. After calling syncState(), call
        ReactorNet.reinitialize() before further integration.

        .. deprecated:: 4.0
           Manual synchronization of reactor state is no longer required. Call
           `ReactorNet.reinitialize` directly to indicate a change in state that
           requires integrator reinitialization.
        """
        self.rbase.syncState()

    @property
    def phase(self) -> "_Solution":
        """
        The `Solution` object representing the reactor's contents.

        .. versionchanged:: 3.2
           Renamed from ``thermo``.
        """
        return self._phase

    @property
    def volume(self) -> float:
        """The volume [m³] of the reactor."""
        return self.rbase.volume()

    @volume.setter
    def volume(self, volume: cython.double) -> None:
        self.rbase.setInitialVolume(volume)

    @property
    def T(self) -> float:
        """The temperature [K] of the reactor's contents."""
        return self.phase.T

    @property
    def density(self) -> float:
        """The density [kg/m³ or kmol/m³] of the reactor's contents."""
        return self.phase.density

    @property
    def mass(self) -> float:
        """The mass of the reactor's contents."""
        return self.phase.density_mass * self.volume

    @property
    def Y(self) -> _Array:
        """The mass fractions of the reactor's contents."""
        return self.phase.Y

    def add_sensitivity_reaction(self, m: cython.int) -> None:
        """
        Specifies that the sensitivity of the state variables with respect to
        reaction ``m`` should be computed. ``m`` is the 0-based reaction index.
        The reactor must be part of a network first. Specifying the same
        reaction more than one time raises an exception.
        """
        self.rbase.addSensitivityReaction(m)

    # Flow devices & walls
    @property
    def inlets(self) -> list["FlowDevice"]:
        """List of flow devices installed as inlets to this reactor"""
        return self._inlets

    @property
    def outlets(self) -> list["FlowDevice"]:
        """List of flow devices installed as outlets to this reactor"""
        return self._outlets

    @property
    def walls(self) -> list["Wall"]:
        """List of walls installed on this reactor"""
        return self._walls

    @property
    def surfaces(self) -> list["ReactorSurface"]:
        """List of reacting surfaces installed on this reactor"""
        return self._surfaces

    def _add_inlet(self, inlet: "FlowDevice") -> None:
        """
        Store a reference to ``inlet`` to prevent it from being prematurely
        garbage collected.
        """
        self._inlets.append(inlet)

    def _add_outlet(self, outlet: "FlowDevice") -> None:
        """
        Store a reference to ``outlet`` to prevent it from being prematurely
        garbage collected.
        """
        self._outlets.append(outlet)

    def _add_wall(self, wall: "Wall") -> None:
        """
        Store a reference to ``wall`` to prevent it from being prematurely
        garbage collected.
        """
        self._walls.append(wall)

    def draw(self, graph: "_Digraph | None" = None, *,
              graph_attr: dict[str, str] | None = None,
              node_attr: dict[str, str] | None = None, print_state: bool = False,
              species: _Literal["X", "Y"] | bool | _Iterable[str] | None = None,
              species_units: _Literal["percent", "ppm"] = "percent") -> "_Digraph":
        """
        Draw as ``graphviz`` ``dot`` node. The node is added to an existing ``graph`` if
        provided. Optionally include current reactor state in the node.

        :param graph:
            ``graphviz.graphs.BaseGraph`` object to which the reactor is added.
            If not provided, a new ``DiGraph`` is created. Defaults to ``None``.
        :param graph_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the general appearance of the drawn network.
            See https://graphviz.org/docs/graph/ for a list of all usable attributes.
        :param node_attr:
            Attributes to be passed to the ``node`` method invoked to draw the reactor.
            See https://graphviz.org/docs/nodes/ for a list of all usable attributes.
        :param print_state:
            Whether state information of the reactor is printed into the node.
            Defaults to ``False``.
        :param species:
            If ``print_state`` is ``True``, define how species are to be printed.
            Options are ``'X'`` and ``'Y'`` for mole and mass fractions of all species,
            respectively, or an iterable that contains the desired species names as
            strings. Defaults to ``None``.
        :param species_units:
            Defines the units the species are displayed in as either ``"percent"`` or
            ``"ppm"``. Defaults to ``"percent"``.
        :return:
            ``graphviz.graphs.BaseGraph`` object with reactor

        .. versionadded:: 3.1
        """
        return draw_reactor(self, graph, graph_attr, node_attr, print_state, species,
                            species_units)

    def __reduce__(self) -> _Never:
        raise NotImplementedError('Reactor object is not picklable')

    def __copy__(self) -> _Never:
        raise NotImplementedError('Reactor object is not copyable')


@cython.cclass
class Reactor(ReactorBase):
    """
    A homogeneous zero-dimensional reactor. By default, they are closed
    (no inlets or outlets), have fixed volume, and have adiabatic,
    chemically-inert walls. These properties may all be changed by adding
    appropriate components such as `Wall`, `MassFlowController` and `Valve`.
    """
    reactor_type: _ClassVar[str] = "Reactor"
    group_name: str

    def __cinit__(self, *args, **kwargs):
        self.reactor = cython.cast(cython.pointer(CxxReactor), self.rbase)

    def __init__(self, phase: "_Solution", *, clone: bool | None = None,
                 name: str = "(none)", energy: _Literal["on", "off"] = "on",
                 group_name: str = "", **kwargs: _Any) -> None:
        """
        :param phase:
            A `Solution` object representing the Reactor contents
        :param clone:
            Determines whether to clone the ``phase`` object used by this reactor so
            that the internal state is independent of the original `Solution` (and any
            `Solution` objects used by other reactors in the network).
        :param name:
            Name string. If not specified, the name initially defaults to ``'(none)'``
            and changes to ``'<reactor_type>_n'`` when `Reactor` objects are installed
            within a `ReactorNet`. For the latter, ``<reactor_type>`` is the type of
            the reactor and *n* is an integer assigned in the order reactors are
            installed.
        :param energy:
            Set to ``'on'`` or ``'off'``. If set to ``'off'``, the energy
            equation is not solved, and the temperature is held at its
            initial value.
        :param node_attr:
            Attributes to be passed to the ``node`` method invoked to draw this reactor.
            See https://graphviz.org/docs/nodes/ for a list of all usable attributes.
        :param group_name:
            Group reactors of the same ``group_name`` when drawn using graphviz.

        .. versionadded:: 3.1
           Added the ``node_attr`` and ``group_name`` parameters.

        .. versionchanged:: 3.2
           Added the ``clone`` parameter with the default value ``None`` indicating that
           contents will not be cloned. Explicitly specifying ``False`` suppresses a
           warning.

        .. versionchanged:: 3.3
           Changed the default value of the ``clone`` to ``True``.

        Some examples showing how to create :class:`Reactor` objects are
        shown below.

        >>> gas = Solution('gri30.yaml')
        >>> r1 = Reactor(gas)

        Arguments may be specified using keywords:

        >>> r2 = Reactor(gas, energy='off',
        ...              name='isothermal_reactor')

        """
        super().__init__(phase, clone=clone, name=name, **kwargs)

        if energy == 'off':
            self.energy_enabled = False
        elif energy != 'on':
            raise ValueError("'energy' must be either 'on' or 'off'")

        self.group_name = group_name

    @property
    def chemistry_enabled(self) -> bool:
        """
        `True` when the reactor composition is allowed to change due to
        chemical reactions in this reactor. When this is `False`, the
        reactor composition is held constant.
        """
        return self.reactor.chemistryEnabled()

    @chemistry_enabled.setter
    def chemistry_enabled(self, value: pybool) -> None:
        self.reactor.setChemistryEnabled(value)

    @property
    def energy_enabled(self) -> bool:
        """
        `True` when the energy equation is being solved for this reactor.
        When this is `False`, the reactor temperature is held constant.
        """
        return self.reactor.energyEnabled()

    @energy_enabled.setter
    def energy_enabled(self, value: pybool) -> None:
        self.reactor.setEnergyEnabled(value)

    def add_sensitivity_species_enthalpy(self, k) -> None:
        """
        Specifies that the sensitivity of the state variables with respect to
        species ``k`` should be computed. The reactor must be part of a network
        first.
        """
        self.reactor.addSensitivitySpeciesEnthalpy(self.phase.species_index(k))

    @property
    def jacobian(self) -> _Array:
        """
        Get the local, reactor-specific Jacobian or an approximation thereof

        .. warning::

            Depending on the particular implementation, this may return an approximate
            Jacobian intended only for use in forming a preconditioner for iterative
            solvers, excluding terms that would generate a fully-dense Jacobian.

        .. warning::

            This method is an experimental part of the Cantera API and may be
            changed or removed without notice.
        """
        return get_from_sparse(self.reactor.jacobian(), self.n_vars, self.n_vars)

    @property
    def finite_difference_jacobian(self) -> _Array:
        """
        Get the reactor-specific Jacobian, calculated using a finite difference method.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        return get_from_sparse(self.reactor.finiteDifferenceJacobian(),
                               self.n_vars, self.n_vars)

    def set_advance_limit(self, name: str, limit: float | None) -> None:
        """
        Limit absolute change of component ``name`` during `ReactorNet.advance`.
        (positive ``limit`` values are considered; negative values disable a
        previously set advance limit for a solution component). Note that
        limits are disabled by default (with individual values set to -1.).
        """
        if limit is None:
            limit = -1.
        self.reactor.setAdvanceLimit(stringify(name), limit)


@cython.cclass
class MoleReactor(Reactor):
    """
    A homogeneous zero-dimensional reactor with a mole based state vector. By default,
    they are closed (no inlets or outlets), have fixed volume, and have adiabatic,
    chemically-inert walls. These properties may all be changed by adding
    appropriate components such as `Wall`, `MassFlowController` and `Valve`.

    .. versionadded:: 3.0
    """
    reactor_type = "MoleReactor"


@cython.cclass
class Reservoir(ReactorBase):
    """
    A reservoir is a reactor with a constant state. The temperature,
    pressure, and chemical composition in a reservoir never change from
    their initial values.
    """
    reactor_type = "Reservoir"


@cython.cclass
class ConstPressureReactor(Reactor):
    """A homogeneous, constant pressure, zero-dimensional reactor. The volume
    of the reactor changes as a function of time in order to keep the
    pressure constant.
    """
    reactor_type = "ConstPressureReactor"


@cython.cclass
class ConstPressureMoleReactor(Reactor):
    """A homogeneous, constant pressure, zero-dimensional reactor with a mole based
    state vector. The volume of the reactor changes as a function of time in order to
    keep the pressure constant.

    .. versionadded:: 3.0
    """
    reactor_type = "ConstPressureMoleReactor"


@cython.cclass
class IdealGasReactor(Reactor):
    """ A constant volume, zero-dimensional reactor for ideal gas mixtures. """
    reactor_type = "IdealGasReactor"


@cython.cclass
class IdealGasMoleReactor(Reactor):
    """
    A constant volume, zero-dimensional reactor for ideal gas mixtures with a mole
    based state vector
    """
    reactor_type = "IdealGasMoleReactor"


@cython.cclass
class IdealGasConstPressureReactor(Reactor):
    """
    A homogeneous, constant pressure, zero-dimensional reactor for ideal gas
    mixtures. The volume of the reactor changes as a function of time in order
    to keep the pressure constant.
    """
    reactor_type = "IdealGasConstPressureReactor"


@cython.cclass
class IdealGasConstPressureMoleReactor(Reactor):
    """
    A homogeneous, constant pressure, zero-dimensional reactor for ideal gas
    mixtures. The volume of the reactor changes as a function of time in order
    to keep the pressure constant. This reactor also uses a mole based state vector.
    """
    reactor_type = "IdealGasConstPressureMoleReactor"


@cython.cclass
class FlowReactor(Reactor):
    """
    A steady-state plug flow reactor with constant cross sectional area.
    Integration follows a fluid element along the length of the reactor.
    The reactor is assumed to be frictionless and adiabatic.
    """
    reactor_type = "FlowReactor"

    @property
    def mass_flow_rate(self) -> _Never:
        """ Mass flow rate [kg/s] """
        raise AttributeError("unreadable attribute 'mass_flow_rate'")

    @mass_flow_rate.setter
    def mass_flow_rate(self, value: cython.double) -> None:
        cython.cast(cython.pointer(CxxFlowReactor), self.reactor).setMassFlowRate(value)

    @property
    def area(self) -> float:
        """
        Get/set the area of the reactor [m²].

        When the area is changed, the flow speed is scaled to keep the total mass flow
        rate constant.
        """
        return cython.cast(cython.pointer(CxxFlowReactor), self.reactor).area()

    @area.setter
    def area(self, area: float) -> None:
        cython.cast(cython.pointer(CxxFlowReactor), self.reactor).setArea(area)

    @property
    def speed(self) -> float:
        """ Speed [m/s] of the flow in the reactor at the current position """
        return cython.cast(cython.pointer(CxxFlowReactor), self.reactor).speed()


@cython.cclass
class ExtensibleReactor(Reactor):
    """
    A base class for a reactor with delegated methods where the base
    functionality corresponds to the `Reactor` class.

    The ``__init__`` method of the derived class should allocate and size any
    internal variables and set the total number of state variables associated with this
    reactor, `n_vars` (if it is different from the base class).

    The following methods of the C++ :ct:`Reactor` class can be modified by a
    Python class which inherits from this class. For each method, the name below
    should be prefixed with ``before_``, ``after_``, or ``replace_``, indicating
    whether this method should be called before, after, or instead of the
    corresponding method from the base class.

    For methods that return a value and have a ``before`` method specified, if
    that method returns a value other than ``None`` that value will be returned
    without calling the base class method; otherwise, the value from the base
    class method will be returned. For methods that return a value and have an
    ``after`` method specified, the returned value wil be the sum of the values
    from the supplied method and the base class method.

    ``initialize(self, t0: double) -> None``
        Responsible for initialization that can only be performed after connecting the
        elements of the reactor network, such as initializing attached walls.

        Called once before the start of time integration.

    ``get_state(self, y : double[:]) -> None``
        Responsible for populating the state vector ``y`` (length `n_vars`)
        with the initial state of the reactor.

    ``update_state(self, y : double[:]) -> None``
        Responsible for setting the state of the reactor object from the
        values in the state vector ``y`` (length `n_vars`)

    ``update_connected(self, update_pressure : bool) -> None``
        Responsible for storing properties which may be accessed by connected
        reactors, and for updating the mass flow rates of connected flow devices.

    ``eval(self, t : double, LHS : double[:], RHS : double[:]) -> None``
        Responsible for calculating the time derivative of the state at time ``t``
        based on the current state of the reactor. For each component ``i`` of the
        state vector, the time derivative ``dy[i]/dt`` is calculated as
        ``LHS[i] * dy[i]/dt = RHS[i]``. ``LHS`` and ``RHS`` are arrays of length
        `n_vars`.

    ``eval_walls(self, t : double) -> None``
        Responsible for calculating the net rate of volume change `expansion_rate`
        and the net rate of heat transfer `heat_rate` caused by walls connected
        to this reactor.

    ``component_name(i : int) -> string``
        Returns the name of the state vector component with index ``i``

    ``component_index(name: string) -> int``
        Returns the index of the state vector component named ``name``

    ``get_jacobian_elements(elements : list) -> None``
        Appends sparse Jacobian elements as ``(row, column, value)`` tuples. Row and
        column indices are global within the containing reactor network.
    """

    reactor_type: _ClassVar[str] = "ExtensibleReactor"

    delegatable_methods: dict[str, tuple[str, str]] = {
        'initialize': ('initialize', 'void(double)'),
        'get_state': ('getState', 'void(double*)'),
        'update_state': ('updateState', 'void(double*)'),
        'update_connected': ('updateConnected', 'void(bool)'),
        'eval': ('eval', 'void(double, double*, double*)'),
        'eval_walls': ('evalWalls', 'void(double)'),
        'component_name': ('componentName', 'string(size_t)'),
        'component_index': ('componentIndex', 'size_t(string)'),
        'get_jacobian_elements': ('getJacobianElements', 'void(SparseTriplets&)'),
    }

    surface_production_rates: _Array

    def __cinit__(self, *args, **kwargs):
        self.accessor = dynamic_cast[CxxReactorAccessorPtr](self.rbase)
        sdot: span[cython.double] = \
            dynamic_cast[CxxReactorAccessorPtr](self.rbase).surfaceProductionRates()
        # Non-owning memoryview over the C++ span data (pure-Python spelling of the
        # .pyx `<double[:sdot.size()]> sdot.data()` sized pointer cast).
        sarr: view.array = view.array(shape=(sdot.size(),),
                                      itemsize=cython.sizeof(cython.double), format="d",
                                      allocate_buffer=False)
        sarr.data = cython.cast(cython.p_char, sdot.data())
        self.surface_production_rates = sarr

    def __init__(self, *args: _Any, **kwargs: _Any) -> None:
        assign_delegates(self, dynamic_cast[CxxDelegatorPtr](self.rbase))
        super().__init__(*args, **kwargs)

    @property
    def n_vars(self) -> int:
        """
        Get/Set the number of state variables in the reactor.
        """
        return self.reactor.neq()

    @n_vars.setter
    def n_vars(self, n: int) -> None:
        self.accessor.setNEq(n)

    @property
    def expansion_rate(self) -> float:
        """
        Get/Set the net rate of volume change (for example, from moving walls) [m³/s]

        .. versionadded:: 3.0
        """
        return self.accessor.expansionRate()

    @expansion_rate.setter
    def expansion_rate(self, vdot: float) -> None:
        self.accessor.setExpansionRate(vdot)

    @property
    def heat_rate(self) -> float:
        """
        Get/Set the net heat transfer rate (for example, through walls) [W]

        .. versionadded:: 3.0
        """
        return self.accessor.heatRate()

    @heat_rate.setter
    def heat_rate(self, qdot: float) -> None:
        self.accessor.setHeatRate(qdot)


@cython.cclass
class ExtensibleIdealGasReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasReactor` class.
    """
    reactor_type = "ExtensibleIdealGasReactor"


@cython.cclass
class ExtensibleConstPressureReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `ConstPressureReactor` class.
    """
    reactor_type = "ExtensibleConstPressureReactor"


@cython.cclass
class ExtensibleIdealGasConstPressureReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasConstPressureReactor` class.
    """
    reactor_type = "ExtensibleIdealGasConstPressureReactor"


@cython.cclass
class ExtensibleMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `MoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleMoleReactor"


@cython.cclass
class ExtensibleIdealGasMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasMoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleIdealGasMoleReactor"


@cython.cclass
class ExtensibleConstPressureMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `ConstPressureMoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleConstPressureMoleReactor"


@cython.cclass
class ExtensibleIdealGasConstPressureMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasConstPressureMoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleIdealGasConstPressureMoleReactor"


@cython.cclass
class ReactorSurface(ReactorBase):
    """
    Represents a reacting surface in contact with the contents of one or more reactors.

    :param phase:
        The `Interface` object representing reactions on this surface.
    :param r:
        A `Reactor` or list of `Reactor` objects that this surface is adjacent to.
    :param clone:
        Determines whether to clone the ``phase`` object used by this reactor so
        that the internal state is independent of the original `Solution` (and any
        `Solution` objects used by other reactors in the network).
    :param name:
        Name string. If not specified, the name initially defaults to ``'(none)'`` and
        changes to ``'ReactorSurface_n'`` when when associated `Reactor` objects are
        installed within a `ReactorNet`. For the latter, *n* is an integer assigned in
        the order reactor surfaces are detected.
    :param A:
        The area of the reacting surface [m²].
    :param node_attr:
        Attributes to be passed to the ``node`` method invoked to draw this surface.
        See https://graphviz.org/docs/nodes/ for a list of all usable attributes.

    .. versionadded:: 3.1
       Added the ``node_attr`` parameter.

    .. versionchanged:: 3.2
       Handle surfaces that are adjacent to multiple reactors.

       Added the ``clone`` parameter with the default value ``None`` indicating that
       contents will not be cloned. Specifying ``False`` suppresses a warning.

    .. versionchanged:: 3.3
       Changed the default value of the ``clone`` to ``True``.
    """
    reactor_type: _ClassVar[str] = "ReactorSurface"

    def __cinit__(self, phase: _SolutionBase, r, *, clone=True, name="(none)",
                  kind=None, **kwargs):
        adj: ReactorBase
        cxx_adj: vector[shared_ptr[CxxReactorBase]]
        if isinstance(r, ReactorBase):
            adj = cython.cast(ReactorBase, r)
            adj._surfaces.append(self)
            cxx_adj.push_back(adj._rbase)
            self._reactors = [r]
        elif hasattr(r, "__len__"):
            self._reactors = r
            for ri in r:
                adj = cython.cast(Reactor, ri)
                adj._surfaces.append(self)
                cxx_adj.push_back(adj._rbase)
        else:
            raise TypeError("Parameter 'r' should be a ReactorBase object or a list "
                            "of ReactorBase objects.")

        adj_span: span[shared_ptr[CxxReactorBase]] = \
            span[shared_ptr[CxxReactorBase]](cxx_adj)

        if kind is None and self.reactor_type != "ReactorSurface":
            kind = self.reactor_type

        if kind is not None:
            self._rbase = CxxNewReactorSurface(stringify(kind), phase._base, adj_span,
                                               clone, stringify(name))
        else:
            self._rbase = CxxNewReactorSurface(phase._base, adj_span, clone,
                                               stringify(name))

        self.rbase = self._rbase.get()
        self.surface = cython.cast(cython.pointer(CxxReactorSurface), self.rbase)

    def __init__(self, phase: _SolutionBase | None = None, r: "Reactor | None" = None,
                 *, kind=None, clone: bool | None = None, name: str = "(none)",
                 A: float | None = None,
                 node_attr: dict[str, str] | None = None) -> None:
        super().__init__(phase, name=name)

        if A is not None:
            self.area = A
        self.node_attr = node_attr or {'shape': 'underline'}

    @property
    def area(self) -> float:
        """Area on which reactions can occur [m²]."""
        return self.surface.area()

    @area.setter
    def area(self, A: float) -> None:
        self.surface.setArea(A)

    @property
    def coverages(self) -> _Array:
        """
        The fraction of sites covered by each surface species.
        """
        if self._phase is None:
            raise CanteraError('No kinetics manager present')
        return self._phase.coverages

    @coverages.setter
    def coverages(self, coverages: _Array) -> None:
        if self._phase is None:
            raise CanteraError("Can't set coverages before assigning kinetics manager.")

        if isinstance(coverages, (dict, str, bytes)):
            self.surface.setCoverages(comp_map(coverages))
            return

        if len(coverages) != self._phase.n_species:
            raise ValueError('Incorrect number of site coverages specified')
        data = np.ascontiguousarray(coverages, dtype=np.double)
        cdata: cython.double[::1] = data
        self.surface.setCoverages(span[double](cython.address(cdata[0]),
                                              cython.cast(cython.size_t, data.size)))

    @property
    def reactor(self) -> "Reactor":
        """
        Return the `Reactor` object the surface is connected to.

        .. versionadded:: 3.1
        """
        if len(self._reactors) > 1:
            warnings.warn("ReactorSurface.reactor: Call is ambiguous because surface is"
                " linked to multiple reactors. Use 'ReactorSurface.reactors' instead.",
                UserWarning)

        return self._reactors[0]

    @property
    def reactors(self) -> list["Reactor"]:
        """
        A list of of `Reactor` objects containing phases that participate in reactions
        on this surface.

        .. versionadded:: 3.2
        """
        return self._reactors

    def draw(self, graph: "_Digraph | None" = None, *,
              graph_attr: dict[str, str] | None = None,
              node_attr: dict[str, str] | None = None,
              surface_edge_attr: dict[str, str] | None = None,
              print_state: bool = False,
              species: _Literal["X", "Y"] | bool | _Iterable[str] | None = None,
              species_units: _Literal["percent", "ppm"] = "percent") -> "_Digraph":
        """
        Draw the surface as a ``graphviz`` ``dot`` node connected to its reactor.
        The node is added to an existing ``graph`` if provided.
        Optionally include current reactor state in the reactor node.

        :param graph:
            ``graphviz.graphs.BaseGraph`` object to which the reactor is added.
            If not provided, a new ``DiGraph`` is created. Defaults to ``None``.
        :param graph_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the general appearance of the drawn network.
            See https://graphviz.org/docs/graph/ for a list of all usable attributes.
        :param node_attr:
            Attributes to be passed to the ``node`` method invoked to draw the reactor.
            See https://graphviz.org/docs/nodes/ for a list of all usable attributes.
        :param surface_edge_attr:
            Attributes to be passed to the ``edge`` method invoked to draw the
            connection between the surface and its reactor.
            See https://graphviz.org/docs/edges/ for a list of all usable attributes.
        :param print_state:
            Whether state information of the reactor is printed into its node.
            Defaults to ``False``
        :param species:
            If ``print_state`` is ``True``, define how species are to be printed.
            Options are ``'X'`` and ``'Y'`` for mole and mass fractions of all species,
            respectively, or an iterable that contains the desired species names as
            strings. Defaults to ``None``.
        :param species_units:
            Defines the units the species are displayed in as either ``"percent"`` or
            ``"ppm"``. Defaults to ``"percent"``.
        :return:
            ``graphviz.graphs.BaseGraph`` object with surface and connected
            reactor.

        .. versionadded:: 3.1
        """
        return draw_surface(self, graph, graph_attr, node_attr, surface_edge_attr,
                            print_state, species, species_units)


@cython.cclass
class FlowReactorSurface(ReactorSurface):
    reactor_type: _ClassVar[str] = "FlowReactorSurface"

    @property
    def initial_atol(self) -> float:
        """
        Get/Set the steady-state tolerances used to determine the initial surface
        species coverages.
        """
        return cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).initialAtol()

    @initial_atol.setter
    def initial_atol(self, atol: float) -> None:
        cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).setInitialAtol(atol)

    @property
    def initial_rtol(self) -> float:
        """
        Get/Set the steady-state tolerances used to determine the initial surface
        species coverages.
        """
        return cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).initialRtol()

    @initial_rtol.setter
    def initial_rtol(self, rtol: float) -> None:
        cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).setInitialRtol(rtol)

    @property
    def initial_max_steps(self) -> int:
        """
        Get/Set the maximum number of integrator steps used to determine the initial
        surface species coverages.
        """
        return cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).initialMaxSteps()

    @initial_max_steps.setter
    def initial_max_steps(self, nsteps: int) -> None:
        cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).setInitialMaxSteps(nsteps)

    @property
    def initial_max_error_failures(self) -> int:
        """
        Get/Set the maximum number of integrator error failures allowed when determining
        the initial surface species coverages.
        """
        return cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).initialMaxErrorFailures()

    @initial_max_error_failures.setter
    def initial_max_error_failures(self, nsteps: int) -> None:
        cython.cast(cython.pointer(CxxFlowReactorSurface), self.surface).setInitialMaxErrorFailures(nsteps)


@cython.cclass
class ExtensibleReactorSurface(ReactorSurface):
    """
    A base class for a reactor surface with delegated methods where the base
    functionality corresponds to the `ReactorSurface` class.

    Methods of the base class can be modified in the same way as in class
    `ExtensibleReactor`. The methods that can be modified are::

    - ``initialize(self, t0: double) -> None``
    - ``get_state(self, y : double[:]) -> None``
    - ``update_state(self, y : double[:]) -> None``
    - ``eval(self, t : double, LHS : double[:], RHS : double[:]) -> None``
    - ``component_name(i : int) -> string``
    - ``component_index(name: string) -> int``
    """

    reactor_type: _ClassVar[str] = "ExtensibleReactorSurface"

    delegatable_methods: dict[str, tuple[str, str]] = {
        'initialize': ('initialize', 'void(double)'),
        'get_state': ('getState', 'void(double*)'),
        'update_state': ('updateState', 'void(double*)'),
        'eval': ('eval', 'void(double, double*, double*)'),
        'component_name': ('componentName', 'string(size_t)'),
        'component_index': ('componentIndex', 'size_t(string)'),
    }

    surface_production_rates: _Array

    def __init__(self, *args: _Any, **kwargs: _Any) -> None:
        assign_delegates(self, dynamic_cast[CxxDelegatorPtr](self.rbase))
        sdot: span[cython.double] = \
            dynamic_cast[CxxReactorAccessorPtr](self.rbase).surfaceProductionRates()
        # Non-owning memoryview over the C++ span data (pure-Python spelling of the
        # .pyx `<double[:sdot.size()]> sdot.data()` sized pointer cast).
        sarr: view.array = view.array(shape=(sdot.size(),),
                                      itemsize=cython.sizeof(cython.double), format="d",
                                      allocate_buffer=False)
        sarr.data = cython.cast(cython.p_char, sdot.data())
        self.surface_production_rates = sarr
        super().__init__(*args, **kwargs)

    @property
    def n_vars(self) -> int:
        """
        Get/Set the number of state variables in the reactor.
        """
        return self.reactor.neq()

    @n_vars.setter
    def n_vars(self, n) -> None:
        dynamic_cast[CxxReactorAccessorPtr](self.rbase).setNEq(n)


@cython.cclass
class ExtensibleMoleReactorSurface(ExtensibleReactorSurface):
    """
    A variant of `ExtensibleReactorSurface` where the base behavior corresponds to the
    :ct:`MoleReactorSurface` class.
    """
    reactor_type: _ClassVar[str] = "ExtensibleMoleReactorSurface"


@cython.cclass
class ConnectorNode:
    """
    Common base class for walls and flow devices.
    """
    node_type: _ClassVar[str] = "none"
    edge_attr: dict[str, str]

    def __cinit__(self, left: ReactorBase = None, right: ReactorBase = None, *,
                  upstream: ReactorBase = None, downstream: ReactorBase = None,
                  name="(none)", **kwargs):
        # ensure that both naming conventions (Wall and FlowDevice) are covered
        r0: ReactorBase = left or upstream
        r1: ReactorBase = right or downstream
        if isinstance(r0, ReactorBase) and isinstance(r1, ReactorBase):
            self._node = newConnectorNode(stringify(self.node_type),
                                         r0._rbase, r1._rbase, stringify(name))
            self.node = self._node.get()
            return
        raise TypeError(f"Invalid reactor types: {r0} and {r1}.")

    @property
    def type(self) -> str:
        """The type of the connector."""
        return pystr(self.node.type())

    @property
    def name(self) -> str:
        """The name of the connector."""
        return pystr(self.node.name())

    @name.setter
    def name(self, name: str) -> None:
        self.node.setName(stringify(name))

    def __reduce__(self) -> _Never:
        raise NotImplementedError('Reactor object is not picklable')

    def __copy__(self) -> _Never:
        raise NotImplementedError('Reactor object is not copyable')


@cython.cclass
class WallBase(ConnectorNode):
    """
    Common base class for walls.
    """
    def __cinit__(self, *args, **kwargs):
        self.wall = cython.cast(cython.pointer(CxxWallBase), self.node)

    def __init__(self, left: ReactorBase, right: ReactorBase, *, name: str = "(none)",
                 A: float | None = None, K: float | None = None, U: float | None = None,
                 Q: _Callable[[float], float] | None = None,
                 velocity: _Callable[[float], float] | None = None,
                 edge_attr: dict[str, str] | None = None) -> None:
        """
        :param left:
            Reactor or reservoir on the left. Required.
        :param right:
            Reactor or reservoir on the right. Required.
        :param name:
            Name string. If not specified, the name initially defaults to ``'(none)'``
            and changes to ``'<wall_type>_n'`` when when associated `Reactor` objects
            are installed within a `ReactorNet`. For the latter, ``<wall_type>`` is
            the type of the wall and *n* is an integer assigned in the order walls are
            detected.
        :param A:
            Wall area [m²]. Defaults to 1.0 m².
        :param K:
            Wall expansion rate parameter [m/s/Pa]. Defaults to 0.0.
        :param U:
            Overall heat transfer coefficient [W/m²/K]. Defaults to 0.0
            (adiabatic wall).
        :param Q:
            Heat flux function :math:`q_0(t)` [W/m²]. Optional. Default:
            :math:`q_0(t) = 0.0`.
        :param velocity:
            Wall velocity function :math:`v_0(t)` [m/s].
            Default: :math:`v_0(t) = 0.0`.
        :param edge_attr:
            Attributes like ``style`` when drawn as a connecting edge using
            graphviz's dot language. Default is ``{}``.

        .. versionadded:: 3.1
           Added the ``edge_attr`` parameter.

        """
        self._velocity_func = None
        self._heat_flux_func = None

        if A is not None:
            self.area = A
        if K is not None:
            self.expansion_rate_coeff = K
        if U is not None:
            self.heat_transfer_coeff = U
        if Q is not None:
            self.heat_flux = Q
        if velocity is not None:
            self.velocity = velocity
        self.edge_attr = edge_attr or {}

        left._add_wall(self)
        right._add_wall(self)
        # Keep references to prevent premature garbage collection
        self._left_reactor = left
        self._right_reactor = right

    @property
    def area(self) -> float:
        """The wall area [m²]."""
        return self.wall.area()

    @area.setter
    def area(self, value: cython.double) -> None:
        self.wall.setArea(value)

    @property
    def left_reactor(self) -> ReactorBase:
        """
        Return the `Reactor` or `Reservoir` object left of the wall.

        .. versionadded:: 3.1
        """
        return self._left_reactor

    @property
    def right_reactor(self) -> ReactorBase:
        """
        Return the `Reactor` or `Reservoir` object right of the wall.

        .. versionadded:: 3.1
        """
        return self._right_reactor

    @property
    def expansion_rate(self) -> float:
        """
        Get the rate of volumetric change [m³/s] associated with the wall at the
        current reactor network time. A positive value corresponds to the left-hand
        reactor volume increasing, and the right-hand reactor volume decreasing.

        .. versionadded:: 3.0
        """
        return self.wall.expansionRate()

    @property
    def heat_rate(self) -> float:
        """
        Get the total heat flux [W] through the wall  at the current reactor network
        time. A positive value corresponds to heat flowing from the left-hand reactor
        to the right-hand one.

        .. versionadded:: 3.0
        """
        return self.wall.heatRate()

    def draw(self, graph: "_Digraph | None" = None, *,
              graph_attr: dict[str, str] | None = None,
              node_attr: dict[str, str] | None = None,
              edge_attr: dict[str, str] | None = None,
              moving_wall_edge_attr: dict[str, str] | None = None,
              show_wall_velocity: bool = True) -> "_Digraph":
        """
        Draw as connection between left and right reactor or reservoir using
        ``graphviz``.

        :param graph:
            ``graphviz.graphs.BaseGraph`` object to which the connection is added.
            If not provided, a new ``DiGraph`` is created. Defaults to ``None``
        :param graph_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the general appearance of the drawn network.
            Has no effect if existing ``graph`` is provided.
            See https://graphviz.org/docs/graph/ for a list of all usable attributes.
        :param node_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the default appearance of any ``node`` (reactors, reservoirs).
            Has no effect if existing ``graph`` is provided.
            See https://graphviz.org/docs/nodes/ for a list of all usable attributes.
        :param edge_attr:
            Attributes to be passed to the ``edge`` method invoked to draw this wall
            connection.
            Default is ``{"color": "red", "style": "dashed"}``.
            See https://graphviz.org/docs/edges/ for a list of all usable attributes.
        :param moving_wall_edge_attr:
            Same as ``edge_attr`` but only applied to edges representing wall movement.
            Default is ``{"arrowtail": "icurveteecurve", "dir": "both",
            "style": "dotted", "arrowhead": "icurveteecurve"}``.
        :param show_wall_velocity:
            If ``True``, wall movement will be indicated by additional arrows with the
            corresponding wall velocity as a label.
        :return:
            A ``graphviz.graphs.BaseGraph`` object depicting the connection.

        .. versionadded:: 3.1
        """
        return draw_walls([self], graph, graph_attr, node_attr, edge_attr,
                          moving_wall_edge_attr, show_wall_velocity)


@cython.cclass
class Wall(WallBase):
    r"""
    A Wall separates two reactors, or a reactor and a reservoir. A wall has a
    finite area, may conduct or radiate heat between the two reactors on either
    side, and may move like a piston.

    Walls are stateless objects in Cantera, meaning that no differential
    equation is integrated to determine any wall property. Since it is the wall
    (piston) velocity that enters the energy equation, this means that it is
    the velocity, not the acceleration or displacement, that is specified.
    The wall velocity is computed from

    .. math:: v = K(P_{\rm left} - P_{\rm right}) + v_0(t),

    where :math:`K` is a non-negative constant, and :math:`v_0(t)` is a
    specified function of time. The velocity is positive if the wall is
    moving to the right.

    The heat flux through the wall is computed from

    .. math::  q = U(T_{\rm left} - T_{\rm right}) + \epsilon\sigma (T_{\rm left}^4 - T_{\rm right}^4) + q_0(t),

    where :math:`U` is the overall heat transfer coefficient for
    conduction/convection, and :math:`\epsilon` is the emissivity. The function
    :math:`q_0(t)` is a specified function of time. The heat flux is positive
    when heat flows from the reactor on the left to the reactor on the right.
    """
    node_type: _ClassVar[str] = "Wall"

    @property
    def expansion_rate_coeff(self) -> float:
        """
        The coefficient *K* [m/s/Pa] that determines the velocity of the wall
        as a function of the pressure difference between the adjacent reactors.
        """
        return cython.cast(cython.pointer(CxxWall), self.wall).getExpansionRateCoeff()

    @expansion_rate_coeff.setter
    def expansion_rate_coeff(self, val: cython.double) -> None:
        cython.cast(cython.pointer(CxxWall), self.wall).setExpansionRateCoeff(val)

    @property
    def heat_transfer_coeff(self) -> float:
        """The overall heat transfer coefficient [W/m²/K]."""
        return cython.cast(cython.pointer(CxxWall), self.wall).getHeatTransferCoeff()

    @heat_transfer_coeff.setter
    def heat_transfer_coeff(self, value: cython.double) -> None:
        cython.cast(cython.pointer(CxxWall), self.wall).setHeatTransferCoeff(value)

    @property
    def emissivity(self) -> float:
        """The emissivity (nondimensional)"""
        return cython.cast(cython.pointer(CxxWall), self.wall).getEmissivity()

    @emissivity.setter
    def emissivity(self, value: cython.double) -> None:
        cython.cast(cython.pointer(CxxWall), self.wall).setEmissivity(value)

    @property
    def velocity(self) -> float:
        """
        The wall velocity [m/s]. May be either set to a constant or an arbitrary
        function of time. See `Func1`.

        .. versionadded:: 3.0
        """
        return cython.cast(cython.pointer(CxxWall), self.wall).velocity()

    @velocity.setter
    def velocity(self, v: _Func1Like) -> None:
        f: Func1
        if isinstance(v, Func1):
            f = v
        else:
            f = Func1(v)

        self._velocity_func = f
        cython.cast(cython.pointer(CxxWall), self.wall).setVelocity(f._func)

    @property
    def heat_flux(self) -> float:
        """
        Heat flux [W/m²] across the wall. May be either set to a constant or
        an arbitrary function of time. See `Func1`.

        .. versionadded:: 3.0
        """
        return cython.cast(cython.pointer(CxxWall), self.wall).heatFlux()

    @heat_flux.setter
    def heat_flux(self, q: _Func1Like) -> None:
        f: Func1
        if isinstance(q, Func1):
            f = q
        else:
            f = Func1(q)

        self._heat_flux_func = f
        cython.cast(cython.pointer(CxxWall), self.wall).setHeatFlux(f._func)


@cython.cclass
class FlowDevice(ConnectorNode):
    """
    Base class for devices that allow flow between reactors.

    FlowDevice objects are assumed to be adiabatic, non-reactive, and have
    negligible internal volume, so that they are internally always in
    steady-state even if the upstream and downstream reactors are not. The
    fluid enthalpy, chemical composition, and mass flow rate are constant
    across a FlowDevice, and the pressure difference equals the difference in
    pressure between the upstream and downstream reactors.
    """
    def __cinit__(self, *args, **kwargs):
        self.dev = cython.cast(cython.pointer(CxxFlowDevice), self.node)

    def __init__(self, upstream: ReactorBase, downstream: ReactorBase, *,
                 name: str = "(none)",
                 edge_attr: dict[str, str] | None = None) -> None:
        assert self.dev != cython.NULL
        self._rate_func = None
        self.edge_attr = edge_attr or {}
        upstream._add_outlet(self)
        downstream._add_inlet(self)
        # Keep references to prevent premature garbage collection
        self._upstream = upstream
        self._downstream = downstream

    @property
    def upstream(self) -> ReactorBase:
        """
        Return the `Reactor` or `Reservoir` object upstream of the flow device.

        .. versionadded:: 3.1
        """
        return self._upstream

    @property
    def downstream(self) -> ReactorBase:
        """
        Return the `Reactor` or `Reservoir` object downstream of the flow device.

        .. versionadded:: 3.1
        """
        return self._downstream

    @property
    def mass_flow_rate(self) -> float:
        """
        Get the mass flow rate [kg/s] through this device at the current reactor
        network time.
        """
        return self.dev.massFlowRate()

    @property
    def pressure_function(self) -> float:
        r"""
        The relationship between mass flow rate and the pressure drop across a flow
        device. The mass flow rate [kg/s] is calculated given the pressure drop [Pa] and
        a coefficient set by a flow device specific function. Unless a user-defined
        pressure function is provided, the function returns the pressure difference
        across the device. The calculation of mass flow rate depends on the flow device.

        >>> f = FlowDevice(res1, reactor1)
        >>> f.pressure_function = lambda dP: dP**2

        where `FlowDevice` is either a `Valve` or `PressureController` object.

        .. versionadded:: 3.0
        """
        return self.dev.evalPressureFunction()

    @pressure_function.setter
    def pressure_function(self, k: _Func1Like) -> None:
        f: Func1
        if isinstance(k, Func1):
            f = k
        else:
            f = Func1(k)
        self._rate_func = f
        self.dev.setPressureFunction(f._func)

    @property
    def time_function(self) -> float:
        r"""
        The time dependence of a flow device. The mass flow rate [kg/s] is calculated
        for a Flow device, and multiplied by a function of time, which returns 1.0
        unless a user-defined function is provided. The calculation of mass flow rate
        depends on the flow device.

        >>> f = FlowDevice(res1, reactor1)
        >>> f.time_function = lambda t: exp(-10 * (t - 0.5)**2)

        where `FlowDevice` is either a `Valve` or `MassFlowController` object.

        .. versionadded:: 3.0
        """
        return self.dev.evalTimeFunction()

    @time_function.setter
    def time_function(self, k: _Func1Like) -> None:
        g: Func1
        if isinstance(k, Func1):
            g = k
        else:
            g = Func1(k)
        self._time_func = g
        self.dev.setTimeFunction(g._func)

    @property
    def device_coefficient(self) -> float:
        """
        Device coefficient (defined by derived class).

        Example:

        >>> v = Valve(res1, reactor1)
        >>> v.device_coefficient = 1e-4  # Set the 'valve coefficient'

        .. versionadded:: 3.2
        """
        return self.dev.deviceCoefficient()

    @device_coefficient.setter
    def device_coefficient(self, value: cython.double) -> None:
        self.dev.setDeviceCoefficient(value)

    def draw(self, graph: "_Digraph | None" = None, *,
              graph_attr: dict[str, str] | None = None,
              node_attr: dict[str, str] | None = None,
              edge_attr: dict[str, str] | None = None) -> "_Digraph":
        """
        Draw as connection between upstream and downstream reactor or reservoir using
        ``graphviz``.

        :param graph:
            ``graphviz.graphs.BaseGraph`` object to which the connection is added.
            If not provided, a new ``DiGraph`` is created. Defaults to ``None``
        :param graph_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the general appearance of the drawn network.
            Has no effect if existing ``graph`` is provided.
            See https://graphviz.org/docs/graph/ for a list of all usable attributes.
        :param node_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the default appearance of any ``node`` (reactors, reservoirs).
            Has no effect if existing ``graph`` is provided.
            See https://graphviz.org/docs/nodes/ for a list of all usable attributes.
        :param edge_attr:
            Attributes to be passed to the ``edge`` method invoked to draw this flow
            controller connection.
            See https://graphviz.org/docs/edges/ for a list of all usable attributes.
        :return:
            A ``graphviz.graphs.BaseGraph`` object depicting the connection.

        .. versionadded:: 3.1
        """
        return draw_flow_controllers([self], graph, graph_attr, node_attr, edge_attr)


@cython.cclass
class MassFlowController(FlowDevice):
    r"""
    A mass flow controller maintains a specified mass
    flow rate independent of upstream and downstream conditions. The equation
    used to compute the mass flow rate is

    .. math:: \dot m = \max(\dot m_0*g(t), 0.),

    where :math:`\dot m_0` is a constant value and :math:`g(t)` is a function of
    time. Both :math:`\dot m_0` and :math:`g(t)` can be set individually by
    properties `mass_flow_coeff` and `time_function`, respectively. The property
    `mass_flow_rate` combines the former into a single interface. Note that if
    :math:`\dot m_0*g(t) < 0`, the mass flow rate will be set to zero, since
    reversal of the flow direction is not allowed.

    Unlike a real mass flow controller, a MassFlowController object will
    maintain the flow even if the downstream pressure is greater than the
    upstream pressure. This allows simple implementation of loops, in which
    exhaust gas from a reactor is fed back into it through an inlet. But note
    that this capability should be used with caution, since no account is
    taken of the work required to do this.
    """
    node_type: _ClassVar[str] = "MassFlowController"

    def __init__(self, upstream: ReactorBase, downstream: ReactorBase, *,
                 name: str = "(none)", mdot: _Func1Like = 1.,
                 edge_attr: dict[str, str] | None = None) -> None:
        super().__init__(upstream, downstream, name=name, edge_attr=edge_attr)
        self.mass_flow_rate = mdot

    @property
    def mass_flow_coeff(self) -> float:
        r"""Set the mass flow rate [kg/s] through the mass flow controller
        as a constant, which may be modified by a function of time, see
        `time_function`.

        >>> mfc = MassFlowController(res1, reactor1)
        >>> mfc.mass_flow_coeff = 1e-4  # Set the flow rate to a constant
        >>> mfc.mass_flow_coeff  # Get the flow rate value
        """
        return cython.cast(cython.pointer(CxxMassFlowController), self.dev).getMassFlowCoeff()

    @mass_flow_coeff.setter
    def mass_flow_coeff(self, value: cython.double) -> None:
        cython.cast(cython.pointer(CxxMassFlowController), self.dev).setMassFlowCoeff(value)

    @property
    def mass_flow_rate(self) -> float:
        r"""
        Set the mass flow rate [kg/s] through this controller to be either
        a constant or an arbitrary function of time. See `Func1`, or get its
        current value.

        Note that depending on the argument type, this method either changes
        the property `mass_flow_coeff` or updates the `time_function` property.

        >>> mfc.mass_flow_rate = 0.3
        >>> mfc.mass_flow_rate = lambda t: 2.5 * exp(-10 * (t - 0.5)**2)
        """
        return self.dev.massFlowRate()

    @mass_flow_rate.setter
    def mass_flow_rate(self, m: _Func1Like) -> None:
        if isinstance(m, _numbers.Real):
            cython.cast(cython.pointer(CxxMassFlowController), self.dev).setMassFlowRate(m)
        else:
            self.mass_flow_coeff = 1.
            self.time_function = m


@cython.cclass
class Valve(FlowDevice):
    r"""
    In Cantera, a `Valve` is a flow device with mass flow rate that is a
    function of the pressure drop across it. The default behavior is linear:

    .. math:: \dot m = K_v*(P_1 - P_2)

    where :math:`K_v` is a constant set using the `valve_coeff` property.
    Note that :math:`P_1` must be greater than :math:`P_2`; otherwise,
    :math:`\dot m = 0`. However, an arbitrary function can also be specified,
    such that

    .. math:: \dot m = K_v*f(P_1 - P_2)

    where :math:`f` is the arbitrary function that multiplies :math:`K_v` given
    a single argument, the pressure differential. Further, a valve opening function
    :math:`g` may be specified using the `time_function` property, such that

    .. math:: \dot m = K_v*g(t)*f(P_1 - P_2)

    See the documentation for the `valve_coeff` property as well as the
    `pressure_function` and `time_function` properties for examples. Note that
    it is never possible for the flow to reverse and go from the downstream to the
    upstream reactor/reservoir through a line containing a `Valve` object.

    `Valve` objects are often used between an upstream reactor and a
    downstream reactor or reservoir to maintain them both at nearly the same
    pressure. By setting the constant :math:`K_v` to a sufficiently large
    value, very small pressure differences will result in flow between the
    reactors that counteracts the pressure difference.
    """
    node_type: _ClassVar[str] = "Valve"

    def __init__(self, upstream: ReactorBase, downstream: ReactorBase, *,
                 name: str = "(none)", K: _Func1Like = 1.,
                 edge_attr: dict[str, str] | None = None) -> None:
        super().__init__(upstream, downstream, name=name, edge_attr=edge_attr)
        if isinstance(K, _numbers.Real):
            self.valve_coeff = K
        else:
            self.valve_coeff = 1.
            self.pressure_function = K

    @property
    def valve_coeff(self) -> float:
        r"""Set valve coefficient, that is, the proportionality constant between mass
        flow rate and pressure drop [kg/s/Pa].

        >>> v = Valve(res1, reactor1)
        >>> v.valve_coeff = 1e-4  # Set the value of K to a constant
        >>> v.valve_coeff  # Get the value of K
        """
        return cython.cast(cython.pointer(CxxValve), self.dev).getValveCoeff()

    @valve_coeff.setter
    def valve_coeff(self, value: cython.double) -> None:
        cython.cast(cython.pointer(CxxValve), self.dev).setValveCoeff(value)


@cython.cclass
class PressureController(FlowDevice):
    r"""
    A PressureController is designed to be used in conjunction with another
    primary flow controller, typically a `MassFlowController`. The primary
    flow controller is installed on the inlet of the reactor, and the
    corresponding `PressureController` is installed on the outlet of the
    reactor. The `PressureController` mass flow rate is equal to the primary
    mass flow rate, plus a small correction dependent on the pressure
    difference:

    .. math:: \dot m = \dot m_{\rm primary} + K_v(P_1 - P_2).

    As an alternative, an arbitrary function of pressure differential can be
    specified using the `pressure_function` property, such that

    .. math:: \dot m = \dot m_{\rm primary} + K_v*f(P_1 - P_2)

    where :math:`f` is the arbitrary function of a single argument.
    """
    node_type: _ClassVar[str] = "PressureController"

    def __init__(self, upstream: ReactorBase, downstream: ReactorBase, *,
                 name: str = "(none)", primary: FlowDevice | None = None,
                 K: _Func1Like = 1., edge_attr: dict[str, str] | None = None) -> None:
        super().__init__(upstream, downstream, name=name, edge_attr=edge_attr)
        if primary is not None:
            self.primary = primary
        if isinstance(K, _numbers.Real):
            self.pressure_coeff = K
        else:
            self.pressure_coeff = 1.
            self.pressure_function = K

    @property
    def pressure_coeff(self) -> float:
        """
        Get/set the proportionality constant :math:`K_v` [kg/s/Pa] between the
        pressure drop and the mass flow rate.
        """
        return cython.cast(cython.pointer(CxxPressureController), self.dev).getPressureCoeff()

    @pressure_coeff.setter
    def pressure_coeff(self, value: cython.double) -> None:
        cython.cast(cython.pointer(CxxPressureController), self.dev).setPressureCoeff(value)

    @property
    def primary(self) -> _Never:
        """
        Primary `FlowDevice` used to compute this device's mass flow rate.

        .. versionadded:: 3.0
        """
        raise NotImplementedError("PressureController.primary")

    @primary.setter
    def primary(self, d: FlowDevice) -> None:
        self.dev.setPrimary(d._node)


@cython.cclass
class ReactorNet:
    """
    Networks of reactors. ReactorNet objects are used to simultaneously
    advance the state of one or more coupled reactors.

    Example:

    >>> r1 = Reactor(gas1)
    >>> r2 = Reactor(gas2)
    >>> <... install walls, inlets, outlets, etc...>

    >>> reactor_network = ReactorNet([r1, r2])
    >>> reactor_network.advance(time)
    """
    def __init__(self, reactors: _Sequence[Reactor] = ()) -> None:
        self._reactors = []  # prevents premature garbage collection
        cxx_reactors: vector[shared_ptr[CxxReactorBase]]
        r: ReactorBase
        for r in reactors:
            self._reactors.append(r)
            cxx_reactors.push_back(r._rbase)
        reactors_span: span[shared_ptr[CxxReactorBase]] = \
            span[shared_ptr[CxxReactorBase]](cxx_reactors)
        self._net = CxxNewReactorNet(reactors_span)
        self.net = self._net.get()

    def advance(self, t: cython.double, apply_limit: pybool = True) -> float:
        """
        Advance the state of the reactor network from the current time/distance towards
        the specified value ``t`` of the independent variable, which depends on the type
        of reactors included in the network.

        The integrator will take as many steps as necessary to reach ``t``. If
        ``apply_limit`` is true and an advance limit is specified, the reactor state at
        the end of the step is estimated prior to advancing. If the difference exceed
        limits, the end value is reduced by half until the projected end state remains
        within specified limits. Returns the time/distance reached at the end of
        integration.
        """
        return self.net.advance(t, apply_limit)

    def step(self) -> float:
        """
        Take a single internal step. The time/distance after taking the step is
        returned.
        """
        return self.net.step()

    def solve_steady(self, loglevel: _LogLevel = 0) -> None:
        """
        Solve directly for the steady-state solution. This approach is generally more
        efficient than time marching to the steady-state (using the
        `advance_to_steady_state` method), but imposes a few limitations:

        - The volume of control volume reactor types (such as `Reactor` and
          `IdealGasMoleReactor`) must be constant; no moving walls can be used.
        - The mass of constant pressure reactor types (such as `ConstPressureReactor`
          and `IdealGasConstPressureReactor`) must be constant; if flow devices are
          used, inlet and outlet flows must be balanced.
        - The solver is currently not compatible with the `ConstPressureMoleReactor` or
          `IdealGasConstPressureMoleReactor` classes.
        - Only ideal gas reactor types can be used for when the energy equation is
          disabled (fixed temperature simulations).
        - Reacting surfaces are not yet supported.

        :param loglevel:
            Print information about solver progress to aid in understanding
            cases where the solver fails to converge. Higher levels are more verbose.

            - 0: No logging.
            - 1: Basic info about each steady-state attempt and round of time stepping.
            - 2: Adds details about each time step and steady-state Newton iteration.
            - 3: Adds details about Newton iterations for each time step.
            - 4: Adds details about state variables that are limiting steady-state
              Newton step sizes.
            - 5: Adds details about state variables that are limiting time-stepping
              Newton step sizes.
            - 6: Print current state vector after different solver stages
            - 7: Print current residual vector after different solver stages

        Uses the hybrid time marching / damped Newton's method solver implemented by
        classes :ct:`SteadyStateSystem` and :ct:`MultiNewton`.

        .. versionadded:: 3.2
        """
        self.net.solveSteady(loglevel)

    def steady_jacobian(self, rdt: cython.float = 0.0) -> _Array:
        """
        Get the Jacobian used by the steady-state solver.

        :param rdt:
            Reciprocal of the pseudo-timestep [1/s]. Default of 0.0 returns the
            steady-state Jacobian.

        .. versionadded:: 3.2
        """
        return get_from_sparse(self.net.steadyJacobian(rdt), self.n_vars, self.n_vars)

    @property
    def jacobian(self) -> _Array:
        """
        Get the analytical preconditioner Jacobian for the reactor network as a
        sparse matrix.

        Collects entries from each reactor's ``getJacobianElements()`` implementation
        and assembles them into a single network-level matrix using global row and
        column indices. Reactors that do not implement ``getJacobianElements()``
        contribute no entries.

        This property is useful for debugging custom Jacobian implementations in
        :py:class:`ExtensibleReactor` subclasses — for example, to verify that
        elements are placed in the correct global positions and have the expected
        values.

        .. warning::

            Depending on the particular implementation, this may return an approximate
            Jacobian intended only for use in forming a preconditioner for iterative
            solvers, excluding terms that would generate a fully-dense Jacobian.

        .. warning::

            This property is an experimental part of the Cantera API and may be
            changed or removed without notice.

        .. versionadded:: 4.0
        """
        return get_from_sparse(self.net.jacobian(), self.n_vars, self.n_vars)

    @property
    def finite_difference_jacobian(self) -> _Array:
        """
        Get the system Jacobian for the reactor network, calculated using finite
        differences, as a sparse matrix.

        Perturbs each element of the network state vector and evaluates the network
        RHS using central differences to estimate the full Jacobian. Perturbation
        step sizes are scaled by the integrator tolerances. This method is intended
        for debugging and validation of analytical Jacobian implementations.

        .. warning::

            This property is an experimental part of the Cantera API and may be
            changed or removed without notice.

        .. versionadded:: 4.0
        """
        return get_from_sparse(self.net.finiteDifferenceJacobian(),
                               self.n_vars, self.n_vars)

    def initialize(self) -> None:
        """
        Force initialization of the integrator after initial setup.
        """
        self.net.initialize()

    def reinitialize(self) -> None:
        """
        Reinitialize the integrator after making changing to the state of the
        system. Changes to Reactor contents will automatically trigger
        reinitialization.
        """
        self.net.reinitialize()

    @property
    def reactors(self) -> list[Reactor]:
        """
        List of all reactors that are part of the reactor network.

        .. versionadded:: 3.1
        """
        return self._reactors

    @property
    def time(self) -> float:
        """
        The current time [s], for reactor networks that are solved in the time domain.
        """
        return self.net.time()

    @property
    def distance(self) -> float:
        """
        The current distance [m] along the length of the reactor network, for reactors
        that are solved as a function of space.
        """
        return self.net.distance()

    @property
    def initial_time(self) -> float:
        """
        The initial time of the integrator. When set, integration is restarted from this
        time using the current state as the initial condition. Default: 0.0 s.

        .. versionadded:: 3.0
        """
        return self.net.getInitialTime()

    @initial_time.setter
    def initial_time(self, t: cython.double) -> None:
        self.net.setInitialTime(t)

    @property
    def max_time_step(self) -> float:
        """
        Get/set the maximum time step *t* [s] that the integrator is
        allowed to use. The default value is set to zero, so that no time
        step maximum is used.
        """
        return self.net.maxTimeStep()

    @max_time_step.setter
    def max_time_step(self, t: cython.double) -> None:
        self.net.setMaxTimeStep(t)

    @property
    def max_err_test_fails(self) -> _Never:
        """
        The maximum number of error test failures permitted by the CVODES integrator
        in a single step. The default is 10.
        """
        raise AttributeError("unreadable attribute 'max_err_test_fails'")

    @max_err_test_fails.setter
    def max_err_test_fails(self, n: int) -> None:
        self.net.setMaxErrTestFails(n)

    @property
    def max_nonlinear_iterations(self) -> int:
        """
        Get/Set the maximum number of nonlinear solver iterations permitted by the
        SUNDIALS solver in one solve attempt. The default value is 4.
        """
        return self.net.integrator().maxNonlinIterations()

    @max_nonlinear_iterations.setter
    def max_nonlinear_iterations(self, n: cython.int) -> None:
        self.net.integrator().setMaxNonlinIterations(n)

    @property
    def max_nonlinear_convergence_failures(self) -> int:
        """
        Get/Set the maximum number of nonlinear solver convergence failures permitted in
        one step of the SUNDIALS integrator. The default value is 10.
        """
        return self.net.integrator().maxNonlinConvFailures()

    @max_nonlinear_convergence_failures.setter
    def max_nonlinear_convergence_failures(self, n: cython.int) -> None:
        self.net.integrator().setMaxNonlinConvFailures(n)

    @property
    def include_algebraic_in_error_test(self) -> bool:
        """
        Get/Set whether to include algebraic variables in the in the local error test.
        Applicable only to DAE systems. The default is `True`.
        """
        return self.net.integrator().algebraicInErrorTest()

    @include_algebraic_in_error_test.setter
    def include_algebraic_in_error_test(self, yesno: pybool) -> None:
        self.net.integrator().includeAlgebraicInErrorTest(yesno)

    @property
    def max_order(self) -> _Literal[1, 2, 3, 4, 5]:
        """
        Get/Set the maximum order of the linear multistep method. The default value and
        maximum is 5.
        """
        return self.net.integrator().maxOrder()

    @max_order.setter
    def max_order(self, n: _Literal[1, 2, 3, 4, 5]) -> None:
        self.net.integrator().setMaxOrder(n)

    @property
    def max_steps(self) -> int:
        """
        The maximum number of internal integration steps that CVODES
        is allowed to take before reaching the next output point.
        """
        return self.net.maxSteps()

    @max_steps.setter
    def max_steps(self, nsteps: int) -> None:
        self.net.setMaxSteps(nsteps)

    @property
    def rtol(self) -> float:
        """
        The relative error tolerance used while integrating the reactor
        equations.
        """
        return self.net.rtol()

    @rtol.setter
    def rtol(self, tol: float) -> None:
        self.net.setRelativeTolerance(tol)

    @property
    def atol(self) -> float:
        """
        The scalar absolute error tolerance used while integrating the reactor
        equations. This value is used for state variables that do not have
        reactor-specific absolute tolerances set by `Reactor.atol`. Setting
        this property to ``None`` clears the scalar override and restores
        reactor-specific default absolute tolerance scaling.
        """
        return self.net.atol()

    @atol.setter
    def atol(self, tol: float | None) -> None:
        if tol is None:
            self.net.clearAbsoluteTolerance()
        else:
            self.net.setAbsoluteTolerance(tol)

    @property
    def rtol_sensitivity(self) -> float:
        """
        The relative error tolerance for sensitivity analysis.
        """
        return self.net.rtolSensitivity()

    @rtol_sensitivity.setter
    def rtol_sensitivity(self, tol: float) -> None:
        self.net.setSensitivityTolerances(tol, -1)

    @property
    def atol_sensitivity(self) -> float:
        """
        The absolute error tolerance for sensitivity analysis.
        """
        return self.net.atolSensitivity()

    @atol_sensitivity.setter
    def atol_sensitivity(self, tol: float) -> None:
        self.net.setSensitivityTolerances(-1, tol)

    @property
    def verbose(self) -> bool:
        """
        If `True`, verbose debug information will be printed during
        integration. The default is `False`.
        """
        return pybool(self.net.verbose())

    @verbose.setter
    def verbose(self, v: pybool) -> None:
        self.net.setVerbose(v)

    def global_component_index(self, name: str, reactor: cython.int) -> int:
        """
        Returns the index of a component named ``name`` of a reactor with index
        ``reactor`` within the global state vector. That is, this determines the
        absolute index of the component, where ``reactor`` is the index of the
        reactor that holds the component. ``name`` is either a species name or the
        name of a reactor state variable, for example, ``'int_energy'``, ``'temperature'``, etc.
        depending on the reactor's equations.
        """
        return self.net.globalComponentIndex(stringify(name), reactor)

    def component_name(self, i: cython.int) -> str:
        """
        Return the name of the i-th component of the global state vector. The
        name returned includes both the name of the reactor and the specific
        component, for example `'reactor1: CH4'`.
        """
        return pystr(self.net.componentName(i))

    def sensitivity(self, component: int | str, p: cython.int,
                     r: cython.int = 0) -> float:
        """
        Returns the sensitivity of the solution variable ``component`` in
        reactor ``r`` with respect to the parameter ``p``. ``component`` can be a
        string or an integer. See `component_index` and `sensitivities` to
        determine the integer index for the variables and the definition of the
        resulting sensitivity coefficient. If it is not given, ``r`` defaults to
        the first reactor. Returns an empty array until the first integration step is
        taken.
        """
        if isinstance(component, int):
            return self.net.sensitivity(component, p)
        elif isinstance(component, (str, bytes)):
            return self.net.sensitivity(stringify(component), p, r)

    def sensitivities(self) -> _Array:
        r"""
        Returns the sensitivities of all of the solution variables with respect
        to all of the registered parameters. The normalized sensitivity
        coefficient :math:`S_{ki}` of the solution variable :math:`y_k` with
        respect to sensitivity parameter :math:`p_i` is defined as:

        .. math:: S_{ki} = \frac{1}{y_k} \frac{\partial y_k}{\partial p_i}

        For reaction sensitivities, the parameter is a multiplier on the forward
        rate constant (and implicitly on the reverse rate constant for
        reversible reactions) which has a nominal value of 1.0, and the sensitivity
        is nondimensional.

        For species enthalpy sensitivities, the parameter is an additive perturbation to
        the molar enthalpy of formation, such that the dimensions of the sensitivity
        are kmol/J.

        The sensitivities are returned in an array with dimensions *(n_vars,
        n_sensitivity_params)*, unless no integration steps have been taken, in which
        case the shape is *(0, n_sensitivity_params)*. The order of the
        variables (that is, rows) is:

        `Reactor` or `IdealGasReactor`:

        - 0  - mass
        - 1  - volume
        - 2  - internal energy or temperature
        - 3+ - mass fractions of the species

        `ConstPressureReactor` or `IdealGasConstPressureReactor`:

        - 0  - mass
        - 1  - enthalpy or temperature
        - 2+ - mass fractions of the species
        """
        data = np.empty((self.n_vars, self.n_sensitivity_params))
        p: cython.int
        k: cython.int
        for p in range(self.n_sensitivity_params):
            for k in range(self.n_vars):
                data[k, p] = self.net.sensitivity(k, p)
        return data

    def sensitivity_parameter_name(self, p: cython.int) -> str:
        """
        Name of the sensitivity parameter with index ``p``.
        """
        return pystr(self.net.sensitivityParameterName(p))

    @property
    def n_sensitivity_params(self) -> int:
        """
        The number of registered sensitivity parameters.
        """
        return self.net.nparams()

    @property
    def n_vars(self) -> int:
        """
        The number of state variables in the system. This is the sum of the
        number of variables for each `Reactor` and `Wall` in the system.
        Equal to:

        `Reactor` and `IdealGasReactor`: `n_species` + 3 (mass, volume,
        internal energy or temperature).

        `ConstPressureReactor` and `IdealGasConstPressureReactor`:
        `n_species` + 2 (mass, enthalpy or temperature).

        `Wall`: number of surface species
        """
        return self.net.neq()

    def get_state(self) -> _Array:
        """
        Get the combined state vector of the reactor network.

        The combined state vector consists of the concatenated state vectors of
        all entities contained.
        """
        if not self.n_vars:
            raise CanteraError('ReactorNet empty or not initialized.')
        y = np.zeros(self.n_vars)
        cy: cython.double[::1] = y
        self.net.getState(span[double](cython.address(cy[0]),
                                      cython.cast(cython.size_t, y.size)))
        return y

    def get_state_dae(self) -> tuple[_Array, _Array]:
        """
        Get the combined state vector and its time derivative for reactor networks
        formulated as DAEs (namely, those containing `FlowReactor`).

        The combined state vector consists of the concatenated state vectors of
        all entities contained.

        .. versionadded:: 4.0
        """
        y = np.zeros(self.n_vars)
        yp = np.zeros(self.n_vars)
        cy: cython.double[::1] = y
        cyp: cython.double[::1] = yp
        self.net.getStateDae(
            span[double](cython.address(cy[0]), cython.cast(cython.size_t, y.size)),
            span[double](cython.address(cyp[0]), cython.cast(cython.size_t, yp.size)))
        return y, yp

    def get_derivative(self, k: int) -> _Array:
        """
        Get the k-th derivative of the state vector of the reactor network with respect
        to the independent integrator variable (time/distance).
        """
        if not self.n_vars:
            raise CanteraError('ReactorNet empty or not initialized.')
        dky = np.zeros(self.n_vars)
        cdky: cython.double[::1] = dky
        self.net.getDerivative(k, span[double](cython.address(cdky[0]),
                                               cython.cast(cython.size_t, dky.size)))
        return dky

    @property
    def advance_limits(self) -> _Array:
        """
        Get or set absolute limits for state changes during `ReactorNet.advance`
        (positive values are considered; negative values disable a previously
        set advance limit for a solution component). Note that limits are
        disabled by default (with individual values set to -1.).
        """
        limits = np.empty(self.n_vars)
        climits: cython.double[::1] = limits
        self.net.getAdvanceLimits(span[double](cython.address(climits[0]),
                                              cython.cast(cython.size_t, limits.size)))
        return limits

    @advance_limits.setter
    def advance_limits(self, limits: _ArrayLike | None) -> None:
        if limits is None:
            limits = -1. * np.ones([self.n_vars])
        elif len(limits) != self.n_vars:
            raise ValueError('array must be of length n_vars')

        data = np.ascontiguousarray(limits, dtype=np.double)
        cdata: cython.double[::1] = data
        self.net.setAdvanceLimits(span[double](cython.address(cdata[0]),
                                              cython.cast(cython.size_t, data.size)))

    @_overload
    def advance_to_steady_state(self, max_steps: int, residual_threshold: float,
                                atol: float,
                                return_residuals: _Literal[False] = False) -> None: ...
    @_overload
    def advance_to_steady_state(self, max_steps: int, residual_threshold: float,
                                atol: float,
                                return_residuals: _Literal[True]) -> _Array: ...
    @_overload
    def advance_to_steady_state(self, max_steps: int = 10000,
                                residual_threshold: float = 0.0, atol: float = 0.0,
                                return_residuals: bool = False) -> _Array | None: ...
    def advance_to_steady_state(self, max_steps: cython.int = 10000,
                                residual_threshold: cython.double = 0.,
                                atol: cython.double = 0.,
                                return_residuals: pybool = False):
        r"""
        Advance the reactor network in time until steady state is reached.

        The steady state is defined by requiring that the state of the system
        only changes below a certain threshold. The residual is computed using
        feature scaling:

        .. math:: r = \left| \frac{x(t + \Delta t) - x(t)}{\text{max}(x) + \text{atol}} \right| \cdot \frac{1}{\sqrt{n_x}}

        :param max_steps:
            Maximum number of steps to be taken
        :param residual_threshold:
            Threshold below which the feature-scaled residual r should drop such
            that the network is defines as steady state. By default,
            residual_threshold is 10 times the solver rtol.
        :param atol:
            The smallest expected value of interest. Used for feature scaling.
            By default, this atol is identical to the solver atol.
        :param return_residuals:
            If set to `True`, this function returns the residual time series
            as a vector with length `max_steps`.

        """
        # get default tolerances:
        if not atol:
            atol = self.atol
        if not residual_threshold:
            residual_threshold = 10. * self.rtol
        if residual_threshold <= self.rtol:
            raise CanteraError('Residual threshold (' + str(residual_threshold) +
                               ') should be below solver rtol (' +
                               str(self.rtol) + ')')
        if return_residuals:
            residuals = np.empty(max_steps)
        # check if system is initialized
        if not self.n_vars:
            self.reinitialize()
        max_state_values = self.get_state()  # denominator for feature scaling
        for step in range(max_steps):
            previous_state = self.get_state()
            # take 10 steps (just to increase speed)
            for n1 in range(10):
                self.step()
            state = self.get_state()
            max_state_values = np.maximum(max_state_values, state)
            # determine feature_scaled residual
            residual = np.linalg.norm((state - previous_state)
                / (max_state_values + atol)) / np.sqrt(self.n_vars)
            if return_residuals:
                residuals[step] = residual
            if residual < residual_threshold:
                break
        if step == max_steps - 1:
            raise CanteraError('Maximum number of steps reached before'
                               ' convergence below maximum residual')
        if return_residuals:
            return residuals[:step + 1]

    def __reduce__(self) -> _Never:
        raise NotImplementedError('ReactorNet object is not picklable')

    def __copy__(self) -> _Never:
        raise NotImplementedError('ReactorNet object is not copyable')

    @property
    def preconditioner(self) -> SystemJacobian:
        """Preconditioner associated with integrator"""
        raise AttributeError("unreadable attribute 'preconditioner'")

    @preconditioner.setter
    def preconditioner(self, precon: SystemJacobian) -> None:
        # set preconditioner
        self.net.setPreconditioner(precon._base)
        # set problem type as default of preconditioner
        self.linear_solver_type = precon.linear_solver_type

    @property
    def linear_solver_type(self) -> _Literal["DENSE", "GMRES", "BAND", "DIAG"]:
        """
            The type of linear solver used in integration.

            Options for this property include:

            - `"DENSE"`
            - `"GMRES"`
            - `"BAND"`
            - `"DIAG"`

        """
        return pystr(self.net.linearSolverType())

    @linear_solver_type.setter
    def linear_solver_type(
            self, linear_solver_type: _Literal["DENSE", "GMRES", "BAND", "DIAG"]
    ) -> None:
        self.net.setLinearSolverType(stringify(linear_solver_type))

    @property
    def solver_stats(self) -> _SolverStats:
        """ODE solver stats from integrator"""
        stats: CxxAnyMap
        stats = self.net.solverStats()
        return anymap_to_py(stats)

    @property
    def derivative_settings(self) -> _Never:
        """
        Apply derivative settings to all reactors in the network.
        See also `Kinetics.derivative_settings`.
        """
        raise AttributeError("unreadable attribute 'derivative_settings'")

    @derivative_settings.setter
    def derivative_settings(self, settings: "_DerivativeSettings") -> None:
        self.net.setDerivativeSettings(py_to_anymap(settings))

    def draw(self, *, graph_attr: dict[str, str] | None = None,
             node_attr: dict[str, str] | None = None,
             edge_attr: dict[str, str] | None = None,
             heat_flow_attr: dict[str, str] | None = None,
             mass_flow_attr: dict[str, str] | None = None,
             moving_wall_edge_attr: dict[str, str] | None = None,
             surface_edge_attr: dict[str, str] | None = None,
             show_wall_velocity: bool = True, print_state: bool = False,
             species: _Literal["X", "Y"] | bool | _Iterable[str] | None = None,
             species_units: _Literal["percent", "ppm"] = "percent") -> "_Digraph":
        """
        Draw as ``graphviz.graphs.DiGraph``. Connecting flow controllers and
        walls are depicted as arrows.

        :param graph_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the general appearance of the drawn network.
            See https://graphviz.org/docs/graph/ for a list of all usable attributes.
        :param node_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the default appearance of any ``node`` (reactors, reservoirs).
            ``node_attr`` defined in the reactor object itself have priority.
            See https://graphviz.org/docs/nodes/ for a list of all usable attributes.
        :param edge_attr:
            Attributes to be passed to the ``graphviz.Digraph`` function that control
            the default appearance of any ``edge`` (flow controllers, walls).
            ``edge_attr`` defined in the connection objects (subclasses of `FlowDevice`
            or walls) themselves have priority.
            See https://graphviz.org/docs/edges/ for a list of all usable attributes.
        :param heat_flow_attr:
            Same as ``edge_attr`` but only applied to edges representing walls.
            Default is ``{"color": "red", "style": "dashed"}``.
        :param mass_flow_attr:
            Same as ``edge_attr`` but only applied to edges representing `FlowDevice`
            objects.
        :param moving_wall_edge_attr:
            Same as ``edge_attr`` but only applied to edges representing wall movement.
        :param surface_edge_attr:
            Same as ``edge_attr`` but only applied to edges representing connections
            between `ReactorSurface` objects and reactors.
            Default is ``{"style": "dotted", "arrowhead": "none"}``.
        :param show_wall_velocity:
            If ``True``, wall movement will be indicated by additional arrows with the
            corresponding wall velocity as a label.
        :param print_state:
            Whether state information of the reactors is printed into each node.
            Defaults to ``False``.
        :param species:
            If ``print_state`` is ``True``, define how species are to be printed.
            Options are ``'X'`` and ``'Y'`` for mole and mass fractions of all species,
            respectively, or an iterable that contains the desired species names as
            strings. Defaults to ``None``.
        :param species_units:
            Defines the units the species are displayed in as either ``"percent"`` or
            ``"ppm"``. Defaults to ``"percent"``.
        :return:
            ``graphviz.graphs.BaseGraph`` object with reactor net.

        .. versionadded:: 3.1
        """
        return draw_reactor_net(self, graph_attr, node_attr, edge_attr,
                                heat_flow_attr, mass_flow_attr, moving_wall_edge_attr,
                                surface_edge_attr, show_wall_velocity, print_state,
                                species, species_units)
