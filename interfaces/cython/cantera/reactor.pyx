# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import warnings
from collections import defaultdict as _defaultdict
import numbers as _numbers
from cython.operator cimport dereference as deref

from .thermo cimport *
from ._utils cimport pystr, stringify, comp_map, py_to_anymap, anymap_to_py
from ._utils import *
from .delegator cimport *
from .drawnetwork import *

_reactor_counts = _defaultdict(int)

cdef class ReactorBase:
    """
    Common base class for reactors and reservoirs.
    """
    reactor_type = "none"
    def __cinit__(self, _SolutionBase contents=None, name=None, *, **kwargs):
        def reactor_name(name):
            if name is not None:
                return name
            _reactor_counts[self.reactor_type] += 1
            return f"{self.reactor_type}_{_reactor_counts[self.reactor_type]}"

        if isinstance(contents, _SolutionBase):
            self._reactor = newReactor(stringify(self.reactor_type),
                                       contents._base, stringify(reactor_name(name)))
        else:
            # deprecated: will raise warnings in C++ layer
            self._reactor = newReactor(stringify(self.reactor_type))
            self._reactor.get().setName(stringify(reactor_name(name)))
        self.rbase = self._reactor.get()

    def __init__(self, _SolutionBase contents=None, name=None, *, volume=None,
                 node_attr=None):
        self._inlets = []
        self._outlets = []
        self._walls = []
        self._surfaces = []
        self._edges = []
        if isinstance(contents, _SolutionBase):
            self.insert(contents)  # leave insert for the time being

        if volume is not None:
            self.volume = volume

        self.node_attr = node_attr or {}

    def insert(self, _SolutionBase solution):
        """
        Set ``solution`` to be the object used to compute thermodynamic
        properties and kinetic rates for this reactor.
        """
        self._thermo = solution
        self.rbase.setSolution(solution._base)

    property type:
        """The type of the reactor."""
        def __get__(self):
            return pystr(self.rbase.type())

    property name:
        """The name of the reactor."""
        def __get__(self):
            return pystr(self.rbase.name())

        def __set__(self, name):
            self.rbase.setName(stringify(name))

    def syncState(self):
        """
        Set the state of the Reactor to match that of the associated
        `ThermoPhase` object. After calling syncState(), call
        ReactorNet.reinitialize() before further integration.
        """
        self.rbase.syncState()

    property thermo:
        """The `ThermoPhase` object representing the reactor's contents."""
        def __get__(self):
            self.rbase.restoreState()
            return self._thermo

    property volume:
        """The volume [m^3] of the reactor."""
        def __get__(self):
            return self.rbase.volume()

        def __set__(self, double value):
            self.rbase.setInitialVolume(value)

    property T:
        """The temperature [K] of the reactor's contents."""
        def __get__(self):
            return self.thermo.T

    property density:
        """The density [kg/m^3 or kmol/m^3] of the reactor's contents."""
        def __get__(self):
            return self.thermo.density

    property mass:
        """The mass of the reactor's contents."""
        def __get__(self):
            return self.thermo.density_mass * self.volume

    property Y:
        """The mass fractions of the reactor's contents."""
        def __get__(self):
            return self.thermo.Y

    # Flow devices & walls
    property inlets:
        """List of flow devices installed as inlets to this reactor"""
        def __get__(self):
            return self._inlets

    property outlets:
        """List of flow devices installed as outlets to this reactor"""
        def __get__(self):
            return self._outlets

    property walls:
        """List of walls installed on this reactor"""
        def __get__(self):
            return self._walls

    property surfaces:
        """List of reacting surfaces installed on this reactor"""
        def __get__(self):
            return self._surfaces

    property edges:
        """List of reacting edges installed on this reactor"""
        def __get__(self):
            return self._edges

    def _add_inlet(self, inlet):
        """
        Store a reference to ``inlet`` to prevent it from being prematurely
        garbage collected.
        """
        self._inlets.append(inlet)

    def _add_outlet(self, outlet):
        """
        Store a reference to ``outlet`` to prevent it from being prematurely
        garbage collected.
        """
        self._outlets.append(outlet)

    def _add_wall(self, wall):
        """
        Store a reference to ``wall`` to prevent it from being prematurely
        garbage collected.
        """
        self._walls.append(wall)

    def draw(self, graph=None, *, graph_attr=None, node_attr=None, print_state=False,
             species=None, species_units="percent"):
        """
        Draw as ``graphviz`` ``dot`` node.
        The node is added to an existing ``graph`` if provided.
        Optionally include current reactor state in the node.

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

    def __reduce__(self):
        raise NotImplementedError('Reactor object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Reactor object is not copyable')


cdef class Reactor(ReactorBase):
    """
    A homogeneous zero-dimensional reactor. By default, they are closed
    (no inlets or outlets), have fixed volume, and have adiabatic,
    chemically-inert walls. These properties may all be changed by adding
    appropriate components such as `Wall`, `MassFlowController` and `Valve`.
    """
    reactor_type = "Reactor"

    def __cinit__(self, *args, **kwargs):
        self.reactor = <CxxReactor*>(self.rbase)

    def __init__(self, contents=None, *, name=None, energy='on', group_name="", **kwargs):
        """
        :param contents:
            Reactor contents. If not specified, the reactor is initially empty.
            In this case, call `insert` to specify the contents. Providing valid
            contents will become mandatory after Cantera 3.1.
        :param name:
            Used only to identify this reactor in output. If not specified,
            defaults to ``'Reactor_n'``, where *n* is an integer assigned in
            the order `Reactor` objects are created.
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

        Some examples showing how to create :class:`Reactor` objects are
        shown below.

        >>> gas = Solution('gri30.yaml')
        >>> r1 = Reactor(gas)

        Arguments may be specified using keywords in any order:

        >>> r2 = Reactor(contents=gas, energy='off',
        ...              name='isothermal_reactor')
        >>> r3 = Reactor(name='adiabatic_reactor', contents=gas)

        """
        super().__init__(contents, name, **kwargs)

        if energy == 'off':
            self.energy_enabled = False
        elif energy != 'on':
            raise ValueError("'energy' must be either 'on' or 'off'")

        self.group_name = group_name

    def insert(self, _SolutionBase solution):
        """
        Set ``solution`` to be the object used to compute thermodynamic
        properties and kinetic rates for this reactor.
        """
        ReactorBase.insert(self, solution)
        self._kinetics = solution

    property kinetics:
        """
        The `Kinetics` object used for calculating kinetic rates in
        this reactor.
        """
        def __get__(self):
            self.rbase.restoreState()
            return self._kinetics

    property chemistry_enabled:
        """
        `True` when the reactor composition is allowed to change due to
        chemical reactions in this reactor. When this is `False`, the
        reactor composition is held constant.
        """
        def __get__(self):
            return self.reactor.chemistryEnabled()

        def __set__(self, pybool value):
            self.reactor.setChemistry(value)

    property energy_enabled:
        """
        `True` when the energy equation is being solved for this reactor.
        When this is `False`, the reactor temperature is held constant.
        """
        def __get__(self):
            return self.reactor.energyEnabled()

        def __set__(self, pybool value):
            self.reactor.setEnergy(int(value))

    def add_sensitivity_reaction(self, m):
        """
        Specifies that the sensitivity of the state variables with respect to
        reaction ``m`` should be computed. ``m`` is the 0-based reaction index.
        The reactor must be part of a network first. Specifying the same
        reaction more than one time raises an exception.
        """
        self.reactor.addSensitivityReaction(m)

    def add_sensitivity_species_enthalpy(self, k):
        """
        Specifies that the sensitivity of the state variables with respect to
        species ``k`` should be computed. The reactor must be part of a network
        first.
        """
        self.reactor.addSensitivitySpeciesEnthalpy(self.thermo.species_index(k))

    def component_index(self, name):
        """
        Returns the index of the component named ``name`` in the system. This determines
        the index of the component in the vector of sensitivity coefficients. ``name``
        is either a species name or the name of a reactor state variable, for example
        ``'int_energy'`` or ``'temperature'``, depending on the reactor's equations.
        """
        k = self.reactor.componentIndex(stringify(name))
        if k == -1:
            raise IndexError('No such component: {!r}'.format(name))
        return k

    def component_name(self, int i):
        """
        Returns the name of the component with index ``i`` within the array of
        variables returned by `get_state`. This is the inverse of
        `component_index`.
        """
        return pystr(self.reactor.componentName(i))

    property n_vars:
        """
        The number of state variables in the reactor.
        Equal to:

        `Reactor` and `IdealGasReactor`: `n_species` + 3 (mass, volume,
        internal energy or temperature).

        `ConstPressureReactor` and `IdealGasConstPressureReactor`:
        `n_species` + 2 (mass, enthalpy or temperature).
        """
        def __get__(self):
            return self.reactor.neq()

    def get_state(self):
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
        cdef np.ndarray[np.double_t, ndim=1] y = np.zeros(self.n_vars)
        self.reactor.getState(&y[0])
        return y

    property jacobian:
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
        def __get__(self):
            return get_from_sparse(self.reactor.jacobian(), self.n_vars, self.n_vars)

    property finite_difference_jacobian:
        """
        Get the reactor-specific Jacobian, calculated using a finite difference method.

        .. warning::

            This property is an experimental part of the Cantera API and
            may be changed or removed without notice.
        """
        def __get__(self):
            return get_from_sparse(self.reactor.finiteDifferenceJacobian(),
                                   self.n_vars, self.n_vars)

    def set_advance_limit(self, name, limit):
        """
        Limit absolute change of component ``name`` during `ReactorNet.advance`.
        (positive ``limit`` values are considered; negative values disable a
        previously set advance limit for a solution component). Note that
        limits are disabled by default (with individual values set to -1.).
        """
        if limit is None:
            limit = -1.
        self.reactor.setAdvanceLimit(stringify(name), limit)

cdef class MoleReactor(Reactor):
    """
    A homogeneous zero-dimensional reactor with a mole based state vector. By default,
    they are closed (no inlets or outlets), have fixed volume, and have adiabatic,
    chemically-inert walls. These properties may all be changed by adding
    appropriate components such as `Wall`, `MassFlowController` and `Valve`.

    .. versionadded:: 3.0
    """
    reactor_type = "MoleReactor"

cdef class Reservoir(ReactorBase):
    """
    A reservoir is a reactor with a constant state. The temperature,
    pressure, and chemical composition in a reservoir never change from
    their initial values.
    """
    reactor_type = "Reservoir"


cdef class ConstPressureReactor(Reactor):
    """A homogeneous, constant pressure, zero-dimensional reactor. The volume
    of the reactor changes as a function of time in order to keep the
    pressure constant.
    """
    reactor_type = "ConstPressureReactor"

cdef class ConstPressureMoleReactor(Reactor):
    """A homogeneous, constant pressure, zero-dimensional reactor with a mole based
    state vector. The volume of the reactor changes as a function of time in order to
    keep the pressure constant.

    .. versionadded:: 3.0
    """
    reactor_type = "ConstPressureMoleReactor"


cdef class IdealGasReactor(Reactor):
    """ A constant volume, zero-dimensional reactor for ideal gas mixtures. """
    reactor_type = "IdealGasReactor"


cdef class IdealGasMoleReactor(Reactor):
    """
    A constant volume, zero-dimensional reactor for ideal gas mixtures with a mole
    based state vector
    """
    reactor_type = "IdealGasMoleReactor"


cdef class IdealGasConstPressureReactor(Reactor):
    """
    A homogeneous, constant pressure, zero-dimensional reactor for ideal gas
    mixtures. The volume of the reactor changes as a function of time in order
    to keep the pressure constant.
    """
    reactor_type = "IdealGasConstPressureReactor"

cdef class IdealGasConstPressureMoleReactor(Reactor):
    """
    A homogeneous, constant pressure, zero-dimensional reactor for ideal gas
    mixtures. The volume of the reactor changes as a function of time in order
    to keep the pressure constant. This reactor also uses a mole based state vector.
    """
    reactor_type = "IdealGasConstPressureMoleReactor"


cdef class FlowReactor(Reactor):
    """
    A steady-state plug flow reactor with constant cross sectional area.
    Integration follows a fluid element along the length of the reactor.
    The reactor is assumed to be frictionless and adiabatic.
    """
    reactor_type = "FlowReactor"

    property mass_flow_rate:
        """ Mass flow rate [kg/s] """
        def __set__(self, double value):
            (<CxxFlowReactor*>self.reactor).setMassFlowRate(value)

    @property
    def area(self):
        """
        Get/set the area of the reactor [m^2].

        When the area is changed, the flow speed is scaled to keep the total mass flow
        rate constant.
        """
        return (<CxxFlowReactor*>self.reactor).area()

    @area.setter
    def area(self, area):
        (<CxxFlowReactor*>self.reactor).setArea(area)

    @property
    def inlet_surface_atol(self):
        """
        Get/Set the steady-state tolerances used to determine the initial surface
        species coverages.
        """
        return (<CxxFlowReactor*>self.reactor).inletSurfaceAtol()

    @inlet_surface_atol.setter
    def inlet_surface_atol(self, atol):
        (<CxxFlowReactor*>self.reactor).setInletSurfaceAtol(atol)

    @property
    def inlet_surface_rtol(self):
        """
        Get/Set the steady-state tolerances used to determine the initial surface
        species coverages.
        """
        return (<CxxFlowReactor*>self.reactor).inletSurfaceRtol()

    @inlet_surface_rtol.setter
    def inlet_surface_rtol(self, rtol):
        (<CxxFlowReactor*>self.reactor).setInletSurfaceRtol(rtol)

    @property
    def inlet_surface_max_steps(self):
        """
        Get/Set the maximum number of integrator steps used to determine the initial
        surface species coverages.
        """
        return (<CxxFlowReactor*>self.reactor).inletSurfaceMaxSteps()

    @inlet_surface_max_steps.setter
    def inlet_surface_max_steps(self, nsteps):
        (<CxxFlowReactor*>self.reactor).setInletSurfaceMaxSteps(nsteps)

    @property
    def inlet_surface_max_error_failures(self):
        """
        Get/Set the maximum number of integrator error failures allowed when determining
        the initial surface species coverages.
        """
        return (<CxxFlowReactor*>self.reactor).inletSurfaceMaxErrorFailures()

    @inlet_surface_max_error_failures.setter
    def inlet_surface_max_error_failures(self, nsteps):
        (<CxxFlowReactor*>self.reactor).setInletSurfaceMaxErrorFailures(nsteps)

    @property
    def surface_area_to_volume_ratio(self):
        """ Get/Set the surface area to volume ratio of the reactor [m^-1] """
        return (<CxxFlowReactor*>self.reactor).surfaceAreaToVolumeRatio()

    @surface_area_to_volume_ratio.setter
    def surface_area_to_volume_ratio(self, sa_to_vol):
        (<CxxFlowReactor*>self.reactor).setSurfaceAreaToVolumeRatio(sa_to_vol)

    @property
    def speed(self):
        """ Speed [m/s] of the flow in the reactor at the current position """
        return (<CxxFlowReactor*>self.reactor).speed()


cdef class ExtensibleReactor(Reactor):
    """
    A base class for a reactor with delegated methods where the base
    functionality corresponds to the `Reactor` class.

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
        Responsible for allocating and setting the sizes of any internal
        variables, initializing attached walls, and setting the total number of
        state variables associated with this reactor, `n_vars`.

        Called once before the start of time integration.

    ``sync_state(self) -> None``
        Responsible for setting the state of the reactor to correspond to the
        state of the associated ThermoPhase object.

    ``get_state(self, y : double[:]) -> None``
        Responsible for populating the state vector ``y`` (length `n_vars`)
        with the initial state of the reactor.

    ``update_state(self, y : double[:]) -> None``
        Responsible for setting the state of the reactor object from the
        values in the state vector ``y`` (length `n_vars`)

    ``update_surface_state(self, y : double[:]) -> None``
        Responsible for setting the state of surface phases in this reactor
        from the values in the state vector ``y``. The length of ``y`` is the
        total number of surface species in all surfaces.

    ``get_surface_initial_conditions(self, y : double[:]) -> None``
        Responsible for populating the state vector ``y`` with the initial
        state of each surface phase in this reactor. The length of ``y`` is the
        total number of surface species in all surfaces.

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

    ``eval_surfaces(LHS : double[:], RHS : double[:], sdot : double[:]) -> None``
        Responsible for calculating the ``LHS`` and ``RHS`` (length: total number of
        surface species in all surfaces) of the ODEs for surface species coverages,
        and the molar production rate of bulk phase species ``sdot`` (length: number
        of bulk phase species).

    ``component_name(i : int) -> string``
        Returns the name of the state vector component with index ``i``

    ``component_index(name: string) -> int``
        Returns the index of the state vector component named ``name``

    ``species_index(name : string) -> int``
        Returns the index of the species named ``name``, in either the bulk
        phase or a surface phase, relative to the start of the species terms in
        the state vector.
    """

    reactor_type = "ExtensibleReactor"

    delegatable_methods = {
        'initialize': ('initialize', 'void(double)'),
        'sync_state': ('syncState', 'void()'),
        'get_state': ('getState', 'void(double*)'),
        'update_state': ('updateState', 'void(double*)'),
        'update_surface_state': ('updateSurfaceState', 'void(double*)'),
        'get_surface_initial_conditions': ('getSurfaceInitialConditions', 'void(double*)'),
        'update_connected': ('updateConnected', 'void(bool)'),
        'eval': ('eval', 'void(double, double*, double*)'),
        'eval_walls': ('evalWalls', 'void(double)'),
        'eval_surfaces': ('evalSurfaces', 'void(double*,double*,double*)'),
        'component_name': ('componentName', 'string(size_t)'),
        'component_index': ('componentIndex', 'size_t(string)'),
        'species_index': ('speciesIndex', 'size_t(string)')
    }

    def __cinit__(self, *args, **kwargs):
        self.accessor = dynamic_cast[CxxReactorAccessorPtr](self.rbase)

    def __init__(self, *args, **kwargs):
        assign_delegates(self, dynamic_cast[CxxDelegatorPtr](self.rbase))
        super().__init__(*args, **kwargs)

    property n_vars:
        """
        Get/Set the number of state variables in the reactor.
        """
        def __get__(self):
            return self.reactor.neq()
        def __set__(self, n):
            self.accessor.setNEq(n)

    @property
    def expansion_rate(self):
        """
        Get/Set the net rate of volume change (for example, from moving walls) [m^3/s]

        .. versionadded:: 3.0
        """
        return self.accessor.expansionRate()

    @expansion_rate.setter
    def expansion_rate(self, vdot):
        self.accessor.setExpansionRate(vdot)

    @property
    def heat_rate(self):
        """
        Get/Set the net heat transfer rate (for example, through walls) [W]

        .. versionadded:: 3.0
        """
        return self.accessor.heatRate()

    @heat_rate.setter
    def heat_rate(self, qdot):
        self.accessor.setHeatRate(qdot)

    def restore_thermo_state(self):
        """
        Set the state of the thermo object to correspond to the state of the
        reactor.
        """
        self.accessor.restoreThermoState()

    def restore_surface_state(self, n):
        """
        Set the state of the thermo object for surface ``n`` to correspond to the
        state of that surface
        """
        self.accessor.restoreSurfaceState(n)


cdef class ExtensibleIdealGasReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasReactor` class.
    """
    reactor_type = "ExtensibleIdealGasReactor"


cdef class ExtensibleConstPressureReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `ConstPressureReactor` class.
    """
    reactor_type = "ExtensibleConstPressureReactor"


cdef class ExtensibleIdealGasConstPressureReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasConstPressureReactor` class.
    """
    reactor_type = "ExtensibleIdealGasConstPressureReactor"


cdef class ExtensibleMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `MoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleMoleReactor"


cdef class ExtensibleIdealGasMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasMoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleIdealGasMoleReactor"


cdef class ExtensibleConstPressureMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `ConstPressureMoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleConstPressureMoleReactor"


cdef class ExtensibleIdealGasConstPressureMoleReactor(ExtensibleReactor):
    """
    A variant of `ExtensibleReactor` where the base behavior corresponds to the
    `IdealGasConstPressureMoleReactor` class.

    .. versionadded:: 3.0
    """
    reactor_type = "ExtensibleIdealGasConstPressureMoleReactor"


cdef class ReactorSurface:
    """
    Represents a surface in contact with the contents of a reactor.

    :param kin:
        The `Kinetics` or `Interface` object representing reactions on this
        surface.
    :param r:
        The `Reactor` into which this surface should be installed.
    :param A:
        The area of the reacting surface [m^2]
    :param node_attr:
        Attributes to be passed to the ``node`` method invoked to draw this surface.
        See https://graphviz.org/docs/nodes/ for a list of all usable attributes.

    .. versionadded:: 3.1
       Added the ``node_attr`` parameter.
    """
    def __cinit__(self):
        self.surface = new CxxReactorSurface()

    def __dealloc__(self):
        del self.surface

    def __init__(self, kin=None, Reactor r=None, *, A=None, node_attr=None):
        if kin is not None:
            self.kinetics = kin
        if r is not None:
            self.install(r)
        if A is not None:
            self.area = A
        self.node_attr = node_attr or {'shape': 'underline'}

    def install(self, Reactor r):
        """
        Add this `ReactorSurface` to the specified `Reactor`
        """
        r._surfaces.append(self)
        r.reactor.addSurface(self.surface)
        self._reactor = r

    property area:
        """ Area on which reactions can occur [m^2] """
        def __get__(self):
            return self.surface.area()
        def __set__(self, A):
            self.surface.setArea(A)

    property kinetics:
        """
        The `InterfaceKinetics` object used for calculating reaction rates on
        this surface.
        """
        def __get__(self):
            self.surface.syncState()
            return self._kinetics
        def __set__(self, Kinetics k):
            self._kinetics = k
            self.surface.setKinetics(self._kinetics.kinetics)

    property coverages:
        """
        The fraction of sites covered by each surface species.
        """
        def __get__(self):
            if self._kinetics is None:
                raise CanteraError('No kinetics manager present')
            self.surface.syncState()
            return self._kinetics.coverages
        def __set__(self, coverages):
            if self._kinetics is None:
                raise CanteraError("Can't set coverages before assigning kinetics manager.")

            if isinstance(coverages, (dict, str, bytes)):
                self.surface.setCoverages(comp_map(coverages))
                return

            if len(coverages) != self._kinetics.n_species:
                raise ValueError('Incorrect number of site coverages specified')
            cdef np.ndarray[np.double_t, ndim=1] data = \
                    np.ascontiguousarray(coverages, dtype=np.double)
            self.surface.setCoverages(&data[0])

    @property
    def reactor(self):
        """
        Return the `Reactor` object the surface is connected to.

        .. versionadded:: 3.1
        """
        return self._reactor

    def draw(self, graph=None, *, graph_attr=None, node_attr=None, surface_edge_attr=None,
             print_state=False, **kwargs):
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
        :param kwargs:
            Additional keywords are passed on to ``~.drawnetwork.draw_reactor``.
        :return:
            ``graphviz.graphs.BaseGraph`` object with surface and connected
            reactor.

        .. versionadded:: 3.1
        """
        return draw_surface(self, graph, graph_attr, node_attr, surface_edge_attr,
                            print_state, **kwargs)

    def add_sensitivity_reaction(self, int m):
        """
        Specifies that the sensitivity of the state variables with respect to
        reaction ``m`` should be computed. ``m`` is the 0-based reaction index.
        The Surface must be installed on a reactor and part of a network first.
        """
        self.surface.addSensitivityReaction(m)

cdef class ReactorEdge:
    """
    Represents an edge in contact with the contents of a reactor.

    :param kin:
        The `Kinetics` or `Interface` object representing reactions on this
        surface.
    :param r:
        The `Reactor` into which this surface should be installed.
    :param L:
        The length of the reacting edge [m]
    :param node_attr:
        Attributes to be passed to the ``node`` method invoked to draw this surface.
        See https://graphviz.org/docs/nodes/ for a list of all usable attributes.

    .. versionadded:: 3.1????
       Added the ``node_attr`` parameter.
    """

    def __cinit__(self):
        self.edge = new CxxReactorEdge()

    def __dealloc__(self):
        del self.edge

    def __init__(self, kin=None, Reactor r=None, *, L=None):
        if kin is not None:
            self.kinetics = kin
        if r is not None:
            self.install(r)
        if L is not None:
            self.length = L
        #self.node_attr = node_attr or {'shape': 'underline'}

    def install(self, Reactor r):
        """
        Add this `ReactorEdge` to the specified `Reactor`
        """
        r._edges.append(self)
        r.reactor.addEdge(self.edge)
        self._reactor = r

    property length:
        """ Length on which reactions can occur [m] """
        def __get__(self):
            return self.edge.length()
        def __set__(self, L):
            self.edge.setLength(L)

    property kinetics:
        """
        The `InterfaceKinetics` object used for calculating reaction rates on
        this edge.
        """
        def __get__(self):
            self.edge.syncState()
            return self._kinetics
        def __set__(self, Kinetics k):
            self._kinetics = k
            self.edge.setKinetics(self._kinetics.kinetics)

    @property
    def reactor(self):
        """
        Return the `Reactor` object the edge is connected to.

        .. versionadded:: 3.1???
        """
        return self._reactor

"""
    def draw(self, graph=None, *, graph_attr=None, node_attr=None, surface_edge_attr=None,
             print_state=False, **kwargs):
        
        Draw the edge as a ``graphviz`` ``dot`` node connected to its reactor.
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
        :param kwargs:
            Additional keywords are passed on to ``~.drawnetwork.draw_reactor``.
        :return:
            ``graphviz.graphs.BaseGraph`` object with edge and connected
            reactor.

        .. versionadded:: 3.1
        
        return draw_edge(self, graph, graph_attr, node_attr, surface_edge_attr,
                            print_state, **kwargs)    
"""

cdef class WallBase:
    """
    Common base class for walls.
    """
    wall_type = "none"
    def __cinit__(self, *args, **kwargs):
        self._wall = newWall(stringify(self.wall_type))
        self.wall = self._wall.get()

    def __init__(self, left, right, *, name=None, A=None, K=None, U=None,
                 Q=None, velocity=None, edge_attr=None):
        """
        :param left:
            Reactor or reservoir on the left. Required.
        :param right:
            Reactor or reservoir on the right. Required.
        :param name:
            Name string. If omitted, the name is ``'Wall_n'``, where ``'n'``
            is an integer assigned in the order walls are created.
        :param A:
            Wall area [m^2]. Defaults to 1.0 m^2.
        :param K:
            Wall expansion rate parameter [m/s/Pa]. Defaults to 0.0.
        :param U:
            Overall heat transfer coefficient [W/m^2]. Defaults to 0.0
            (adiabatic wall).
        :param Q:
            Heat flux function :math:`q_0(t)` [W/m^2]. Optional. Default:
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

        self._install(left, right)
        if name is not None:
            self.name = name
        else:
            _reactor_counts['Wall'] += 1
            n = _reactor_counts['Wall']
            self.name = 'Wall_{0}'.format(n)

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

    def _install(self, ReactorBase left, ReactorBase right):
        """
        Install this Wall between two `Reactor` objects or between a
        `Reactor` and a `Reservoir`.
        """
        left._add_wall(self)
        right._add_wall(self)
        self.wall.install(deref(left.rbase), deref(right.rbase))
        # Keep references to prevent premature garbage collection
        self._left_reactor = left
        self._right_reactor = right

    property type:
        """The type of the wall."""
        def __get__(self):
            return pystr(self.wall.type())

    property area:
        """ The wall area [m^2]. """
        def __get__(self):
            return self.wall.area()
        def __set__(self, double value):
            self.wall.setArea(value)

    @property
    def left_reactor(self):
        """
        Return the `Reactor` or `Reservoir` object left of the wall.

        .. versionadded:: 3.1
        """
        return self._left_reactor

    @property
    def right_reactor(self):
        """
        Return the `Reactor` or `Reservoir` object right of the wall.

        .. versionadded:: 3.1
        """
        return self._right_reactor

    @property
    def expansion_rate(self):
        """
        Get the rate of volumetric change [m^3/s] associated with the wall at the
        current reactor network time. A positive value corresponds to the left-hand
        reactor volume increasing, and the right-hand reactor volume decreasing.

        .. versionadded:: 3.0
        """
        return self.wall.expansionRate()

    @property
    def heat_rate(self):
        """
        Get the total heat flux [W] through the wall  at the current reactor network
        time. A positive value corresponds to heat flowing from the left-hand reactor
        to the right-hand one.

        .. versionadded:: 3.0
        """
        return self.wall.heatRate()


    def draw(self, graph=None, *, graph_attr=None, node_attr=None, edge_attr=None,
             moving_wall_edge_attr=None, show_wall_velocity=True):
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


cdef class Wall(WallBase):
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
    wall_type = "Wall"

    property expansion_rate_coeff:
        """
        The coefficient *K* [m/s/Pa] that determines the velocity of the wall
        as a function of the pressure difference between the adjacent reactors.
        """
        def __get__(self):
            return (<CxxWall*>(self.wall)).getExpansionRateCoeff()
        def __set__(self, double val):
            (<CxxWall*>(self.wall)).setExpansionRateCoeff(val)

    property heat_transfer_coeff:
        """the overall heat transfer coefficient [W/m^2/K]"""
        def __get__(self):
            return (<CxxWall*>(self.wall)).getHeatTransferCoeff()
        def __set__(self, double value):
            (<CxxWall*>(self.wall)).setHeatTransferCoeff(value)

    property emissivity:
        """The emissivity (nondimensional)"""
        def __get__(self):
            return (<CxxWall*>(self.wall)).getEmissivity()
        def __set__(self, double value):
            (<CxxWall*>(self.wall)).setEmissivity(value)

    @property
    def velocity(self):
        """
        The wall velocity [m/s]. May be either set to a constant or an arbitrary
        function of time. See `Func1`.

        .. versionadded:: 3.0
        """
        return (<CxxWall*>(self.wall)).velocity()

    @velocity.setter
    def velocity(self, v):
        cdef Func1 f
        if isinstance(v, Func1):
            f = v
        else:
            f = Func1(v)

        self._velocity_func = f
        (<CxxWall*>(self.wall)).setVelocity(f.func)

    @property
    def heat_flux(self):
        """
        Heat flux [W/m^2] across the wall. May be either set to a constant or
        an arbitrary function of time. See `Func1`.

        .. versionadded:: 3.0
        """
        return (<CxxWall*>(self.wall)).heatFlux()

    @heat_flux.setter
    def heat_flux(self, q):
        cdef Func1 f
        if isinstance(q, Func1):
            f = q
        else:
            f = Func1(q)

        self._heat_flux_func = f
        (<CxxWall*>self.wall).setHeatFlux(f.func)


cdef class FlowDevice:
    """
    Base class for devices that allow flow between reactors.

    FlowDevice objects are assumed to be adiabatic, non-reactive, and have
    negligible internal volume, so that they are internally always in
    steady-state even if the upstream and downstream reactors are not. The
    fluid enthalpy, chemical composition, and mass flow rate are constant
    across a FlowDevice, and the pressure difference equals the difference in
    pressure between the upstream and downstream reactors.
    """
    flowdevice_type = "none"
    def __cinit__(self, *args, **kwargs):
        self._dev = newFlowDevice(stringify(self.flowdevice_type))
        self.dev = self._dev.get()

    def __init__(self, upstream, downstream, *, name=None, edge_attr=None):
        assert self.dev != NULL
        self._rate_func = None

        if name is not None:
            self.name = name
        else:
            _reactor_counts[self.__class__.__name__] += 1
            n = _reactor_counts[self.__class__.__name__]
            self.name = '{0}_{1}'.format(self.__class__.__name__, n)

        self.edge_attr = edge_attr or {}

        self._install(upstream, downstream)

    property type:
        """The type of the flow device."""
        def __get__(self):
            return pystr(self.dev.type())

    def _install(self, ReactorBase upstream, ReactorBase downstream):
        """
        Install the device between the ``upstream`` (source) and ``downstream``
        (destination) reactors or reservoirs.
        """
        upstream._add_outlet(self)
        downstream._add_inlet(self)
        self.dev.install(deref(upstream.rbase), deref(downstream.rbase))
        # Keep references to prevent premature garbage collection
        self._upstream = upstream
        self._downstream = downstream

    @property
    def upstream(self):
        """
        Return the `Reactor` or `Reservoir` object upstream of the flow device.

        .. versionadded:: 3.1
        """
        return self._upstream

    @property
    def downstream(self):
        """
        Return the `Reactor` or `Reservoir` object downstream of the flow device.

        .. versionadded:: 3.1
        """
        return self._downstream

    property mass_flow_rate:
        """
        Get the mass flow rate [kg/s] through this device at the current reactor
        network time.
        """
        def __get__(self):
            return self.dev.massFlowRate()

    @property
    def pressure_function(self):
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
    def pressure_function(self, k):
        cdef Func1 f
        if isinstance(k, Func1):
            f = k
        else:
            f = Func1(k)
        self._rate_func = f
        self.dev.setPressureFunction(f.func)

    @property
    def time_function(self):
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
    def time_function(self, k):
        cdef Func1 g
        if isinstance(k, Func1):
            g = k
        else:
            g = Func1(k)
        self._time_func = g
        self.dev.setTimeFunction(g.func)


    def draw(self, graph=None, *, graph_attr=None, node_attr=None, edge_attr=None):
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


cdef class MassFlowController(FlowDevice):
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
    flowdevice_type = "MassFlowController"

    def __init__(self, upstream, downstream, *, name=None, mdot=1., **kwargs):
        super().__init__(upstream, downstream, name=name, **kwargs)
        self.mass_flow_rate = mdot

    property mass_flow_coeff:
        r"""Set the mass flow rate [kg/s] through the mass flow controller
        as a constant, which may be modified by a function of time, see
        `time_function`.

        >>> mfc = MassFlowController(res1, reactor1)
        >>> mfc.mass_flow_coeff = 1e-4  # Set the flow rate to a constant
        >>> mfc.mass_flow_coeff  # Get the flow rate value
        """
        def __get__(self):
            return (<CxxMassFlowController*>self.dev).getMassFlowCoeff()
        def __set__(self, double value):
            (<CxxMassFlowController*>self.dev).setMassFlowCoeff(value)

    property mass_flow_rate:
        r"""
        Set the mass flow rate [kg/s] through this controller to be either
        a constant or an arbitrary function of time. See `Func1`, or get its
        current value.

        Note that depending on the argument type, this method either changes
        the property `mass_flow_coeff` or updates the `time_function` property.

        >>> mfc.mass_flow_rate = 0.3
        >>> mfc.mass_flow_rate = lambda t: 2.5 * exp(-10 * (t - 0.5)**2)
        """
        def __get__(self):
            return self.dev.massFlowRate()

        def __set__(self, m):
            if isinstance(m, _numbers.Real):
                (<CxxMassFlowController*>self.dev).setMassFlowRate(m)
            else:
                self.mass_flow_coeff = 1.
                self.time_function = m


cdef class Valve(FlowDevice):
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
    flowdevice_type = "Valve"

    def __init__(self, upstream, downstream, *, name=None, K=1., **kwargs):
        super().__init__(upstream, downstream, name=name, **kwargs)
        if isinstance(K, _numbers.Real):
            self.valve_coeff = K
        else:
            self.valve_coeff = 1.
            self.pressure_function = K

    property valve_coeff:
        r"""Set valve coefficient, that is, the proportionality constant between mass
        flow rate and pressure drop [kg/s/Pa].

        >>> v = Valve(res1, reactor1)
        >>> v.valve_coeff = 1e-4  # Set the value of K to a constant
        >>> v.valve_coeff  # Get the value of K
        """
        def __get__(self):
            return (<CxxValve*>self.dev).getValveCoeff()
        def __set__(self, double value):
            (<CxxValve*>self.dev).setValveCoeff(value)


cdef class PressureController(FlowDevice):
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
    flowdevice_type = "PressureController"

    def __init__(self, upstream, downstream, *,
            name=None, primary=None, K=1., **kwargs):
        if "master" in kwargs:
            warnings.warn(
                "PressureController: The 'master' keyword argument is deprecated; "
                "use 'primary' instead.", DeprecationWarning)
            primary = kwargs["master"]
        super().__init__(upstream, downstream, name=name, **kwargs)
        if primary is not None:
            self.primary = primary
        if isinstance(K, _numbers.Real):
            self.pressure_coeff = K
        else:
            self.pressure_coeff = 1.
            self.pressure_function = K

    property pressure_coeff:
        """
        Get/set the proportionality constant :math:`K_v` [kg/s/Pa] between the
        pressure drop and the mass flow rate.
        """
        def __get__(self):
            return (<CxxPressureController*>self.dev).getPressureCoeff()
        def __set__(self, double value):
            (<CxxPressureController*>self.dev).setPressureCoeff(value)

    @property
    def primary(self):
        """
        Primary `FlowDevice` used to compute this device's mass flow rate.

        .. versionadded:: 3.0
        """
        raise NotImplementedError("PressureController.primary")

    @primary.setter
    def primary(self, FlowDevice d):
        (<CxxPressureController*>self.dev).setPrimary(d.dev)


cdef class ReactorNet:
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
    def __init__(self, reactors=()):
        self._reactors = []  # prevents premature garbage collection
        for R in reactors:
            self.add_reactor(R)

    def add_reactor(self, Reactor r):
        """Add a reactor to the network."""
        self._reactors.append(r)
        self.net.addReactor(deref(r.reactor))

    def advance(self, double t, pybool apply_limit=True):
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

    def step(self):
        """
        Take a single internal step. The time/distance after taking the step is
        returned.
        """
        return self.net.step()

    def initialize(self):
        """
        Force initialization of the integrator after initial setup.
        """
        self.net.initialize()

    def reinitialize(self):
        """
        Reinitialize the integrator after making changing to the state of the
        system. Changes to Reactor contents will automatically trigger
        reinitialization.
        """
        self.net.reinitialize()

    @property
    def reactors(self):
        """
        List of all reactors that are part of the reactor network.

        .. versionadded:: 3.1
        """
        return self._reactors

    @property
    def time(self):
        """
        The current time [s], for reactor networks that are solved in the time domain.
        """
        return self.net.time()

    @property
    def distance(self):
        """
        The current distance[ m] along the length of the reactor network, for reactors
        that are solved as a function of space.
        """
        return self.net.distance()

    @property
    def initial_time(self):
        """
        The initial time of the integrator. When set, integration is restarted from this
        time using the current state as the initial condition. Default: 0.0 s.

        .. versionadded:: 3.0
        """
        return self.net.getInitialTime()

    @initial_time.setter
    def initial_time(self, double t):
        self.net.setInitialTime(t)

    property max_time_step:
        """
        Get/set the maximum time step *t* [s] that the integrator is
        allowed to use. The default value is set to zero, so that no time
        step maximum is used.
        """
        def __get__(self):
            return self.net.maxTimeStep()

        def __set__(self, double t):
            self.net.setMaxTimeStep(t)

    property max_err_test_fails:
        """
        The maximum number of error test failures permitted by the CVODES integrator
        in a single step. The default is 10.
        """
        def __set__(self, n):
            self.net.setMaxErrTestFails(n)

    @property
    def max_nonlinear_iterations(self):
        """
        Get/Set the maximum number of nonlinear solver iterations permitted by the
        SUNDIALS solver in one solve attempt. The default value is 4.
        """
        return self.net.integrator().maxNonlinIterations()

    @max_nonlinear_iterations.setter
    def max_nonlinear_iterations(self, int n):
        self.net.integrator().setMaxNonlinIterations(n)

    @property
    def max_nonlinear_convergence_failures(self):
        """
        Get/Set the maximum number of nonlinear solver convergence failures permitted in
        one step of the SUNDIALS integrator. The default value is 10.
        """
        return self.net.integrator().maxNonlinConvFailures()

    @max_nonlinear_convergence_failures.setter
    def max_nonlinear_convergence_failures(self, int n):
        self.net.integrator().setMaxNonlinConvFailures(n)

    @property
    def include_algebraic_in_error_test(self):
        """
        Get/Set whether to include algebraic variables in the in the local error test.
        Applicable only to DAE systems. The default is `True`.
        """
        return self.net.integrator().algebraicInErrorTest()

    @include_algebraic_in_error_test.setter
    def include_algebraic_in_error_test(self, pybool yesno):
        self.net.integrator().includeAlgebraicInErrorTest(yesno)

    @property
    def max_order(self):
        """
        Get/Set the maximum order of the linear multistep method. The default value and
        maximum is 5.
        """
        return self.net.integrator().maxOrder()

    @max_order.setter
    def max_order(self, int n):
        self.net.integrator().setMaxOrder(n)

    property max_steps:
        """
        The maximum number of internal integration steps that CVODES
        is allowed to take before reaching the next output point.
        """
        def __set__(self, nsteps):
            self.net.setMaxSteps(nsteps)
        def __get__(self):
            return self.net.maxSteps()

    property rtol:
        """
        The relative error tolerance used while integrating the reactor
        equations.
        """
        def __get__(self):
            return self.net.rtol()
        def __set__(self, tol):
            self.net.setTolerances(tol, -1)

    property atol:
        """
        The absolute error tolerance used while integrating the reactor
        equations.
        """
        def __get__(self):
            return self.net.atol()
        def __set__(self, tol):
            self.net.setTolerances(-1, tol)

    property rtol_sensitivity:
        """
        The relative error tolerance for sensitivity analysis.
        """
        def __get__(self):
            return self.net.rtolSensitivity()
        def __set__(self, tol):
            self.net.setSensitivityTolerances(tol, -1)

    property atol_sensitivity:
        """
        The absolute error tolerance for sensitivity analysis.
        """
        def __get__(self):
            return self.net.atolSensitivity()
        def __set__(self, tol):
            self.net.setSensitivityTolerances(-1, tol)

    property verbose:
        """
        If `True`, verbose debug information will be printed during
        integration. The default is `False`.
        """
        def __get__(self):
            return pybool(self.net.verbose())
        def __set__(self, pybool v):
            self.net.setVerbose(v)

    def global_component_index(self, name, int reactor):
        """
        Returns the index of a component named ``name`` of a reactor with index
        ``reactor`` within the global state vector. That is, this determines the
        absolute index of the component, where ``reactor`` is the index of the
        reactor that holds the component. ``name`` is either a species name or the
        name of a reactor state variable, for example, ``'int_energy'``, ``'temperature'``, etc.
        depending on the reactor's equations.
        """
        return self.net.globalComponentIndex(stringify(name), reactor)

    def component_name(self, int i):
        """
        Return the name of the i-th component of the global state vector. The
        name returned includes both the name of the reactor and the specific
        component, for example `'reactor1: CH4'`.
        """
        return pystr(self.net.componentName(i))

    def sensitivity(self, component, int p, int r=0):
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

    def sensitivities(self):
        r"""
        Returns the sensitivities of all of the solution variables with respect
        to all of the registered parameters. The normalized sensitivity
        coefficient :math:`S_{ki}` of the solution variable :math:`y_k` with
        respect to sensitivity parameter :math:`p_i` is defined as:

        .. math:: S_{ki} = \frac{p_i}{y_k} \frac{\partial y_k}{\partial p_i}

        For reaction sensitivities, the parameter is a multiplier on the forward
        rate constant (and implicitly on the reverse rate constant for
        reversible reactions).

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
        cdef np.ndarray[np.double_t, ndim=2] data = \
                np.empty((self.n_vars, self.n_sensitivity_params))
        cdef int p,k
        for p in range(self.n_sensitivity_params):
            for k in range(self.n_vars):
                data[k,p] = self.net.sensitivity(k,p)
        return data

    def sensitivity_parameter_name(self, int p):
        """
        Name of the sensitivity parameter with index ``p``.
        """
        return pystr(self.net.sensitivityParameterName(p))

    property n_sensitivity_params:
        """
        The number of registered sensitivity parameters.
        """
        def __get__(self):
            return self.net.nparams()

    property n_vars:
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
        def __get__(self):
            return self.net.neq()

    def get_state(self):
        """
        Get the combined state vector of the reactor network.

        The combined state vector consists of the concatenated state vectors of
        all entities contained.
        """
        if not self.n_vars:
            raise CanteraError('ReactorNet empty or not initialized.')
        cdef np.ndarray[np.double_t, ndim=1] y = np.zeros(self.n_vars)
        self.net.getState(&y[0])
        return y

    def get_derivative(self, k):
        """
        Get the k-th derivative of the state vector of the reactor network with respect
        to the independent integrator variable (time/distance).
        """
        if not self.n_vars:
            raise CanteraError('ReactorNet empty or not initialized.')
        cdef np.ndarray[np.double_t, ndim = 1] dky = np.zeros(self.n_vars)
        self.net.getDerivative(k, & dky[0])
        return dky

    property advance_limits:
        """
        Get or set absolute limits for state changes during `ReactorNet.advance`
        (positive values are considered; negative values disable a previously
        set advance limit for a solution component). Note that limits are
        disabled by default (with individual values set to -1.).
        """
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] limits = np.empty(self.n_vars)
            self.net.getAdvanceLimits(&limits[0])
            return limits

        def __set__(self, limits):
            if limits is None:
                limits = -1. * np.ones([self.n_vars])
            elif len(limits) != self.n_vars:
                raise ValueError('array must be of length n_vars')

            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(limits, dtype=np.double)
            self.net.setAdvanceLimits(&data[0])

    def advance_to_steady_state(self, int max_steps=10000,
                                double residual_threshold=0., double atol=0.,
                                pybool return_residuals=False):
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

    def __reduce__(self):
        raise NotImplementedError('ReactorNet object is not picklable')

    def __copy__(self):
        raise NotImplementedError('ReactorNet object is not copyable')

    property preconditioner:
        """Preconditioner associated with integrator"""
        def __set__(self, PreconditionerBase precon):
            # set preconditioner
            self.net.setPreconditioner(precon.pbase)
            # set problem type as default of preconditioner
            self.linear_solver_type = precon.precon_linear_solver_type

    property linear_solver_type:
        """
            The type of linear solver used in integration.

            Options for this property include:

              - `"DENSE"`
              - `"GMRES"`
              - `"BAND"`
              - `"DIAG"`

        """
        def __set__(self, linear_solver_type):
            self.net.setLinearSolverType(stringify(linear_solver_type))

        def __get__(self):
            return pystr(self.net.linearSolverType())


    property solver_stats:
        """ODE solver stats from integrator"""
        def __get__(self):
            cdef CxxAnyMap stats
            stats = self.net.solverStats()
            return anymap_to_py(stats)

    property derivative_settings:
        """
        Apply derivative settings to all reactors in the network.
        See also `Kinetics.derivative_settings`.
        """
        def __set__(self, settings):
            self.net.setDerivativeSettings(py_to_anymap(settings))

    def draw(self, *, graph_attr=None, node_attr=None, edge_attr=None,
             heat_flow_attr=None, mass_flow_attr=None, moving_wall_edge_attr=None,
             surface_edge_attr=None, print_state=False, show_wall_velocity=True,
             **kwargs):
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
        :param print_state:
            Whether state information of the reactors is printed into each node.
            Defaults to ``False``.
        :param kwargs:
            Additional keywords are passed on to each call of
            `~.drawnetwork.draw_reactor`, `~.drawnetwork.draw_surface`,
            `~.drawnetwork.draw_flow_controllers`, and `~.drawnetwork.draw_walls`.
        :return:
            ``graphviz.graphs.BaseGraph`` object with reactor net.

        .. versionadded:: 3.1
        """
        return draw_reactor_net(self, graph_attr, node_attr, edge_attr, heat_flow_attr,
                 mass_flow_attr, moving_wall_edge_attr, surface_edge_attr, print_state,
                 **kwargs)
