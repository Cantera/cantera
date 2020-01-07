# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import warnings
from collections import defaultdict as _defaultdict
import numbers as _numbers

_reactor_counts = _defaultdict(int)

# Need a pure-python class to store weakrefs to
class _WeakrefProxy:
    pass

cdef class ReactorBase:
    """
    Common base class for reactors and reservoirs.
    """
    reactor_type = "None"
    def __cinit__(self, *args, **kwargs):
        self.rbase = newReactor(stringify(self.reactor_type))

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, ThermoPhase contents=None, name=None, *, volume=None):
        self._weakref_proxy = _WeakrefProxy()
        self._inlets = []
        self._outlets = []
        self._walls = []
        if isinstance(contents, ThermoPhase):
            self.insert(contents)

        if name is not None:
            self.name = name
        else:
            _reactor_counts[self.reactor_type] += 1
            n = _reactor_counts[self.reactor_type]
            self.name = '{0}_{1}'.format(self.reactor_type, n)

        if volume is not None:
            self.volume = volume

    def __dealloc__(self):
        del self.rbase

    def insert(self, _SolutionBase solution):
        """
        Set *solution* to be the object used to compute thermodynamic
        properties and kinetic rates for this reactor.
        """
        self._thermo = solution
        self._thermo._references[self._weakref_proxy] = True
        self.rbase.setThermoMgr(deref(solution.thermo))

    property type:
        """The type of the reactor."""
        def __get__(self):
            return pystr(self.rbase.typeStr())

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

    def _add_inlet(self, inlet):
        """
        Store a reference to *inlet* to prevent it from being prematurely
        garbage collected.
        """
        self._inlets.append(inlet)

    def _add_outlet(self, outlet):
        """
        Store a reference to *outlet* to prevent it from being prematurely
        garbage collected.
        """
        self._outlets.append(outlet)

    def _add_wall(self, wall):
        """
        Store a reference to *wall* to prevent it from being prematurely
        garbage collected.
        """
        self._walls.append(wall)

    def __reduce__(self):
        raise NotImplementedError('Reactor object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Reactor object is not copyable')


cdef class Reactor(ReactorBase):
    """
    A homogeneous zero-dimensional reactor. By default, they are closed
    (no inlets or outlets), have fixed volume, and have adiabatic,
    chemically-inert walls. These properties may all be changed by adding
    appropriate components, e.g. `Wall`, `MassFlowController` and `Valve`.
    """
    reactor_type = "Reactor"

    def __cinit__(self, *args, **kwargs):
        self.reactor = <CxxReactor*>(self.rbase)

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, contents=None, *, name=None, energy='on', **kwargs):
        """
        :param contents:
            Reactor contents. If not specified, the reactor is initially empty.
            In this case, call `insert` to specify the contents.
        :param name:
            Used only to identify this reactor in output. If not specified,
            defaults to ``'Reactor_n'``, where *n* is an integer assigned in
            the order `Reactor` objects are created.
        :param energy:
            Set to ``'on'`` or ``'off'``. If set to ``'off'``, the energy
            equation is not solved, and the temperature is held at its
            initial value..

        Some examples showing how to create :class:`Reactor` objects are
        shown below.

        >>> gas = Solution('gri30.xml')
        >>> r1 = Reactor(gas)

        This is equivalent to:

        >>> r1 = Reactor()
        >>> r1.insert(gas)

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

    def insert(self, _SolutionBase solution):
        ReactorBase.insert(self, solution)
        self._kinetics = solution
        if solution.kinetics != NULL:
            self.reactor.setKineticsMgr(deref(solution.kinetics))

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
        *True* when the reactor composition is allowed to change due to
        chemical reactions in this reactor. When this is *False*, the
        reactor composition is held constant.
        """
        def __get__(self):
            return self.reactor.chemistryEnabled()

        def __set__(self, pybool value):
            self.reactor.setChemistry(value)

    property energy_enabled:
        """
        *True* when the energy equation is being solved for this reactor.
        When this is *False*, the reactor temperature is held constant.
        """
        def __get__(self):
            return self.reactor.energyEnabled()

        def __set__(self, pybool value):
            self.reactor.setEnergy(int(value))

    def add_sensitivity_reaction(self, m):
        """
        Specifies that the sensitivity of the state variables with respect to
        reaction *m* should be computed. *m* is the 0-based reaction index.
        The reactor must be part of a network first. Specifying the same
        reaction more than one time raises an exception.
        """
        self.reactor.addSensitivityReaction(m)

    def add_sensitivity_species_enthalpy(self, k):
        """
        Specifies that the sensitivity of the state variables with respect to
        species *k* should be computed. The reactor must be part of a network
        first.
        """
        self.reactor.addSensitivitySpeciesEnthalpy(self.thermo.species_index(k))

    def component_index(self, name):
        """
        Returns the index of the component named *name* in the system. This
        determines the (relative) index of the component in the vector of
        sensitivity coefficients. *name* is either a species name or the name of
        a reactor state variable, e.g. 'int_energy', 'temperature', depending on
        the reactor's equations.
        """

        k = self.reactor.componentIndex(stringify(name))
        if k == CxxNpos:
            raise IndexError('No such component: {!r}'.format(name))
        return k

    def component_name(self, int i):
        """
        Returns the name of the component with index *i* within the array of
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

        The order of the variables (i.e. rows) is:

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

    def set_advance_limit(self, name, limit):
        """
        Limit absolute change of component *name* during `ReactorNet.advance`.
        (positive *limit* values are considered; negative values disable a
        previously set advance limit for a solution component). Note that
        limits are disabled by default (with individual values set to -1.).
        """
        if limit is None:
            limit = -1.
        self.reactor.setAdvanceLimit(stringify(name), limit)


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


cdef class IdealGasReactor(Reactor):
    """ A constant volume, zero-dimensional reactor for ideal gas mixtures. """
    reactor_type = "IdealGasReactor"


cdef class IdealGasConstPressureReactor(Reactor):
    """
    A homogeneous, constant pressure, zero-dimensional reactor for ideal gas
    mixtures. The volume of the reactor changes as a function of time in order
    to keep the pressure constant.
    """
    reactor_type = "IdealGasConstPressureReactor"


cdef class FlowReactor(Reactor):
    """
    A steady-state plug flow reactor with constant cross sectional area.
    Time integration follows a fluid element along the length of the reactor.
    The reactor is assumed to be frictionless and adiabatic.
    """
    reactor_type = "FlowReactor"

    property mass_flow_rate:
        """ Mass flow rate per unit area [kg/m^2*s] """
        def __set__(self, double value):
            (<CxxFlowReactor*>self.reactor).setMassFlowRate(value)

    property speed:
        """ Speed [m/s] of the flow in the reactor at the current position """
        def __get__(self):
            return (<CxxFlowReactor*>self.reactor).speed()

    property distance:
        """ The distance of the fluid element from the inlet of the reactor."""
        def __get__(self):
            return (<CxxFlowReactor*>self.reactor).distance()


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
    """
    def __cinit__(self):
        self.surface = new CxxReactorSurface()

    def __dealloc__(self):
        del self.surface

    def __init__(self, kin=None, Reactor r=None, *, A=None):
        if kin is not None:
            self.kinetics = kin
        if r is not None:
            self.install(r)
        if A is not None:
            self.area = A

    def install(self, Reactor r):
        r.reactor.addSurface(self.surface)

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
            self.surface.syncCoverages()
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

    def add_sensitivity_reaction(self, int m):
        """
        Specifies that the sensitivity of the state variables with respect to
        reaction *m* should be computed. *m* is the 0-based reaction index.
        The Surface must be installed on a reactor and part of a network first.
        """
        self.surface.addSensitivityReaction(m)


cdef class WallBase:
    """
    Common base class for walls.
    """
    wall_type = "None"
    def __cinit__(self, *args, **kwargs):
        self.wall = newWall(stringify(self.wall_type))

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, left, right, *, name=None, A=None, K=None, U=None,
                 Q=None, velocity=None):
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
        """
        self.left_surface = WallSurface(self, 0)
        self.right_surface = WallSurface(self, 1)

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
            self.set_heat_flux(Q)
        if velocity is not None:
            self.set_velocity(velocity)

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

    property left:
        """ The left surface of this wall. """
        def __get__(self):
            return self.left_surface

    property right:
        """ The right surface of this wall. """
        def __get__(self):
            return self.right_surface

    property area:
        """ The wall area [m^2]. """
        def __get__(self):
            return self.wall.area()
        def __set__(self, double value):
            self.wall.setArea(value)

    def vdot(self, double t):
        """
        The rate of volumetric change [m^3/s] associated with the wall
        at time *t*. A positive value corresponds to the left-hand reactor
        volume increasing, and the right-hand reactor volume decreasing.
        """
        return self.wall.vdot(t)

    def qdot(self, double t):
        """
        Total heat flux [W] through the wall at time *t*. A positive value
        corresponds to heat flowing from the left-hand reactor to the
        right-hand one.
        """
        return self.wall.Q(t)


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

    def set_velocity(self, v):
        """
        The wall velocity [m/s]. May be either a constant or an arbitrary
        function of time. See `Func1`.
        """
        cdef Func1 f
        if isinstance(v, Func1):
            f = v
        else:
            f = Func1(v)

        self._velocity_func = f
        (<CxxWall*>(self.wall)).setVelocity(f.func)

    def set_heat_flux(self, q):
        """
        Heat flux [W/m^2] across the wall. May be either a constant or
        an arbitrary function of time. See `Func1`.
        """
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
    flowdevice_type = "None"
    def __cinit__(self, *args, **kwargs):
        self.dev = newFlowDevice(stringify(self.flowdevice_type))

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, upstream, downstream, *, name=None):
        assert self.dev != NULL
        self._rate_func = None

        if name is not None:
            self.name = name
        else:
            _reactor_counts[self.__class__.__name__] += 1
            n = _reactor_counts[self.__class__.__name__]
            self.name = '{0}_{1}'.format(self.__class__.__name__, n)

        self._install(upstream, downstream)

    def __dealloc__(self):
        del self.dev

    property type:
        """The type of the flow device."""
        def __get__(self):
            return pystr(self.dev.typeStr())

    def _install(self, ReactorBase upstream, ReactorBase downstream):
        """
        Install the device between the *upstream* (source) and *downstream*
        (destination) reactors or reservoirs.
        """
        upstream._add_outlet(self)
        downstream._add_inlet(self)
        self.dev.install(deref(upstream.rbase), deref(downstream.rbase))
        # Keep references to prevent premature garbage collection
        self._upstream = upstream
        self._downstream = downstream

    def mdot(self, double t):
        """
        The mass flow rate [kg/s] through this device at time *t* [s].
        """
        return self.dev.massFlowRate(t)

    def set_pressure_function(self, k):
        r"""
        Set the relationship between mass flow rate and the pressure drop across a
        flow device. The mass flow rate [kg/s] is calculated given the pressure
        drop [Pa] and a coefficient set by a flow device specific function.
        The calculation of mass flow rate depends to the flow device.

        >>> F = FlowDevice(res1, reactor1)
        >>> F.set_pressure_function(lambda dP: dP**2)

        where FlowDevice is either a Valve or PressureController object.
        """
        cdef Func1 f
        if isinstance(k, Func1):
            f = k
        else:
            f = Func1(k)
        self._rate_func = f
        self.dev.setPressureFunction(f.func)

    def set_time_function(self, k):
        r"""
        Set the time dependence of a flow device. The mass flow rate [kg/s] is
        calculated for a flow device, and multiplied by a function of time.
        The calculation of mass flow rate depends to the flow device.

        >>> F = FlowDevice(res1, reactor1)
        >>> F.set_time_function(lambda t: exp(-10 * (t - 0.5)**2))

        where FlowDevice is either a Valve or MassFlowController object.
        """
        cdef Func1 g
        if isinstance(k, Func1):
            g = k
        else:
            g = Func1(k)
        self._time_func = g
        self.dev.setTimeFunction(g.func)


cdef class MassFlowController(FlowDevice):
    r"""
    A mass flow controller maintains a specified mass
    flow rate independent of upstream and downstream conditions. The equation
    used to compute the mass flow rate is

    .. math:: \dot m = \max(\dot m_0*g(t), 0.),

    where :math:`\dot m_0` is a constant value and :math:`g(t)` is a function of
    time. Both :math:`\dot m_0` and :math:`g(t)` can be set individually by
    the property `mass_flow_coeff` and the method `set_time_function`,
    respectively. The property `mass_flow_rate` combines the former
    into a single interface. Note that if :math:`\dot m_0*g(t) < 0`, the mass flow
    rate will be set to zero, since reversal of the flow direction is not allowed.

    Unlike a real mass flow controller, a MassFlowController object will
    maintain the flow even if the downstream pressure is greater than the
    upstream pressure. This allows simple implementation of loops, in which
    exhaust gas from a reactor is fed back into it through an inlet. But note
    that this capability should be used with caution, since no account is
    taken of the work required to do this.
    """
    flowdevice_type = "MassFlowController"

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, upstream, downstream, *, name=None, mdot=1.):
        super().__init__(upstream, downstream, name=name)
        self.mass_flow_rate = mdot

    property mass_flow_coeff:
        r"""Set the mass flow rate [kg/s] through the mass flow controller
        as a constant, which may be modified by a function of time, see
        `set_time_function`.

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
        a constant or an arbitrary function of time. See `Func1`.

        Note that depending on the argument type, this method either changes
        the property `mass_flow_coeff` or calls the `set_time_function` method.

        >>> mfc.mass_flow_rate = 0.3
        >>> mfc.mass_flow_rate = lambda t: 2.5 * exp(-10 * (t - 0.5)**2)
        """
        def __set__(self, m):
            if isinstance(m, _numbers.Real):
                (<CxxMassFlowController*>self.dev).setMassFlowRate(m)
            else:
                self.mass_flow_coeff = 1.
                self.set_time_function(m)

    def set_mass_flow_rate(self, m):
        r"""
        Set the mass flow rate [kg/s] through this controller to be either
        a constant or an arbitrary function of time. See `Func1`.

        Note that depending on the argument type, this method either changes
        the property `mass_flow_coeff` or calls the `set_time_function` method.

        >>> mfc.set_mass_flow_rate(0.3)
        >>> mfc.set_mass_flow_rate(lambda t: 2.5 * exp(-10 * (t - 0.5)**2))

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Replaced by property `mass_flow_rate`.
        """
        warnings.warn("To be removed after Cantera 2.5. "
                      "Replaced by property 'mass_flow_rate'", DeprecationWarning)

        self.mass_flow_rate = m


cdef class Valve(FlowDevice):
    r"""
    In Cantera, a `Valve` is a flow device with mass flow rate that is a
    function of the pressure drop across it. The default behavior is linear:

    .. math:: \dot m = K_v*(P_1 - P_2)

    where :math:`K_v` is a constant set by the `set_valve_coeff` method.
    Note that :math:`P_1` must be greater than :math:`P_2`; otherwise,
    :math:`\dot m = 0`. However, an arbitrary function can also be specified,
    such that

    .. math:: \dot m = K_v*f(P_1 - P_2)

    where :math:`f` is the arbitrary function that multiplies :math:`K_v` given
    a single argument, the pressure differential. Further, a valve opening function
    :math:`g` may be specified using the method `set_time_function`, such that

    .. math:: \dot m = K_v*g(t)*f(P_1 - P_2)

    See the documentation for the `valve_coeff` property as well as the
    `set_pressure_function` and `set_time_function` methods for examples. Note that
    it is never possible for the flow to reverse and go from the downstream to the
    upstream reactor/reservoir through a line containing a `Valve` object.

    `Valve` objects are often used between an upstream reactor and a
    downstream reactor or reservoir to maintain them both at nearly the same
    pressure. By setting the constant :math:`K_v` to a sufficiently large
    value, very small pressure differences will result in flow between the
    reactors that counteracts the pressure difference.
    """
    flowdevice_type = "Valve"

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, upstream, downstream, *, name=None, K=1.):
        super().__init__(upstream, downstream, name=name)
        if isinstance(K, _numbers.Real):
            self.valve_coeff = K
        else:
            self.valve_coeff = 1.
            self.set_pressure_function(K)

    property valve_coeff:
        r"""Set valve coefficient, i.e. the proportionality constant between mass
        flow rate and pressure drop [kg/s/Pa].

        >>> V = Valve(res1, reactor1)
        >>> V.valve_coeff = 1e-4  # Set the value of K to a constant
        >>> V.valve_coeff  # Get the value of K
        """
        def __get__(self):
            return (<CxxValve*>self.dev).getValveCoeff()
        def __set__(self, double value):
            (<CxxValve*>self.dev).setValveCoeff(value)

    def set_valve_function(self, k):
        r"""
        Set the relationship between mass flow rate and the pressure drop across the
        valve. The mass flow rate [kg/s] is calculated given the pressure drop [Pa].

        >>> V = Valve(res1, reactor1)
        >>> V.set_valve_function(lambda dP: (1e-5 * dP)**2)

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Renamed to `set_pressure_function`.
        """
        warnings.warn("To be removed after Cantera 2.5. "
                      "Renamed to 'set_pressure_function' instead", DeprecationWarning)

        self.set_pressure_function(k)

    def set_valve_coeff(self, k):
        """
        Set the relationship between mass flow rate and the pressure drop across
        the valve. If a number is given, it is the proportionality constant
        [kg/s/Pa]. If a function is given, it should compute the mass flow
        rate [kg/s] given the pressure drop [Pa].

        >>> V = Valve(res1, reactor1)
        >>> V.set_valve_coeff(1e-4)  # Set the value of K to a constant
        >>> V.set_valve_coeff(lambda dP: (1e-5 * dP)**2)  # Set to a function

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Functionality is now handled by property `valve_coeff` and
             `set_pressure_function`.
        """
        warnings.warn("To be removed after Cantera 2.5. "
                      "Use property 'valve_coeff' and/or function "
                      "'set_pressure_function' instead.", DeprecationWarning)

        if isinstance(k, _numbers.Real):
            self.valve_coeff = k
        else:
            self.valve_coeff = 1.
            self.set_pressure_function(k)


cdef class PressureController(FlowDevice):
    r"""
    A PressureController is designed to be used in conjunction with another
    'master' flow controller, typically a `MassFlowController`. The master
    flow controller is installed on the inlet of the reactor, and the
    corresponding `PressureController` is installed on the outlet of the
    reactor. The `PressureController` mass flow rate is equal to the master
    mass flow rate, plus a small correction dependent on the pressure
    difference:

    .. math:: \dot m = \dot m_{\rm master} + K_v(P_1 - P_2).

    As an alternative, an arbitrary function of pressure differential can be
    specified using the method `set_pressure_function`, such that

    .. math:: \dot m = \dot m_{\rm master} + K_v*f(P_1 - P_2)

    where :math:`f` is the arbitrary function of a single argument.
    """
    flowdevice_type = "PressureController"

    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, upstream, downstream, *, name=None, master=None, K=1.):
        super().__init__(upstream, downstream, name=name)
        if master is not None:
            self.set_master(master)
        if isinstance(K, _numbers.Real):
            self.pressure_coeff = K
        else:
            self.pressure_coeff = 1.
            self.set_pressure_function(K)

    property pressure_coeff:
        """
        Get/set the proportionality constant :math:`K_v` [kg/s/Pa] between the
        pressure drop and the mass flow rate.
        """
        def __get__(self):
            return (<CxxPressureController*>self.dev).getPressureCoeff()
        def __set__(self, double value):
            (<CxxPressureController*>self.dev).setPressureCoeff(value)

    def set_pressure_coeff(self, double k):
        """
        Set the proportionality constant :math:`K_v` [kg/s/Pa] between the pressure
        drop and the mass flow rate.

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Replaced by property `pressure_coeff`.
        """
        warnings.warn("To be removed after Cantera 2.5. "
                      "Use property 'pressure_coeff' instead", DeprecationWarning)
        (<CxxPressureController*>self.dev).setPressureCoeff(k)

    def set_master(self, FlowDevice d):
        """
        Set the "master" `FlowDevice` used to compute this device's mass flow
        rate.
        """
        (<CxxPressureController*>self.dev).setMaster(d.dev)


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
        Advance the state of the reactor network in time from the current time
        towards time *t* [s], taking as many integrator timesteps as necessary.
        If *apply_limit* is true and an advance limit is specified, the reactor
        state at the end of the timestep is estimated prior to advancing. If
        the difference exceed limits, the end time is reduced by half until
        the projected end state remains within specified limits.
        Returns the time reached at the end of integration.
        """
        return self.net.advance(t, apply_limit)

    def step(self):
        """
        Take a single internal time step. The time after taking the step is
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

    property time:
        """The current time [s]."""
        def __get__(self):
            return self.net.time()

    def set_initial_time(self, double t):
        """
        Set the initial time. Restarts integration from this time using the
        current state as the initial condition. Default: 0.0 s.
        """
        self.net.setInitialTime(t)

    property max_time_step:
        """
        Get/set the maximum time step *t* [s] that the integrator is
        allowed to use. The default value is set to zero, i.e. no time
        step maximum is used.
        """
        def __get__(self):
            return self.net.maxTimeStep()

        def __set__(self, double t):
            self.net.setMaxTimeStep(t)

    def set_max_time_step(self, double t):
        """
        Set the maximum time step *t* [s] that the integrator is allowed
        to use.

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Replaced by property `max_time_step`.
        """
        warnings.warn("To be removed after Cantera 2.5. "
                      "Use property 'max_time_step' instead", DeprecationWarning)
        self.net.setMaxTimeStep(t)

    property max_err_test_fails:
        """
        The maximum number of error test failures permitted by the CVODES
        integrator in a single time step.
        """
        def __set__(self, n):
            self.net.setMaxErrTestFails(n)

    property max_steps:
        """
        The maximum number of internal integration time-steps that CVODES
        is allowed to take before reaching the next output time.
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
        If *True*, verbose debug information will be printed during
        integration. The default is *False*.
        """
        def __get__(self):
            return pybool(self.net.verbose())
        def __set__(self, pybool v):
            self.net.setVerbose(v)

    def global_component_index(self, name, int reactor):
        """
        Returns the index of a component named *name* of a reactor with index
        *reactor* within the global state vector. I.e. this determines the
        (absolute) index of the component, where *reactor* is the index of the
        reactor that holds the component. *name* is either a species name or the
        name of a reactor state variable, e.g. 'int_energy', 'temperature', etc.
        depending on the reactor's equations.
        """
        return self.net.globalComponentIndex(stringify(name), reactor)

    def component_name(self, int i):
        """
        Return the name of the i-th component of the global state vector. The
        name returned includes both the name of the reactor and the specific
        component, e.g. `'reactor1: CH4'`.
        """
        return pystr(self.net.componentName(i))

    def sensitivity(self, component, int p, int r=0):
        """
        Returns the sensitivity of the solution variable *component* in
        reactor *r* with respect to the parameter *p*. *component* can be a
        string or an integer. See `component_index` and `sensitivities` to
        determine the integer index for the variables and the definition of the
        resulting sensitivity coefficient. If it is not given, *r* defaults to
        the first reactor. Returns an empty array until the first time step is
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
        n_sensitivity_params)*, unless no timesteps have been taken, in which
        case the shape is *(0, n_sensitivity_params)*. The order of the
        variables (i.e., rows) is:

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
        Name of the sensitivity parameter with index *p*.
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
        Get the k-th time derivative of the state vector of the reactor network.
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
            atol = self.rtol
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
