from collections import defaultdict
import numbers

reactor_counts = defaultdict(int)

cdef class ReactorBase:
    """
    Common base class for reactors and reservoirs.
    """
    reactor_type = "None"
    def __cinit__(self, *args, **kwargs):
        self.rbase = newReactor(stringify(self.reactor_type))

    def __init__(self, ThermoPhase contents=None, name=None, **kwargs):
        self._inlets = []
        self._outlets = []
        self._walls = []
        if isinstance(contents, ThermoPhase):
            self.insert(contents)

        if name is not None:
            self.name = name
        else:
            reactor_counts[self.reactor_type] += 1
            n = reactor_counts[self.reactor_type]
            self.name = '{0}_{1}'.format(self.reactor_type, n)

    def __dealloc__(self):
        del self.rbase

    def insert(self, _SolutionBase solution):
        """
        Set *solution* to be the object used to compute thermodynamic
        properties and kinetic rates for this reactor.
        """
        self._thermo = solution
        self.rbase.setThermoMgr(deref(solution.thermo))

    property name:
        """The name of the reactor."""
        def __get__(self):
            return pystr(self.rbase.name())

        def __set__(self, name):
            self.rbase.setName(stringify(name))

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
        super().__init__(contents, **kwargs)

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
        self.reactor.addSensitivityReaction(m)

    def component_index(self, name):
        """
        Returns the index of the component named *name* in the system. This
        determines the (relative) index of the component in the vector of
        sensitivity coefficients. *name* is either a species name or the name
        of a reactor state variable, e.g. 'U', 'T', depending on the reactor's
        equations.
        """

        k = self.reactor.componentIndex(stringify(name))
        if k == CxxNpos:
            raise IndexError('No such component: {!r}'.format(name))
        return k


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


cdef class WallSurface:
    """
    Represents a wall surface in contact with the contents of a reactor.
    """
    def __cinit__(self, Wall wall, int side):
        self.wall = wall
        self.cxxwall = &wall.wall
        self.side = side
        self._kinetics = None

    property kinetics:
        """
        The `InterfaceKinetics` object used for calculating reaction
        rates on this wall surface.
        """
        def __get__(self):
            return self._kinetics
        def __set__(self, Kinetics k):
            self._kinetics = k
            self.wall._set_kinetics()

    property coverages:
        """
        The fraction of sites covered by each surface species.
        """
        def __get__(self):
            if self._kinetics is None:
                raise Exception('No kinetics manager present')
            self.cxxwall.syncCoverages(self.side)
            return self._kinetics.coverages
        def __set__(self, coverages):
            if self._kinetics is None:
                raise Exception("Can't set coverages before assigning kinetics manager.")
            if len(coverages) != self._kinetics.n_species:
                raise ValueError('Incorrect number of site coverages specified')
            cdef np.ndarray[np.double_t, ndim=1] data = \
                    np.ascontiguousarray(coverages, dtype=np.double)
            self.cxxwall.setCoverages(self.side, &data[0])

    def add_sensitivity_reaction(self, int m):
        self.cxxwall.addSensitivityReaction(self.side, m)


cdef class Wall:
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

    A heterogeneous reaction mechanism may be specified for one or both of the
    wall surfaces. The mechanism object (typically an instance of class
    `Interface`) must be constructed so that it is properly linked to
    the object representing the fluid in the reactor the surface in question
    faces. The surface temperature on each side is taken to be equal to the
    temperature of the reactor it faces.
    """

    def __init__(self, left, right, *, name=None, A=None, K=None, U=None,
                 Q=None, velocity=None, kinetics=(None,None)):
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
        :param kinetics:
            Surface reaction mechanisms for the left-facing and right-facing
            surface, respectively. These must be instances of class Kinetics,
            or of a class derived from Kinetics, such as Interface. If
            chemistry occurs on only one side, enter ``None`` for the
            non-reactive side.
        """
        self.left_surface = WallSurface(self, 0)
        self.right_surface = WallSurface(self, 1)

        self._velocity_func = None
        self._heat_flux_func = None

        self._install(left, right)
        if name is not None:
            self.name = name
        else:
            reactor_counts['Wall'] += 1
            n = reactor_counts['Wall']
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
        if kinetics[0] is not None:
            self.left_surface.kinetics = kinetics[0]
        if kinetics[1] is not None:
            self.right_surface.kinetics = kinetics[1]

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

    property left:
        """ The left surface of this wall. """
        def __get__(self):
            return self.left_surface

    property right:
        """ The right surface of this wall. """
        def __get__(self):
            return self.right_surface

    property expansion_rate_coeff:
        """
        The coefficient *K* [m/s/Pa] that determines the velocity of the wall
        as a function of the pressure difference between the adjacent reactors.
        """
        def __get__(self):
            return self.wall.getExpansionRateCoeff()
        def __set__(self, double val):
            self.wall.setExpansionRateCoeff(val)

    property area:
        """ The wall area [m^2]. """
        def __get__(self):
            return self.wall.area()
        def __set__(self, double value):
            self.wall.setArea(value)

    property heat_transfer_coeff:
        """the overall heat transfer coefficient [W/m^2/K]"""
        def __get__(self):
            return self.wall.getHeatTransferCoeff()
        def __set__(self, double value):
            self.wall.setHeatTransferCoeff(value)

    property emissivity:
        """The emissivity (nondimensional)"""
        def __get__(self):
            return self.wall.getEmissivity()
        def __set__(self, double value):
            self.wall.setEmissivity(value)

    def set_velocity(self, v):
        """
        The wall velocity [m/s]. May be either a constant or an arbirary
        function of time. See `Func1`.
        """
        cdef Func1 f
        if isinstance(v, Func1):
            f = v
        else:
            f = Func1(v)

        self._velocity_func = f
        self.wall.setVelocity(f.func)

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
        self.wall.setHeatFlux(f.func)

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

    def _set_kinetics(self):
        cdef CxxKinetics* L = (self.left_surface._kinetics.kinetics
                               if self.left_surface._kinetics else NULL)
        cdef CxxKinetics* R = (self.right_surface._kinetics.kinetics
                               if self.right_surface._kinetics else NULL)
        self.wall.setKinetics(L, R)


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
    def __cinit__(self, *args, **kwargs):
        # Children of this abstract class are responsible for allocating dev
        self.dev = NULL

    def __init__(self, upstream, downstream, *, name=None):
        assert self.dev != NULL
        self._rate_func = None

        if name is not None:
            self.name = name
        else:
            reactor_counts[self.__class__.__name__] += 1
            n = reactor_counts[self.__class__.__name__]
            self.name = '{0}_{1}'.format(self.__class__.__name__, n)

        self._install(upstream, downstream)

    def __dealloc__(self):
        del self.dev

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


cdef class MassFlowController(FlowDevice):
    r"""
    A mass flow controller maintains a specified mass
    flow rate independent of upstream and downstream conditions. The equation
    used to compute the mass flow rate is

    .. math::

        \dot m = \max(\dot m_0, 0.0),

    where :math:`\dot m_0` is either a constant value or a function of time.
    Note that if :math:`\dot m_0 < 0`, the mass flow rate will be set to zero,
    since reversal of the flow direction is not allowed.

    Unlike a real mass flow controller, a MassFlowController object will
    maintain the flow even if the downstream pressure is greater than the
    upstream pressure.  This allows simple implementation of loops, in which
    exhaust gas from a reactor is fed back into it through an inlet. But note
    that this capability should be used with caution, since no account is
    taken of the work required to do this.
    """
    def __cinit__(self, *args, **kwargs):
        self.dev = new CxxMassFlowController()

    def __init__(self, upstream, downstream, *, name=None, mdot=None):
        super().__init__(upstream, downstream, name=name)
        if mdot is not None:
            self.set_mass_flow_rate(mdot)

    def set_mass_flow_rate(self, m):
        """
        Set the mass flow rate [kg/s] through this controller to be either
        a constant or an arbitrary function of time. See `Func1`.

        >>> mfc.set_mass_flow_rate(0.3)
        >>> mfc.set_mass_flow_rate(lambda t: 2.5 * exp(-10 * (t - 0.5)**2))
        """
        cdef Func1 f
        if isinstance(m, Func1):
            f = m
        else:
            f = Func1(m)

        self._rate_func = f
        self.dev.setFunction(f.func)


cdef class Valve(FlowDevice):
    r"""
    In Cantera, a `Valve` is a flow devices with mass flow rate that is a
    function of the pressure drop across it. The default behavior is linear:

    .. math:: \dot m = K_v (P_1 - P_2)

    if :math:`P_1 > P_2.` Otherwise, :math:`\dot m = 0`.
    However, an arbitrary function can also be specified, such that

    .. math:: \dot m = F(P_1 - P_2)

    if :math:`P_1 > P_2`, or :math:`\dot m = 0` otherwise.
    It is never possible for the flow to reverse and go from the downstream
    to the upstream reactor/reservoir through a line containing a Valve object.

    :class:`Valve` objects are often used between an upstream reactor and a
    downstream reactor or reservoir to maintain them both at nearly the same
    pressure. By setting the constant :math:`K_v` to a sufficiently large
    value, very small pressure differences will result in flow between the
    reactors that counteracts the pressure difference.
    """
    def __cinit__(self, *args, **kwargs):
        self.dev = new CxxValve()

    def __init__(self, upstream, downstream, *, name=None, K=None):
        super().__init__(upstream, downstream, name=name)
        if K is not None:
            self.set_valve_coeff(K)

    def set_valve_coeff(self, k):
        """
        Set the relationship betwen mass flow rate and the pressure drop across
        the valve. If a number is given, it is the proportionality constant
        [kg/s/Pa]. If a function is given, it should compute the mass flow
        rate [kg/s] given the pressure drop [Pa].

        >>> V = Valve(res1, reactor1)
        >>> V.set_valve_coeff(1e-4)
        >>> V.set_valve_coeff(lambda dP: (1e-5 * dP)**2)
        """
        cdef double kv
        cdef Func1 f
        if isinstance(k, numbers.Real):
            kv = k
            self.dev.setParameters(1, &kv)
            return

        if isinstance(k, Func1):
            f = k
        else:
            f = Func1(k)
        self._rate_func = f
        self.dev.setFunction(f.func)


cdef class PressureController(FlowDevice):
    r"""
    A PressureController is designed to be used in conjunction with another
    'master' flow controller, typically a `MassFlowController`. The master
    flow controller is installed on the inlet of the reactor, and the
    corresponding `PressureController` is installed on on outlet of the
    reactor. The `PressureController` mass flow rate is equal to the master
    mass flow rate, plus a small correction dependent on the pressure
    difference:

    .. math:: \dot m = \dot m_{\rm master} + K_v(P_1 - P_2).
    """
    def __cinit__(self, *args, **kwargs):
        self.dev = new CxxPressureController()

    def __init__(self, upstream, downstream, *, name=None, master=None, K=None):
        super().__init__(upstream, downstream, name=name)
        if master is not None:
            self.set_master(master)
        if K is not None:
            self.set_pressure_coeff(K)

    def set_pressure_coeff(self, double k):
        """
        Set the proportionality constant *k* [kg/s/Pa] between the pressure
        drop and the mass flow rate.
        """
        self.dev.setParameters(1, &k)

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

    def add_reactor(self, ReactorBase r):
        """Add a reactor to the network."""
        self._reactors.append(r)
        self.net.addReactor(r.rbase)

    def advance(self, double t):
        """
        Advance the state of the reactor network in time from the current
        time to time *t* [s], taking as many integrator timesteps as necessary.
        """
        self.net.advance(t)

    def step(self, double t):
        """
        Take a single internal time step toward time *t* [s]. The time after
        taking the step is returned.
        """
        return self.net.step(t)

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

    def set_max_time_step(self, double t):
        """
        Set the maximum time step *t* [s] that the integrator is allowed
        to use.
        """
        self.net.setMaxTimeStep(t)

    property max_err_test_fails:
        """
        The maximum number of error test failures permitted by the CVODES
        integrator in a single time step.
        """
        def __set__(self, n):
            self.net.setMaxErrTestFails(n)

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
            return pybool(self.verbose())
        def __set__(self, pybool v):
            self.net.setVerbose(v)

    def sensitivity(self, species, int p, int r=0):
        if isinstance(species, int):
            return self.net.sensitivity(species,p)
        elif isinstance(species, (str, unicode)):
            return self.net.sensitivity(stringify(species), p, r)

    def sensitivities(self):
        cdef np.ndarray[np.double_t, ndim=2] data = \
                np.empty((self.n_vars, self.n_sensitivity_params))
        cdef int p,k
        for p in range(self.n_sensitivity_params):
            for k in range(self.n_vars):
                data[k,p] = self.net.sensitivity(k,p)
        return data

    def sensitivity_parameter_name(self, int p):
        return pystr(self.net.sensitivityParameterName(p))

    property n_sensitivity_params:
        def __get__(self):
            return self.net.nparams()

    property n_vars:
        def __get__(self):
            return self.net.neq()
