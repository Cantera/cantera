from collections import defaultdict
import numbers

reactor_counts = defaultdict(int)

cdef class ReactorBase:
    reactorType = "None"
    def __cinit__(self, *args, **kwargs):
        self.rbase = newReactor(stringify(self.reactorType))

    def __init__(self, ThermoPhase contents=None, name=None, **kwargs):
        self._inlets = []
        self._outlets = []
        self._walls = []
        if isinstance(contents, ThermoPhase):
            self.insert(contents)

        if name is not None:
            self.name = name
        else:
            reactor_counts[self.reactorType] += 1
            n = reactor_counts[self.reactorType]
            self.name = '{0}_{1}'.format(self.reactorType, n)

    def __dealloc__(self):
        del self.rbase

    def insert(self, _SolutionBase solution):
        self._thermo = solution
        self.rbase.setThermoMgr(deref(solution.thermo))

    property name:
        def __get__(self):
            return pystr(self.rbase.name())

        def __set__(self, name):
            self.rbase.setName(stringify(name))

    property thermo:
        def __get__(self):
            self.rbase.restoreState()
            return self._thermo

    property volume:
        def __get__(self):
            return self.rbase.volume()

        def __set__(self, double value):
            self.rbase.setInitialVolume(value)

    property T:
        def __get__(self):
            return self.thermo.T

    property density:
        def __get__(self):
            return self.thermo.density

    property Y:
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

    def _addInlet(self, inlet):
        """
        Store a reference to *inlet* to prevent it from being prematurely
        garbage collected.
        """
        self._inlets.append(inlet)

    def _addOutlet(self, outlet):
        """
        Store a reference to *outlet* to prevent it from being prematurely
        garbage collected.
        """
        self._outlets.append(outlet)

    def _addWall(self, wall):
        """
        Store a reference to *wall* to prevent it from being prematurely
        garbage collected.
        """
        self._walls.append(wall)


cdef class Reactor(ReactorBase):
    reactorType = "Reactor"

    def __cinit__(self, *args, **kwargs):
        self.reactor = <CxxReactor*>(self.rbase)

    def __init__(self, contents=None, *, name=None, energy='on', **kwargs):
        super().__init__(contents, **kwargs)

        if energy == 'off':
            self.energyEnabled = False
        elif energy != 'on':
            raise ValueError("'energy' must be either 'on' or 'off'")

    def insert(self, _SolutionBase solution):
        ReactorBase.insert(self, solution)
        self._kinetics = solution
        self.reactor.setKineticsMgr(deref(solution.kinetics))

    property kinetics:
        def __get__(self):
            self.rbase.restoreState()
            return self._kinetics

    property energyEnabled:
        def __get__(self):
            return self.reactor.energyEnabled()

        def __set__(self, pybool value):
            self.reactor.setEnergy(int(value))


cdef class Reservoir(ReactorBase):
    reactorType = "Reservoir"


cdef class ConstPressureReactor(Reactor):
    reactorType = "ConstPressureReactor"


cdef class FlowReactor(Reactor):
    reactorType = "FlowReactor"

    property massFlowRate:
        def __set__(self, double value):
            (<CxxFlowReactor*>self.reactor).setMassFlowRate(value)

    property speed:
        def __get__(self):
            return (<CxxFlowReactor*>self.reactor).speed()

    property distance:
        def __get__(self):
            return (<CxxFlowReactor*>self.reactor).distance()


cdef class Wall:
    def __cinit__(self, *args, **kwargs):
        self.wall = new CxxWall()

    def __init__(self, left, right, *, name=None, A=None, K=None, U=None,
                 Q=None, velocity=None, kinetics=(None,None)):
        self._velocityFunc = None
        self._heatFluxFunc = None
        self._leftKinetics = None
        self._rightKinetics = None

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
            self.expansionRateCoeff = K
        if U is not None:
            self.heatTransferCoeff = U
        if Q is not None:
            self.setHeatFlux(Q)
        if velocity is not None:
            self.setVelocity(velocity)
        if kinetics[0] is not None:
            self.leftKinetics = kinetics[0]
        if kinetics[1] is not None:
            self.leftKinetics = kinetics[1]

    def _install(self, ReactorBase left, ReactorBase right):
        left._addWall(self)
        right._addWall(self)
        self.wall.install(deref(left.rbase), deref(right.rbase))

    property expansionRateCoeff:
        def __get__(self):
            return self.wall.getExpansionRateCoeff()
        def __set__(self, double val):
            self.wall.setExpansionRateCoeff(val)

    property area:
        def __get__(self):
            return self.wall.area()
        def __set__(self, double value):
            self.wall.setArea(value)

    property heatTransferCoeff:
        def __get__(self):
            return self.wall.getHeatTransferCoeff()
        def __set__(self, double value):
            self.wall.setHeatTransferCoeff(value)

    property emissivity:
        def __get__(self):
            return self.wall.getEmissivity()
        def __set__(self, double value):
            self.wall.setEmissivity(value)

    def setVelocity(self, v):
        cdef Func1 f
        if isinstance(v, Func1):
            f = v
        else:
            f = Func1(v)

        self._velocityFunc = f
        self.wall.setVelocity(f.func)

    def setHeatFlux(self, q):
        cdef Func1 f
        if isinstance(q, Func1):
            f = q
        else:
            f = Func1(q)

        self._heatFluxFunc = f
        self.wall.setHeatFlux(f.func)

    def vdot(self, double t):
        return self.wall.vdot(t)

    def qdot(self, double t):
        return self.wall.Q(t)

    property leftKinetics:
        def __get__(self):
            return self._leftKinetics
        def __set__(self, Kinetics k):
            self._leftKinetics = k
            self._setKinetics()

    property rightKinetics:
        def __get__(self):
            return self._rightKinetics
        def __set__(self, Kinetics k):
            self._rightKinetics = k
            self._setKinetics()

    def _setKinetics(self):
        cdef CxxKinetics* L = (self._leftKinetics.kinetics
                               if self._leftKinetics else NULL)
        cdef CxxKinetics* R = (self._rightKinetics.kinetics
                               if self._rightKinetics else NULL)
        self.wall.setKinetics(L, R)

    property leftCoverages:
        def __get__(self):
            if self._leftKinetics is None:
                raise Exception('No kinetics manager present')
            self.wall.syncCoverages(0)
            return self.leftKinetics.coverages
        def __set__(self, coverages):
            if self._leftKinetics is None:
                raise Exception("Can't set coverages before assigning kinetics manager.")
            if len(coverages) != self._leftKinetics.nSpecies:
                raise ValueError('Incorrect number of site coverages specified')
            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(coverages, dtype=np.double)
            self.wall.setCoverages(0, &data[0])

    property rightCoverages:
        def __get__(self):
            if self._rightKinetics is None:
                raise Exception('No kinetics manager present')
            self.wall.syncCoverages(1)
            return self._rightKinetics.coverages
        def __set__(self, coverages):
            if self._rightKinetics is None:
                raise Exception("Can't set coverages before assigning kinetics manager.")
            if len(coverages) != self._rightKinetics.nSpecies:
                raise ValueError('Incorrect number of site coverages specified')
            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(coverages, dtype=np.double)
            self.wall.setCoverages(1, &data[0])


cdef class FlowDevice:
    def __cinit__(self, *args, **kwargs):
        # Children of this abstract class are responsible for allocating dev
        self.dev = NULL

    def __init__(self, upstream, downstream, *, name=None):
        assert self.dev != NULL
        self._rateFunc = None

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
        upstream._addOutlet(self)
        downstream._addInlet(self)
        self.dev.install(deref(upstream.rbase), deref(downstream.rbase))

    def mdot(self, double t):
        return self.dev.massFlowRate(t)


cdef class MassFlowController(FlowDevice):
    def __cinit__(self, *args, **kwargs):
        self.dev = new CxxMassFlowController()

    def __init__(self, upstream, downstream, *, name=None, mdot=None):
        super().__init__(upstream, downstream, name=name)
        if mdot is not None:
            self.setMassFlowRate(mdot)

    def setMassFlowRate(self, m):
        cdef Func1 f
        if isinstance(m, Func1):
            f = m
        else:
            f = Func1(m)

        self._rateFunc = f
        self.dev.setFunction(f.func)


cdef class Valve(FlowDevice):
    def __cinit__(self, *args, **kwargs):
        self.dev = new CxxValve()

    def __init__(self, upstream, downstream, *, name=None, K=None):
        super().__init__(upstream, downstream, name=name)
        if K is not None:
            self.setValveCoeff(K)

    def setValveCoeff(self, k):
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
        self._rateFunc = f
        self.dev.setFunction(f.func)


cdef class PressureController(FlowDevice):
    def __cinit__(self, *args, **kwargs):
        self.dev = new CxxPressureController()

    def __init__(self, upstream, downstream, *, name=None, master=None, K=None):
        super().__init__(upstream, downstream, name=name)
        if master is not None:
            self.setMaster(master)
        if K is not None:
            self.setPressureCoeff(K)

    def setPressureCoeff(self, double k):
        self.dev.setParameters(1, &k)

    def setMaster(self, FlowDevice d):
        (<CxxPressureController*>self.dev).setMaster(d.dev)


cdef class ReactorNet:
    cdef CxxReactorNet* net
    cdef list _reactors

    def __cinit__(self, *args, **kwargs):
        self.net = new CxxReactorNet()

    def __init__(self, reactors=()):
        self._reactors = []  # prevents premature garbage collection
        for R in reactors:
            self.addReactor(R)

    def addReactor(self, ReactorBase r):
        self._reactors.append(r)
        self.net.addReactor(r.rbase)

    def advance(self, double t):
        self.net.advance(t)

    def step(self, double t):
        return self.net.step(t)

    property time:
        def __get__(self):
            return self.net.time()

    def setInitialTime(self, double t):
        self.net.setInitialTime(t)

    def setMaxTimeStep(self, double t):
        self.net.setMaxTimeStep(t)

    property rtol:
        def __get__(self):
            return self.net.rtol()
        def __set__(self, tol):
            self.net.setTolerances(tol, -1)

    property atol:
        def __get__(self):
            return self.net.atol()
        def __set__(self, tol):
            self.net.setTolerances(-1, tol)

    property verbose:
        def __get__(self):
            return pybool(self.verbose())
        def __set__(self, pybool v):
            self.net.setVerbose(v)
