cdef class ReactorBase:
    reactorType = "None"
    def __cinit__(self, *args, **kwargs):
        self.rbase = newReactor(stringify(self.reactorType))

    def __init__(self, *args, **kwargs):
        if args and isinstance(args[0], _SolutionBase):
            self.insert(args[0])

    def __dealloc__(self):
        del self.rbase

    def insert(self, _SolutionBase solution):
        self.rbase.setThermoMgr(deref(solution.thermo))

    def volume(self):
        return self.rbase.volume()


cdef class Reactor(ReactorBase):
    reactorType = "Reactor"
    cdef CxxReactor* reactor

    def __cinit__(self, *args, **kwargs):
        self.reactor = <CxxReactor*>(self.rbase)

    def insert(self, _SolutionBase solution):
        ReactorBase.insert(self, solution)
        self.reactor.setKineticsMgr(deref(solution.kinetics))


cdef class Wall:
    def __cinit__(self, *args, **kwargs):
        self.wall = new CxxWall()

    def __init__(self, *args, **kwargs):
        self._expansionRateCoeff = 0.0

    def install(self, ReactorBase left, ReactorBase right):
        self.wall.install(deref(left.rbase), deref(right.rbase))

    property expansionRateCoeff:
        def __get__(self):
            return self._expansionRateCoeff
        def __set__(self, double val):
            self._expansionRateCoeff = val
            self.wall.setExpansionRateCoeff(val)

    property area:
        def __get__(self):
            return self.wall.area()
        def __set__(self, double value):
            self.wall.setArea(value)

cdef class ReactorNet:
    cdef CxxReactorNet* net
    def __cinit__(self, *args, **kwargs):
        self.net = new CxxReactorNet()

    def addReactor(self, ReactorBase r):
        self.net.addReactor(r.rbase)

    def initialize(self, double t=0.0):
        self.net.initialize(t)

    def advance(self, double t):
        self.net.advance(t)

    def step(self, double t):
        return self.net.step(t)
