from onedim import *
from Cantera.num import array, zeros

class StagnationFlow(Stack):
    """An axisymmetric flow impinging on a surface at normal incidence."""

    def __init__(self, gas = None, surfchem = None, grid = None):
        """
        :param gas:
            object to use to evaluate all gas properties and reaction
            rates. Required.
        :param surfchem:
            object used to evaluate surface reaction rates. If omitted,
            surface will be treated as inert.
        :param grid:
            array of initial grid points

        A domain of type :class:`.AxisymmetricFlow` named ``flow`` will be
        created to represent the flow, and one of type :class:`.Surface` named
        ``surface`` will be created to represent the surface. The three domains
        comprising the stack are stored as ``self.inlet``, ``self.flow``,
        and ``self.surface``.
        """
        self.inlet = Inlet('inlet')
        self.gas = gas
        self.surfchem = surfchem
        self.inlet.set(temperature = gas.temperature())
        self.surface = Surface(id = 'surface', surface_mech = surfchem)
        self.flow = AxisymmetricFlow('flow',gas = gas)
        self.flow.setupGrid(grid)
        Stack.__init__(self, [self.inlet, self.flow, self.surface])
        self.setRefineCriteria()


    def init(self, products = 'inlet'):
        """Set the initial guess for the solution. If products = 'equil',
        then the equilibrium composition at the adiabatic flame temperature
        will be used to form the initial guess. Otherwise the inlet composition
        will be used."""
        self.getInitialSoln()
        gas = self.gas
        nsp = gas.nSpecies()
        yin = zeros(nsp, 'd')
        for k in range(nsp):
            yin[k] = self.inlet.massFraction(k)
        gas.setState_TPY(self.inlet.temperature(), self.flow.pressure(), yin)
        u0 = self.inlet.mdot()/gas.density()
        t0 = self.inlet.temperature()
        V0 = 0.0

        tsurf = self.surface.temperature()

        zz = self.flow.grid()
        dz = zz[-1] - zz[0]

        if products == 'equil':
            gas.equilibrate('HP')
            teq = gas.temperature()
            yeq = gas.massFractions()
            locs = array([0.0, 0.3, 0.7, 1.0],'d')
            self.setProfile('T', locs, [t0, teq, teq, tsurf])
            for n in range(nsp):
                self.setProfile(gas.speciesName(n), locs, [yin[n], yeq[n], yeq[n], yeq[n]])
        else:
            locs = array([0.0, 1.0],'d')
            self.setProfile('T', locs, [t0, tsurf])
            for n in range(nsp):
                self.setProfile(gas.speciesName(n), locs, [yin[n], yin[n]])

        locs = array([0.0, 1.0],'d')
        self.setProfile('u', locs, [u0, 0.0])
        self.setProfile('V', locs, [V0, V0])

        self._initialized = 1


    def solve(self, loglevel = 1, refine_grid = 1):
        if not self._initialized: self.init()
        Stack.solve(self, loglevel = loglevel, refine_grid = refine_grid)


    def setRefineCriteria(self, ratio = 10.0, slope = 0.8,
                          curve = 0.8, prune = 0.0):
        Stack.setRefineCriteria(self, domain = self.flow,
                                ratio = ratio, slope = slope, curve = curve,
                                prune = prune)

    def setGridMin(self, gridmin):
        Stack.setGridMin(self, self.flow, gridmin)

    def setProfile(self, component, locs, vals):
        self._initialized = 1
        Stack.setProfile(self, self.flow, component, locs, vals)

    def set(self, tol = None, energy = '', tol_time = None):
        """Set parameters.

        :param tol:
            (rtol, atol) for steady-state
        :param tol_time:
            (rtol, atol) for time stepping
        :param energy:
            'on' or 'off' to enable or disable the energy equation
        """
        if tol:
            self.flow.setTolerances(default = tol)
        if tol_time:
            self.flow.setTolerances(default = tol_time, time = 1)
        if energy:
            self.flow.set(energy = energy)

    def T(self, point = -1):
        """The temperature [K]"""
        return self.solution('T', point)

    def u(self, point = -1):
        """The axial velocity [m/s]"""
        return self.solution('u', point)

    def V(self, point = -1):
        """The radial velocity divided by radius [s^-1]"""
        return self.solution('V', point)

    def solution(self, component = '', point = -1):
        """The solution for one specified component. If a point number
        is given, return the value of component *component* at this
        point. Otherwise, return the entire profile for this
        component."""
        if point >= 0: return self.value(self.flow, component, point)
        else: return self.profile(self.flow, component)

    def coverages(self):
        """The coverages of the surface species."""
        nsurf = self.surfchem.nSpecies()
        cov = zeros(nsurf,'d')
        for n in range(nsurf):
            nm = self.surfchem.speciesName(n)
            cov[n] = self.value(self.surface, nm, 0)
        return cov

    def setGasState(self, j):
        """Set the state of the object representing the gas to the
        current solution at grid point *j*."""
        nsp = self.gas.nSpecies()
        y = zeros(nsp, 'd')
        for n in range(nsp):
            nm = self.gas.speciesName(n)
            y[n] = self.solution(nm, j)
        self.gas.setState_TPY(self.T(j), self.flow.pressure(), y)

fix_docs(StagnationFlow)
