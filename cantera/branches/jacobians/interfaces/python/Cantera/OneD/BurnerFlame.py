from onedim import *
from Cantera.num import array, zeros

class BurnerFlame(Stack):
    """A burner-stabilized flat flame."""

    def __init__(self, gas = None, burner = None, outlet = None, grid = None):
        """
        :param gas:
            object to use to evaluate all gas properties and reaction
            rates. Required
        :param burner:
            Inlet object representing the burner. Optional;
            if not supplied, one will be created with name ``burner``
        :param outlet:
            Outlet object representing the outlet. Optional;
            if not supplied, one will be created with name ``outlet``
        :param grid:
            array of initial grid points

        A domain of type :class:`.AxisymmetricFlow` named ``flame`` will be
        created to represent the flame. The three domains comprising the stack
        are stored as ``self.burner``, ``self.flame``, and ``self.outlet``.
        """

        if burner:
            self.burner = burner
        else:
            self.burner = Inlet('burner')
        self.gas = gas
        self.burner.set(temperature = gas.temperature())
        if outlet:
            self.outlet = outlet
        else:
            self.outlet = Outlet('outlet')
        self.flame = AxisymmetricFlow('flame',gas = gas)
        self.flame.setupGrid(grid)
        Stack.__init__(self, [self.burner, self.flame, self.outlet])
        self.setRefineCriteria()


    def init(self):
        """Set the initial guess for the solution. The adiabatic flame
        temperature and equilibrium composition are computed for the
        burner gas composition. The temperature profile rises linearly
        in the first 20% of the flame to Tad, then is flat. The mass
        fraction profiles are set similarly.
        """
        self.getInitialSoln()
        gas = self.gas
        nsp = gas.nSpecies()
        yin = zeros(nsp, 'd')
        for k in range(nsp):
            yin[k] = self.burner.massFraction(k)
        gas.setState_TPY(self.burner.temperature(), self.flame.pressure(), yin)
        u0 = self.burner.mdot()/gas.density()
        t0 = self.burner.temperature()

        # get adiabatic flame temperature and composition
        gas.equilibrate('HP',solver=1)
        teq = gas.temperature()
        yeq = gas.massFractions()
        u1 = self.burner.mdot()/gas.density()

        z1 = 0.2
        locs = array([0.0, z1, 1.0],'d')
        self.setProfile('u', locs, [u0, u1, u1])
        self.setProfile('T', locs, [t0, teq, teq])
        for n in range(nsp):
            self.setProfile(gas.speciesName(n), locs, [yin[n], yeq[n], yeq[n]])
        self._initialized = 1


    def solve(self, loglevel = 1, refine_grid = 1):
        if not self._initialized: self.init()
        Stack.solve(self, loglevel = loglevel, refine_grid = refine_grid)


    def setRefineCriteria(self, ratio = 10.0, slope = 0.8,
                          curve = 0.8, prune = 0.0):
        Stack.setRefineCriteria(self, domain = self.flame,
                                ratio = ratio, slope = slope, curve = curve,
                                prune = prune)

    def setGridMin(self, gridmin):
        Stack.setGridMin(self, self.flame, gridmin)

    def setProfile(self, component, locs, vals):
        self._initialized = 1
        Stack.setProfile(self, self.flame, component, locs, vals)

    def set(self, tol = None, energy = '', tol_time = None):
        """Set parameters.

        :param tol:
            (rtol, atol) for steady-state
        :param tol_time:
            (rtol, atol) for time stepping
        :param energy:
            ``'on'`` or ``'off'`` to enable or disable the energy equation
        """
        if tol:
            self.flame.setTolerances(default = tol)
        if tol_time:
            self.flame.setTolerances(default = tol_time, time = 1)
        if energy:
            self.flame.set(energy = energy)

    def T(self, point = -1):
        """Temperature profile or value at one point."""
        return self.solution('T', point)

    def u(self, point = -1):
        """Axial velocity profile or value at one point."""
        return self.solution('u', point)

    def V(self, point = -1):
        """Radial velocity profile or value at one point."""
        return self.solution('V', point)

    def solution(self, component = '', point = -1):
        """Solution component at one point, or full profile if no
        point specified."""
        if point >= 0: return self.value(self.flame, component, point)
        else: return self.profile(self.flame, component)

    def setGasState(self, j):
        """Set the state of the object representing the gas to the
        current solution at grid point *j*."""
        nsp = self.gas.nSpecies()
        y = zeros(nsp, 'd')
        for n in range(nsp):
            nm = self.gas.speciesName(n)
            y[n] = self.solution(nm, j)
        self.gas.setState_TPY(self.T(j), self.flame.pressure(), y)

fix_docs(BurnerFlame)
