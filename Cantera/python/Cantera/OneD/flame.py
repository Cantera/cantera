from onedim import *
import Numeric

class BurnerFlame(Stack):
    """A burner-stabilized flat flame."""
    
    def __init__(self, gas = None, burner = None, outlet = None, grid = None):
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
        self.pressure = gas.pressure()
        self.flame = AxisymmetricFlow('flame',gas = gas)
        self.flame.setupGrid(grid)
        Stack.__init__(self, [self.burner, self.flame, self.outlet])
        self.setRefineCriteria()
        self._initialized = 0


    def init(self):
        """Set the initial guess for the solution."""
        gas = self.gas
        nsp = gas.nSpecies()
        yin = Numeric.zeros(nsp, 'd')
        for k in range(nsp):
            yin[k] = self.burner.massFraction(k)
        gas.setState_TPY(self.burner.temperature(), self.pressure, yin)
        u0 = self.burner.mdot()/gas.density()
        t0 = self.burner.temperature()
        
        # get adiabatic flame temperature and composition
        gas.equilibrate('HP')
        teq = gas.temperature()
        yeq = gas.massFractions()
        u1 = self.burner.mdot()/gas.density()

        z1 = 0.2
        locs = Numeric.array([0.0, z1, 1.0],'d')
        self.setProfile('u', locs, [u0, u1, u1])
        self.setProfile('T', locs, [t0, teq, teq])
        for n in range(nsp):
            self.setProfile(gas.speciesName(n), locs, [yin[n], yeq[n], yeq[n]])
        self._initialized = 1


    def solve(self, loglevel = 1, refine_grid = 1):
        if not self._initialized: self.init()
        Stack.solve(self, loglevel = loglevel, refine_grid = refine_grid)


    def setRefineCriteria(self, ratio = 10.0, slope = 0.8, curve = 0.8, prune = 0.0):
        Stack.setRefineCriteria(self, domain = self.flame, ratio = ratio, slope = slope, curve = curve,
                                prune = prune)

    def setProfile(self, component, locs, vals):
        self._initialized = 1
        Stack.setProfile(self, self.flame, component, locs, vals)

    def set(self, tol = None, energy = '', tol_time = None):
        if tol:
            self.flame.setTolerances(default = tol)
        if tol_time:
            self.flame.setTolerances(default = tol_time, time = 1)
        if energy:
            self.flame.set(energy = energy)
        
    def T(self, point = -1):
        return self.solution('T', point)

    def u(self, point = -1):
        return self.solution('u', point)        
        
    def V(self, point = -1):
        return self.solution('V', point)                

    def solution(self, component = '', point = -1):
        if point >= 0: return self.value(self.flame, component, point)
        else: return self.profile(self.flame, component)        

    def setGasState(self, j):
        nsp = self.gas.nSpecies()
        y = Numeric.zeros(nsp, 'd')
        for n in range(nsp):
            nm = self.gas.speciesName(n)
            y[n] = self.solution(nm, j)
        self.gas.setState_TPY(self.T(j), self.flame.pressure(), y)
        
                

        
        
