
from Cantera import OneAtm
from Cantera.exceptions import CanteraError
from Cantera.Flow import Flow1D
from Cantera.boundaries1D import Inlet1D, Surf1D, Symm1D
from Numeric import array, zeros, arrayrange

from Cantera.gases import IdealGasMix, GRI30
from Cantera.solve import solve
#from Cantera.esolve import esolve
from Cantera.OneDim import OneDim
from Cantera.FlowBoundary import Inlet, Outlet, SymmPlane
from Cantera import stoich
import math

class BurnerFlame:
    """One-dimensional flat, premixed flames.

    flame = BurnerFlame(gas, domain, fuel, oxidizer, inert, grid, pressure)

    arguments:

    gas       ---  an object representing the gas mixture
    domain    ---  [zmin, zmax]
    fuel      ---  a string specifying the fuel stream composition
    oxidizer  ---  a string specifying the oxidizer stream composition
    inert     ---  a string specifying the composition of an inert
                   stream (optional)
    grid      ---  a sequence defining the initial grid. The first point
                   should be zmin, and the last one zmax. If omitted,
                   a default grid will be used.
    pressure  --- the pressure, which is treated as constant.
    
    example:
       flame = BurnerFlame(gas = GRI30(),
                           domain = [0.0, 10.0*units.cm],
                           fuel = 'CH4:1',
                           oxidizer = 'O2:1,N2:3.76',
                           grid = [0.0, 0.01, 0.03, 0.06, 0.1],
                           pressure = OneAtm)
    
    """
    
    def __init__(self, gas = None, domain = None,
                 fuel = '', oxidizer = '', inert = '',
                 grid = None, pressure = -1.0):

        # check that all required inputs have been specified
        if not gas or not domain or not fuel or not oxidizer or not pressure:
            raise self.__doc__
        
        self.gas = gas
        self.p = pressure

        dx = (domain[1] - domain[0])

        # if no grid specified, use this one that concentrates points
        # near the burner
        if grid == None:
            grid = dx * array([0.0, 0.01, 0.03, 0.1, 0.3, 0.6, 1.0])

        self.__flow = Flow1D(flow_type = 'OneDim', gas = gas,
                             grid = grid, pressure = self.p)


        #------ these methods are deprecated, but still needed for now.
        self.inlet = Inlet(gas)
        self.outlet = Outlet(gas)
                           
        self.__flow.setBoundaries(left = self.inlet, right = self.outlet)
        #----------------------------------------------------------------

        # The container contains only the Flow1D object. Should be
        # modified at some point to contain an Inlet1D and an Outlet1D
        # object.
        self.__container = OneDim([self.__flow])
        self.start = 0


        # get the compositions of the fuel and oxidizer streams, and
        # calculate the fuel/oxidizer ratio for stoichiometric
        # combustion
        gas.setMoleFractions(fuel)
        self._xfuel = gas.moleFractions()

        gas.setMoleFractions(oxidizer)
        self._xox = gas.moleFractions()

        if inert:
            gas.setMoleFractions(inert)
            self._xinert = gas.moleFractions()
        else:
            self._xinert = zeros(gas.nSpecies(),'d')

        self._stoich_FO = stoich.stoich_fuel_to_oxidizer(gas, fuel, oxidizer)

        
    # TODO: accout for inert stream
    def setEquivRatio(self, phi):
        """Set the equivalence ratio."""
        f_flow = self._stoich_FO * phi
        comp = f_flow * self._xfuel + self._xox
        self.gas.setState_PX(self.p, comp)
        self.inlet.set(X = self.gas.moleFractions())

        
    def setEquilProducts(self):
        """Generate a starting estimate for the flame state.

        The following procedure is used:
        
        1) At the burner, the composition is the specified inlet
           composition;
        2) The last 80% of the domain has constant composition and
           temperature corresponding to the adiabatic equilibrium
           solution;
        3) In the initial 20%, the composition and temperature vary linearly
           from the inlet values to the equilibrium values.

        """
        x0 = self.inlet.X
        self.gas.setState_TPX(self.inlet.T, self.p, x0)
        rho0 = self.gas.density()
        mdot = self.inlet.mdot
        self.gas.equilibrate('HP')
        xp = self.gas.moleFractions()
        xinit = {}
        z0 = 0.2
        teq = self.gas.temperature()
        rhoeq = self.gas.density()
        xinit['T'] = [(0.0, self.inlet.T), (z0, teq), (1.0, teq)]
        xinit['u'] = [(0.0, mdot/rho0), (z0, mdot/rhoeq), (1.0, mdot/rhoeq)]
        for k in range(self.gas.nSpecies()):
            nm = self.gas.speciesName(k)
            x = [(0.0, x0[k]), (z0, xp[k]), (1.0, xp[k])]
            xinit[nm] = x
        self.__flow.setInitialProfiles(xinit)


    def plot(self, plotfile = '', title = '', fmt = 'TECPLOT',
             zone = 'c0', append = 0):
        """Plot the current solution."""
        self.__flow.plotter.plot(fname = plotfile, title = title,
                          fmt = fmt, zone = zone, append=append)        

        
    def setInitialProfiles(self, **init):
        """Specify estimates for the initial profiles.

        """
        self.__flow.setInitialProfiles(init)
        self.start = 1
        
    def restore(self, src = '', solution = ''):
        """Start from a previously-saved solution."""
        self.__container.restore(0, src, solution)
        self.start = 1

    def setTolerances(self, V = None, T = None, Y = None):
        """Set tolerances for convergence for velocity, temperature,
        and mass fractions."""
        self.__flow.setTolerances( V, V, T, Y)

    def prune(self, loglevel = 2):
        """Remove unneeded grid points.

        This method attempts to remove each grid point one by one, and
        calls 'refine' each time to see whether it puts it back. If it does,
        the point is not removed, otherwise it is.
        """
        self.__container.prune(loglevel)

    def refine(self, loglevel = 2):
        """Refine the grid using the current grid refinement parameters."""
        self.__container.refine(loglevel)        

    def show(self):
        """Print a summary of the current solution to the screen."""
        self.__flow.show()

    def stretch(self, factor):
        """Stretch the grid by 'factor'"""
        self.__flow.setGrid(factor*self.__flow.z)
        
    def set(self, **opt):
        """Set options.

        The options that may be set are:

        energy   --- 'on' or 'off'. If 'on', the energy equation is
        solved; otherwise, the temperature is held to the specified
        profile.

        pressure      --- the pressure in Pa.
        
        mdot          --- the inlet mass flow rate per unit area.
        
        equiv_ratio   --- the equivalence ratio

        T_burner      --- the burner surface temperature [K].

        refine        --- a triplet specifying the refinement criteria.
                          See refine.py for more information.

        tol           --- error tolerances for u, V, T, and Y.

        max_jac_age   --- the maximum number of times to use a Jacobian
                          before recomputing it.

        timesteps     --- number and duration of time steps to take
                          when Newton iteration fails. The format is
                          ( number_sequence, initial_stepsize )
                          
        These parameters can be changed as the solution proceeds."""

        # TODO: is this necessary? 
        if self.__container == None:
            self.__container = OneDim([self.__flow,])        
        
        for o in opt.keys():
            v = opt[o]
            if o == 'energy':
                self.__flow.setEnergyEqn(v,loglevel=1)
            elif o == 'pressure':
                self.p = v
                self.__flow.setPressure(v)
            elif o == 'mdot':
                self.inlet.set(mdot = v)
            elif o == 'equiv_ratio':
                self.setEquivRatio(v)
            elif o == 'T_burner':
                self.inlet.set(T = v)
            elif o == 'refine':
                self.__flow.refiner.delta = v
            elif o == 'tol':
                self.__flow.setTolerances(u = v, V = v, T = v, Y = v)
            elif o == 'max_jac_age':
                self.__container.setOptions(max_jac_age = v)
            elif o == 'timesteps':
                self.__container.setOptions(nsteps = v[0], timestep = v[1])


    def solve(self, loglevel = 0):
        """ Solve the flame equations.

        If no starting estimate has been given, setEquilProducts()
        is called to generate one.

        """
        if not self.start:
            self.setEquilProducts()
            self.start = 1
        solve(self.__container, loglevel = loglevel, refine_grid = 1)


##     def esolve(self, loglevel = 0, efactor = 1.0e4):
##         if not self.start:
##             self.setEquilProducts()
##             self.start = 1
##         esolve(self.__container, efactor = efactor, loglevel = loglevel, refine_grid = 1) 
        

    def save(self, soln, desc, file = 'flame.xml'):
        """Save the current solution.

        soln    --- string to identify this solution in the file.
        desc    --- descriptive text string.
        file    --- file name.
        """
        self.__container.save(file, soln, desc)

    def showStatistics(self):
        """Show numerical statistics."""
        self.__container.showStatistics()
    
    



class StagnationFlame:
    """Axisymmetric premixed stagnation-point flames.

    flame = StagnationFlame(gas, domain, fuel, oxidizer, inert, grid, pressure)

    example:
       flame = BurnerFlame(gas = GRI30(),
                           domain = [0.0, 10.0*units.cm],
                           fuel = 'CH4:1',
                           oxidizer = 'O2:1,N2:3.76',
                           grid = [0.0, 0.01, 0.03, 0.06, 0.1],
                           pressure = OneAtm)
    
    """
    
    def __init__(self, gas = None, domain = None,
                 fuel = '', oxidizer = '', inert = '',
                 grid = None, pressure = -1.0):

        if not gas or not domain or not fuel or not oxidizer or not pressure:
            raise self.__doc__
        
        self.gas = gas
        self.p = pressure

        dx = (domain[1] - domain[0])
        self.dx = dx
        
        if grid == None:
            grid = dx * array([0.0, 0.01, 0.03, 0.1, 0.3, 0.6, 1.0])

        self.__flow = Flow1D(flow_type = 'Stag', gas = gas,
                             grid = grid, pressure = self.p)

        self.__left = Inlet1D()

        self.__right = Surf1D()

        self.__container = OneDim([self.__left, self.__flow, self.__right])
        self.start = 0

        # get the compositions of the fuel and oxidizer streams, and
        # calculate the fuel/oxidizer ratio for stoichiometric
        # combustion
        
        gas.setMoleFractions(fuel)
        self._xfuel = gas.moleFractions()

        gas.setMoleFractions(oxidizer)
        self._xox = gas.moleFractions()

        if inert:
            gas.setMoleFractions(inert)
            self._xinert = gas.moleFractions()
        else:
            self._xinert = zeros(gas.nSpecies(),'d')

        self._stoich_FO = stoich.stoich_fuel_to_oxidizer(gas, fuel, oxidizer)


    def nPoints(self):
        return len(self.__flow.z)
    
    def setEquivRatio(self, phi):
        """Set the equivalence ratio."""
        f_flow = self._stoich_FO * phi
        comp = f_flow * self._xfuel + self._xox
        self.gas.setState_PX(self.p, comp)
        self.__left.set(X = self.gas.moleFractions())

        
    def setEquilProducts(self):
        """Set the flame state to chemical equilibrium.

        This is useful to generate a starting estimate.
        """

        x0 = self.__left.X
        self.gas.setState_TPX(self.__left.T, self.p, x0)
        rho0 = self.gas.density()
        
        mdot = self.__left.mdot
        self.gas.equilibrate('HP')
        xp = self.gas.moleFractions()

        xinit = {}
        z0 = 0.2
        teq = self.gas.temperature()
        rhoeq = self.gas.density()

        re = self.dx * mdot / self.gas.viscosity()
        z1 = 1.0 - 1.0/math.sqrt(re)

        tw = self.__right.T
        self.gas.setState_TPX(tw, self.p, x0)
        self.gas.equilibrate('TP')
        x1 = self.gas.moleFractions()
        rho1 = self.gas.density()
        
        xinit['T'] = [(0.0, self.__left.T), (z0, teq), (z1, teq),
                      (1.0, tw)]
        xinit['u'] = [(0.0, mdot/rho0), (1.0, 0.0)]
        xinit['V'] = [(0.0, 0.0), (z1, mdot/(rhoeq*z1*self.dx)), (1.0, 0.0)]
        for k in range(self.gas.nSpecies()):
            nm = self.gas.speciesName(k)
            x = [(0.0, x0[k]), (z0, xp[k]), (z1, xp[k]), (1.0, x1[k])]
            xinit[nm] = x
        self.__flow.setInitialProfiles(xinit)


    def plot(self, plotfile = '', title = '', fmt = 'TECPLOT',
             zone = 'c0', append = 0):
        self.__flow.plotter.plot(fname = plotfile, title = title,
                          fmt = fmt, zone = zone, append=append)        
        
    def setInitialProfiles(self, **init):
        self.__flow.setInitialProfiles(init)
        self.start = 1

    def resid(self):
        return self.__container.resid(1)
        
    def restore(self, src = '', solution = ''):
        self.__container.restore(1,src, solution)
        self.start = 1

    def setTolerances(self, V = None, T = None, Y = None):
        self.__flow.setTolerances( V, V, T, Y)

    def show(self):
        self.__flow.show()

    def stretch(self, factor):
        self.__flow.setGrid(factor*self.__flow.z)

    def enableEnergy(self, pt):
        self.__flow.setEnergyEqn('on',loglevel=1,pt=pt)

    def prune(self, loglevel = 2):
        self.__container.prune(loglevel)

    def refine(self, loglevel = 2):
        self.__container.refine(loglevel)        

    def set(self, **opt):

        if self.__container == None:
            self.__container = OneDim([self.__flow,])        
        
        for o in opt.keys():
            v = opt[o]
            if o == 'energy':
                self.__flow.setEnergyEqn(v,loglevel=1)
            elif o == 'pressure':
                self.p = v
                self.__flow.setPressure(v)
            elif o == 'mdot':
                self.__left.set(mdot = v)
            elif o == 'equiv_ratio':
                self.setEquivRatio(v)
            elif o == 'T_burner':
                self.__left.set(T = v)
            elif o == 'spreadingRate':
                self.__left.set(V = v)                
            elif o == 'T_surface':
                self.__right.set(T = v)                
            elif o == 'refine':
                self.__flow.refiner.delta = v
            elif o == 'efactor':
                self.__flow.setEnergyFactor(v)
            elif o == 'tol':
                self.__flow.setTolerances(u = v, V = v, T = v, Y = v)
            elif o == 'max_jac_age':
                self.__container.setOptions(max_jac_age = v)
            elif o == 'jac_age':
                self.__container.setOptions(max_jac_age = v[0])
                self.__container.setOptions(ts_jac_age = v[1])                                
            elif o == 'timesteps':
                self.__container.setOptions(nsteps = v[0], timestep = v[1])
            else:
                raise CanteraError("unknown option: "+o)

    def solve(self, loglevel = 0):
        if not self.start:
            self.setEquilProducts()
            self.start = 1
        solve(self.__container, loglevel = loglevel, refine_grid = 1)


    def save(self, soln, desc, file = 'flame.xml'):
        self.__container.save(file, soln, desc)

    def showStatistics(self):
        self.__container.showStatistics()
    
    








