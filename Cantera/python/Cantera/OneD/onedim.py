from Cantera import *
from Cantera import _cantera
import Numeric

_onoff = {'on':1, 'yes':1, 'off':0, 'no':0}

class Domain1D:
    """One-dimensional domains."""
    
    def __init__(self):
        self._hndl = 0
        
    def __del__(self):
        _cantera.domain_del(self._hndl)

    def domain_hndl(self):
        """Integer used to reference the kernel object."""
        return self._hndl
    
    def type(self):
        """Domain type. Integer."""
        return _cantera.domain_type(self._hndl)
    
    def index(self):
        """Index of this domain in a stack. Returns -1 if this domain
        is not part of a stack."""
        return _cantera.domain_index(self._hndl)
    
    def nComponents(self):
        """Number of solution components at each grid point."""
        return _cantera.domain_nComponents(self._hndl)

    def nPoints(self):
        """Number of grid points belonging to this domain."""
        return _cantera.domain_nPoints(self._hndl)

    def componentName(self, n):
        """Name of the nth component."""
        return _cantera.domain_componentName(self._hndl, n)

    def componentNames(self):
        names = []
        for n in range(self.nComponents()):
            names.append(self.componentName(n))
        return names
    
    def componentIndex(self, name):
        """Index of the component with name 'name'"""
        return _cantera.domain_componentIndex(self._hndl, name)
    
    def setBounds(self, **bounds):        
        """Set the lower and upper bounds on the solution.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (lower_bound, upper_bound)
        tuples as the values.  The keyword 'default' may be used to
        specify default bounds for all unspecified components. The
        keyword 'Y' can be used to stand for all species mass
        fractions in flow domains.  """
        
        d = {}
        if bounds.has_key('default'):
            for n in range(self.nComponents()):
                d[self.componentName(n)] = bounds['default']
            del bounds['default']
        
        for b in bounds.keys():
            if b == 'Y':
                if self.type >= 50:
                    nc = self.nComponents()
                    for n in range(4, nc):
                        d[self.componentName(n)] = bounds[b]
                else:
                    raise CanteraError('Y can only be specified in flow domains.')
            else:
                d[b] = bounds[b]
        for b in d.keys():
            n = self.componentIndex(b)
            _cantera.domain_setBounds(self._hndl, n, d[b][0], d[b][1])

    def bounds(self, component):
        ic = self.componentIndex(component)
        lower = _cantera.domain_lowerBound(self._hndl, ic)
        upper = _cantera.domain_upperBound(self._hndl, ic)
        return (lower, upper)

    def tolerances(self, component):
        ic = self.componentIndex(component)
        r = _cantera.domain_rtol(self._hndl, ic)
        a = _cantera.domain_atol(self._hndl, ic)        
        return (r, a)
    
    def setTolerances(self, **tol):
        """Set the error tolerances. If 'time' is present and
        non-zero, then the values entered will apply to the transient
        problem. Otherwise, they will apply to the steady-state
        problem. 

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (rtol, atol) tuples as the
        values.  The keyword 'default' may be used to specify default
        bounds for all unspecified components. The keyword 'Y' can be
        used to stand for all species mass fractions in flow domains.
        """
        
        d = {}
        if tol.has_key('default'):
            for n in range(self.nComponents()):
                d[self.componentName(n)] = tol['default']
            del tol['default']
                
        itime = 0
        for b in tol.keys():
            if b == 'time': itime = -1
            elif b == 'steady': itime = 1
            elif b == 'Y':
                if self.type >= 50:
                    nc = self.nComponents()
                    for n in range(4, nc):
                        d[self.componentName(n)] = tol[b]
                else:
                    raise CanteraError('Y can only be specified in flow domains.')
            else:
                d[b] = tol[b]
        for b in d.keys():
            n = self.componentIndex(b)
            _cantera.domain_setTolerances(self._hndl, n, d[b][0], d[b][1], itime)

    
    def setupGrid(self, grid):
        """Specify the grid."""
        return _cantera.domain_setupGrid(self._hndl, Numeric.asarray(grid))
    
    def setID(self, id):
        return _cantera.domain_setID(self._hndl, id)
    
    def setDesc(self, desc):
        return _cantera.domain_setDesc(self._hndl, desc)
    
    def grid(self, n = -1):
        if n >= 0:
            return _cantera.domain_grid(self._hndl, n)
        else:
            g = zeros(self.nPoints(),'d')
            for j in range(len(g)):
                g[j] = _cantera.domain_grid(self._hndl, j)
            return g

    def set(self, **options):
        self._set(options)

    def _set(self, options):
        for opt in options.keys():
            v = options[opt]
            if opt == 'grid':
                self.setupGrid(v)
            elif opt == 'name':
                self.setID(v)
            elif opt == 'desc':
                self.setDesc(v)
            elif opt == 'bounds':
                lower, upper = self._dict2arrays(v)
                self.setBounds(lower,upper)
            elif opt == 'tol':
                self.setTolerances(v[0],v[1])
            else:
                raise CanteraError('unknown attribute: '+opt)

    def _dict2arrays(self, d = None, array1 = None, array2 = None):
        nc = self.nComponents()
        if d.has_key('default'):
            a1 = zeros(nc,'d') + d['default'][0]            
            a2 = zeros(nc,'d') + d['default'][1]
            del d['default']
        else:
            if array1: a1 = array(array1)
            else: a1 = zeros(nc,'d')
            if array2: a2 = array(array2)
            else: a2 = zeros(nc,'d')        

        for k in d.keys():
            c = self.componentIndex(k)
            if c >= 0:
                a1[self.componentIndex(k)] = d[k][0]
                a2[self.componentIndex(k)] = d[k][1]
            else:
                raise CanteraError('unknown component '+k)
        return (a1, a2)
    


class Bdry1D(Domain1D):
    
    def __init__(self):
        Domain1D.__init__(self)
    
    def setMdot(self, mdot):
        _cantera.bdry_setMdot(self._hndl, mdot)
        
    def setTemperature(self, t):
        _cantera.bdry_setTemperature(self._hndl, t)

    def setMoleFractions(self, x):
        _cantera.bdry_setMoleFractions(self._hndl, x)

    def temperature(self):
        return _cantera.bdry_temperature(self._hndl)

    def massFraction(self, k):
        return _cantera.bdry_massFraction(self._hndl, k)

    def mdot(self):
        return _cantera.bdry_mdot(self._hndl)

    def set(self, **options):
        for opt in options.keys():
            v = options[opt]
            if opt == 'mdot' or opt == 'massflux':
                self.setMdot(v)
                del options[opt]
            elif opt == 'temperature' or opt == 'T':
                self.setTemperature(v)
                del options[opt]                
            elif opt == 'mole_fractions' or opt == 'X':
                self.setMoleFractions(v)
                del options[opt]                
        self._set(options)


class Inlet(Bdry1D):
    """A one-dimensional inlet.
    Note that an inlet can only be a terminal domain - it must be
    either the leftmost or rightmost domain in a stack.
   """
    def __init__(self, id = 'inlet'):
        Bdry1D.__init__(self)
        self._hndl = _cantera.inlet_new()
        if id: self.setID(id)
        
    def setSpreadRate(self, V0 = 0.0):
        _cantera.inlet_setSpreadRate(self._hndl, V0)
        
        
class Outlet(Bdry1D):
    def __init__(self, id = 'outlet'):
        Bdry1D.__init__(self)        
        self._hndl = _cantera.outlet_new()
        if id: self.setID(id)        
    
class SymmPlane(Bdry1D):
    def __init__(self, id = 'symmetry_plane'):
        Bdry1D.__init__(self)        
        self._hndl = _cantera.symm_new()
        if id: self.setID(id)        

class Surface(Bdry1D):
    def __init__(self, id = 'surface', surface_mech = None):
        Bdry1D.__init__(self)        
        if surface_mech:
            self._hndl = _cantera.reactingsurf_new()
            self.setKineticsMgr(surface_mech)
        else:
            self._hndl = _cantera.surf_new()
        if id: self.setID(id)


    def setKineticsMgr(self, kin):
        _cantera.reactingsurf_setkineticsmgr(self._hndl,
                                             kin.kinetics_hndl())
        
    def setCoverageEqs(self, onoff='on'):
        if onoff == 'on':
            _cantera.reactingsurf_enableCoverageEqs(self._hndl, 1)
        else:
            _cantera.reactingsurf_enableCoverageEqs(self._hndl, 0)

        
class AxisymmetricFlow(Domain1D):
    """An axisymmetric flow domain.

    In an axisymmetric flow domain, the equations solved are the
    similarity equations for the flow in a finite-height gap of
    infinite radial extent. The solution variables are
      u       -- axial velocity 
      V       -- radial velocity divided by radius
      T       -- temperature
      lambda  -- (1/r)(dP/dr)
      Y_k     -- species mass fractions
      
    It may be shown that if the boundary conditions on these variables
    are independent of radius, then a similarity solution to the exact
    governing equations exists in which these variables are all
    independent of radius. This solution holds only in in
    low-Mach-number limit, in which case (dP/dz) = 0, and lambda is a
    constant. (Lambda is treated as a spatially-varying solution
    variable for numerical reasons, but in the final solution it is
    always independent of z.) As implemented here, the governing
    equations assume an ideal gas mixture.  Arbitrary chemistry is
    allowed, as well as arbitrary variation of the transport
    properties.
    """
    def __init__(self, id = 'axisymmetric_flow', gas = None):
        Domain1D.__init__(self)
        iph = gas.thermo_hndl()
        ikin = gas.kinetics_hndl()
        itr = gas.transport_hndl()
        self._hndl = _cantera.stflow_new(iph, ikin, itr)
        if id: self.setID(id)
        self._p = -1.0
        self.setPressure(gas.pressure())
        self.solveEnergyEqn()
        
    def setPressure(self, p):
        """Set the pressure [Pa]. The pressure is a constant, since
        the governing equations are those for the low-Mach-number limit."""
        _cantera.stflow_setPressure(self._hndl, p)
        self._p = p

    def setTransportModel(self, transp):
        itr = transp.transport_hndl()
        _cantera.stflow_setTransport(self._hndl, itr)

    def pressure(self):
        return self._p
        
    def setFixedTempProfile(self, temp):
        """Set the fixed temperature profile.  This profile is used
        whenever the energy equation is disabled.  """
        return _cantera.stflow_setFixedTempProfile(self._hndl, temp)
    
    def solveSpeciesEqs(self, flag = 1):
        """Enable or disable solving the species equations. If invoked
        with no arguments or with a non-zero argument, the species
        equations will be solved. If invoked with a zero argument,
        they will not be, and instead the species profiles will be
        held at their initial values. Default: species equations
        enabled."""
        return _cantera.stflow_solveSpeciesEqs(self._hndl, flag)
    
    def solveEnergyEqn(self, flag = 1):
        """Enable or disable solving the energy equation. If invoked
        with no arguments or with a non-zero argument, the energy
        equations will be solved. If invoked with a zero argument,
        it will not be, and instead the temperature profiles will be
        held to the one specified by the call to setFixedTempProfile.
        Default: energy equation enabled."""        
        return _cantera.stflow_solveEnergyEqn(self._hndl, flag)

    def set(self, **opt):
        for o in opt.keys():
            v = opt[o]
            if o == 'P' or o == 'pressure':
                self.setPressure(v)
                del opt[o]
            elif o == 'energy':
                self.solveEnergyEqn(flag = _onoff[v])
            else:
                self._set(opt)
                

class Stack:
    
    """ Class Stack is a container for one-dimensional domains. It
    also holds the multi-domain solution vector, and controls the
    process of finding the solution.  

    Domains are ordered left-to-right, with domain number 0 at the left.
    """
    def __init__(self, domains = None):
        nd = len(domains)
        hndls = Numeric.zeros(nd,'i')
        for n in range(nd):
            hndls[n] = domains[n].domain_hndl()
        self._hndl = _cantera.sim1D_new(hndls)
        self._domains = domains
        
    def __del__(self):
        _cantera.sim1D_del(self._hndl)

    def setValue(self, dom, comp, localPoint, value):
        idom = dom.domain_hndl()
        _cantera.sim1D_setValue(self._hndl, idom,
                                comp, localPoint, value)
    
    def setProfile(self, dom, comp, pos, v):
        idom = dom.index()
        icomp = dom.componentIndex(comp)
        _cantera.sim1D_setProfile(self._hndl, idom, icomp,
                                  Numeric.asarray(pos), Numeric.asarray(v))
        
    def setFlatProfile(self, dom, comp, v):
        idom = dom.index()
        icomp = dom.componentIndex(comp)        
        _cantera.sim1D_setFlatProfile(self._hndl, idom, icomp, v)
    
    def showSolution(self, fname='-'):
        _cantera.sim1D_showSolution(self._hndl, fname)
        
    def setTimeStep(self, stepsize, nsteps):
        _cantera.sim1D_setTimeStep(self._hndl, stepsize,
                                   Numeric.asarray(nsteps))
        
    def getInitialSoln(self):
        _cantera.sim1D_getInitialSoln(self._hndl)
            
    def solve(self, loglevel=1, refine_grid=1):
        return _cantera.sim1D_solve(self._hndl, loglevel, refine_grid)
    
    def refine(self, loglevel=1):
        return _cantera.sim1D_refine(self._hndl, loglevel)
    
    def setRefineCriteria(self, domain = None, ratio = 10.0, slope = 0.8,
                          curve = 0.8, prune = 0.05):
        idom = domain.index()
        return _cantera.sim1D_setRefineCriteria(self._hndl,
                                                idom, ratio, slope, curve, prune)
    def save(self, file = 'soln.xml', id = 'solution', desc = 'none'):
        return _cantera.sim1D_save(self._hndl, file, id, desc)
    
    def restore(self, file = 'soln.xml', id = 'solution'):
        return _cantera.sim1D_restore(self._hndl, file, id)

    def showStats(self):
        return _cantera.sim1D_writeStats(self._hndl)

    def domainIndex(self, name):
        return _cantera.sim1D_domainIndex(self._hndl, name)

    def value(self, domain, component, localPoint):
        icomp = domain.componentIndex(component)
        idom = domain.index()
        return _cantera.sim1D_value(self._hndl, idom, icomp, localPoint)

    def profile(self, domain, component):
        np = domain.nPoints()
        x = zeros(np,'d')
        for n in range(np):
            x[n] = self.value(domain, component, n)
        return x

    def workValue(self, dom, icomp, localPoint):
        idom = dom.index()        
        return _cantera.sim1D_workValue(self._hndl, idom, icomp, localPoint)
    
    def eval(self, rdt, count=1):
        return _cantera.sim1D_eval(self._hndl, rdt, count)
    
    def setMaxJacAge(self, ss_age, ts_age):
        return _cantera.sim1D_setMaxJacAge(self._hndl, ss_age, ts_age)
    
    def timeStepFactor(self, tfactor):
        return _cantera.sim1D_timeStepFactor(self._hndl, tfactor)
    
    def setTimeStepLimits(self, tsmin, tsmax):
        return _cantera.sim1D_setTimeStepLimits(self._hndl, tsmin, tsmax)

def clearDomains():
    _cantera.domain_clear()

def clearSim1D():
    _cantera.sim1D_clear()
