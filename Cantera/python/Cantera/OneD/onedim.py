from Cantera import *
from Cantera import _cantera
import Numeric

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
        """Domain type."""
        return _cantera.domain_type(self._hndl)
    
    def index(self):
        """Index of this domain in a stack."""
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
    
    def componentIndex(self, name):
        """Index of the component named 'name'"""
        return _cantera.domain_componentIndex(self._hndl, name)
    
    def setBounds(self, lower, upper):
        """Set the lower and upper bounds on the solution."""
        return _cantera.domain_setBounds(self._hndl,
                                         Numeric.asarray(lower),
                                         Numeric.asarray(upper))
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
    
    def setTolerances(self, rtol, atol, time=0):
        """Set the error tolerances. If 'time' is present and
        non-zero, then the values entered will apply to the transient
        problem. Otherwise, they will apply to the steady-state
        problem.  """
        return _cantera.domain_setTolerances(self._hndl,
                                             Numeric.asarray(rtol),
                                             Numeric.asarray(atol), itime)
    
    def setupGrid(self, grid):
        """Specify the grid."""
        return _cantera.domain_setupGrid(self._hndl, Numeric.asarray(grid))
    
    def setID(self, id):
        print 'id = ',id
        return _cantera.domain_setID(self._hndl, id)
    
    def setDesc(self, desc):
        return _cantera.domain_setDesc(self._hndl, desc)
    
    def grid(self, n):
        return _cantera.domain_grid(self._hndl, n)

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

    def _dict2array(self, d):
        a = zeros(self.nComponents(),'d')
        if d.has_key('default'):
            a += d['default']
        for k in d.keys():
            a[self.componentIndex(k)] = d[k]
        print a
        return a

    def _dict2arrays(self, d):
        a1 = zeros(self.nComponents(),'d')
        a2 = zeros(self.nComponents(),'d')        
        if d.has_key('default'):
            a1 += d['default'][0]
            a2 += d['default'][1]
            del d['default']
        for k in d.keys():
            a1[self.componentIndex(k)] = d[k][0]
            a2[self.componentIndex(k)] = d[k][1]            
        print a1, a2
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
            if opt == 'mdot':
                self.setMdot(v)
            elif opt == 'temperature':
                self.setTemperature(v)
            elif opt == 'mole_fractions':
                self.setMoleFractions(v)
            else:
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
        else:
            self._hndl = _cantera.surf_new()
        if id: self.setID(id)


    def setKineticsMgr(self, kin):
        _cantera.reactingsurf_setkineticsmgr(self._hndl,
                                             kin.kinetics_hndl())
    def enableCoverageEqs(self, onoff=1):
        _cantera.reactingsurf_enableCoverageEqs(self._hndl, onoff)


        
class AxisymmetricFlow(Domain1D):
    """An axisymmetric flow"""
    def __init__(self, id = 'axisymmetric_flow', gas = 'None'):
        Domain1D.__init__(self)
        iph = gas.thermo_hndl()
        ikin = gas.kinetics_hndl()
        itr = gas.transport_hndl()
        self._hndl = _cantera.stflow_new(iph, ikin, itr)
        if id: self.setID(id)

        self.setPressure(gas.pressure())
        
    def setPressure(self, p):
        """Set the pressure [Pa]. The pressure is a constant, since
        the governing equations are those appropriate for the
        low-Mach-number limit."""
        _cantera.stflow_setPressure(self._hndl, p)
        
    def setFixedTempProfile(self, temp):
        """Set the fixed temperature profile.
        This profile is used whenever the energy equation is disabled.
        """
        return _cantera.stflow_setFixedTempProfile(self._hndl, temp)
    
    def solveSpeciesEqs(self, flag):
        return _cantera.stflow_solveSpeciesEqs(self._hndl, flag)
    
    def solveEnergyEqn(self, flag):
        return _cantera.stflow_solveEnergyEqn(self._hndl, flag)


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
        
    def solve(self, loglevel=1, refine_grid=1):
        return _cantera.sim1D_solve(self._hndl, loglevel, refine_grid)
    
    def refine(self, loglevel=1):
        return _cantera.sim1D_refine(self._hndl, loglevel)
    
    def setRefineCriteria(self, dom, ratio = 10.0, slope = 0.8,
                          curve = 0.8, prune = 0.05):
        idom = dom.index()
        return _cantera.sim1D_setRefineCriteria(self._hndl,
                                                idom, ratio, slope, curve, prune)
    def save(self, fname, id, desc):
        return _cantera.sim1D_save(self._hndl, fname, id, desc)
    
    def restore(self, fname, id):
        return _cantera.sim1D_restore(self._hndl, fname, id)

    def writeStats(self):
        return _cantera.sim1D_writeStats(self._hndl)

    def domainIndex(self, name):
        return _cantera.sim1D_domainIndex(self._hndl, name)

    def value(self, dom, icomp, localPoint):
        idom = dom.index()
        return _cantera.sim1D_value(self._hndl, idom, icomp, localPoint)

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
