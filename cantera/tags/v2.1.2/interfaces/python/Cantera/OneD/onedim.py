from Cantera import *
from Cantera import _cantera
from Cantera.num import asarray, zeros

_onoff = {'on':1, 'yes':1, 'off':0, 'no':0, 1:1, 0:0}

class Domain1D:
    """Base class for one-dimensional domains."""

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
        """List of the names of all components of this domain."""
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
        tuples as the values.  The keyword *default* may be used to
        specify default bounds for all unspecified components. The
        keyword *Y* can be used to stand for all species mass
        fractions in flow domains.

        >>> d.setBounds(default=(0, 1),
        ...             Y=(-1.0e-5, 2.0))
        """

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
        """Return the (lower, upper) bounds for a solution component.

        >>> d.bounds('T')
        (200.0, 5000.0)
        """

        ic = self.componentIndex(component)
        lower = _cantera.domain_lowerBound(self._hndl, ic)
        upper = _cantera.domain_upperBound(self._hndl, ic)
        return (lower, upper)

    def tolerances(self, component):
        """Return the (relative, absolute) error tolerances for
        a solution component.

        >>> (r, a) = d.tolerances('u')
        """
        ic = self.componentIndex(component)
        r = _cantera.domain_rtol(self._hndl, ic)
        a = _cantera.domain_atol(self._hndl, ic)
        return (r, a)

    def setTolerances(self, **tol):
        """Set the error tolerances. If *time* is present and
        non-zero, then the values entered will apply to the transient
        problem. Otherwise, they will apply to the steady-state
        problem.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (rtol, atol) tuples as the
        values.  The keyword *default* may be used to specify default
        bounds for all unspecified components. The keyword *Y* can be
        used to stand for all species mass fractions in flow domains.

        >>> d.setTolerances(Y=(1.0e-5, 1.0e-9),
        ...                 default=(1.0e-7, 1.0e-12),
        ...                 time=1)
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
            # print 'setting tol for ',b,' itime = ',itime
            _cantera.domain_setTolerances(self._hndl, n, d[b][0], d[b][1], itime)


    def setupGrid(self, grid):
        """Specify the grid.

        >>> d.setupGrid([0.0, 0.1, 0.2])

        """
        return _cantera.domain_setupGrid(self._hndl, asarray(grid))

    def setID(self, id):
        return _cantera.domain_setID(self._hndl, id)

    def setDesc(self, desc):
        """Set the description of this domain."""
        return _cantera.domain_setDesc(self._hndl, desc)

    def grid(self, n = -1):
        """ If *n* >= 0, return the value of the nth grid point
        from the left in this domain. If n is not supplied, return
        the entire grid.

        >>> z4 = d.grid(4)
        >>> z_array = d.grid()

        """
        if n >= 0:
            return _cantera.domain_grid(self._hndl, n)
        else:
            g = zeros(self.nPoints(),'d')
            for j in range(len(g)):
                g[j] = _cantera.domain_grid(self._hndl, j)
            return g

    def set(self, **options):
        """
        convenient function to invoke other methods.
        Parameters that can be set:

        grid, name, desc

        >>> d.set(name='flame', grid=z)
        """
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
            #elif opt == 'bounds':
            #    lower, upper = self._dict2arrays(v)
            #    self.setBounds(lower,upper)
            #elif opt == 'tol':
            #    self.setTolerances(v[0],v[1])
            #else:
            #    raise CanteraError('unknown attribute: '+opt)

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
    """Base class for boundary domains."""

    def __init__(self):
        Domain1D.__init__(self)

    def setMdot(self, mdot):
        """Set the mass flow rate per unit area [kg/m2]."""
        _cantera.bdry_setMdot(self._hndl, mdot)

    def setTemperature(self, t):
        """Set the temperature [K]"""
        _cantera.bdry_setTemperature(self._hndl, t)

    def setMoleFractions(self, x):
        """set the mole fraction values. """
        _cantera.bdry_setMoleFractions(self._hndl, x)

    def temperature(self):
        """Set the temperature [K]."""
        return _cantera.bdry_temperature(self._hndl)

    def massFraction(self, k):
        """The mass fraction of species k."""
        return _cantera.bdry_massFraction(self._hndl, k)

    def mdot(self):
        """The mass flow rate per unit area [kg/m2/s"""
        return _cantera.bdry_mdot(self._hndl)

    def set(self, **options):
        """Set parameters:
        mdot or massflux
        temperature or T
        mole_fractions or X
        """
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
        """Set the spead rate, defined as the value of V = v/r at the inlet."""
        _cantera.inlet_setSpreadRate(self._hndl, V0)


class Outlet(Bdry1D):
    """A one-dimensional outlet. An outlet imposes a
    zero-gradient boundary condition on the flow."""

    def __init__(self, id = 'outlet'):
        Bdry1D.__init__(self)
        self._hndl = _cantera.outlet_new()
        if id: self.setID(id)

class OutletRes(Bdry1D):
    """A one-dimensional outlet into a reservoir."""

    def __init__(self, id = 'outletres'):
        Bdry1D.__init__(self)
        self._hndl = _cantera.outletres_new()
        if id: self.setID(id)


class SymmPlane(Bdry1D):
    """A symmetry plane."""
    def __init__(self, id = 'symmetry_plane'):
        Bdry1D.__init__(self)
        self._hndl = _cantera.symm_new()
        if id: self.setID(id)

class Surface(Bdry1D):
    """A surface (possibly reacting)."""
    def __init__(self, id = 'surface', surface_mech = None):
        Bdry1D.__init__(self)
        if surface_mech:
            self._hndl = _cantera.reactingsurf_new()
            self.setKineticsMgr(surface_mech)
        else:
            self._hndl = _cantera.surf_new()
        if id: self.setID(id)


    def setKineticsMgr(self, kin):
        """Set the kinetics manager (surface reaction mechanism object)."""
        _cantera.reactingsurf_setkineticsmgr(self._hndl,
                                             kin.kinetics_hndl())

    def setCoverageEqs(self, onoff='on'):
        """Turn solving the surface coverage equations on or off."""
        if onoff == 'on':
            _cantera.reactingsurf_enableCoverageEqs(self._hndl, 1)
        else:
            _cantera.reactingsurf_enableCoverageEqs(self._hndl, 0)


class AxisymmetricFlow(Domain1D):
    """An axisymmetric flow domain.

    In an axisymmetric flow domain, the equations solved are the
    similarity equations for the flow in a finite-height gap of
    infinite radial extent. The solution variables are

    *u*
        axial velocity
    *V*
        radial velocity divided by radius
    *T*
        temperature
    *lambda*
        (1/r)(dP/dr)
    *Y_k*
        species mass fractions

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
    def __init__(self, id = 'axisymmetric_flow', gas = None, type = 1):
        Domain1D.__init__(self)
        iph = gas.thermo_hndl()
        ikin = gas.kinetics_hndl()
        itr = gas.transport_hndl()
        self._hndl = _cantera.stflow_new(iph, ikin, itr, type)
        if id: self.setID(id)
        self.setPressure(gas.pressure())
        self.solveEnergyEqn()

    def setPressure(self, p):
        """Set the pressure [Pa]. The pressure is a constant, since
        the governing equations are those for the low-Mach-number limit."""
        _cantera.stflow_setPressure(self._hndl, p)

    def setTransportModel(self, transp, withSoret = 0):
        """Set the transport model. The argument must be a transport
        manager for the 'gas' object."""
        itr = transp.transport_hndl()
        _cantera.stflow_setTransport(self._hndl, itr, withSoret)

    def enableSoret(self, withSoret = 1):
        """Include or exclude thermal diffusion (Soret effect) when computing
        diffusion velocities. If withSoret is not supplied or is positive,
        thermal diffusion is enabled; otherwise it is disabled."""
        _cantera.stflow_enableSoret(self._hndl, withSoret)

    def pressure(self):
        """Pressure [Pa]."""
        return _cantera.stflow_pressure(self._hndl)

    def setFixedTempProfile(self, pos, temp):
        """Set the fixed temperature profile.  This profile is used
        whenever the energy equation is disabled.

        :param pos:
            arrray of relative positions from 0 to 1
        :param temp:
            array of temperature values

        >>> d.setFixedTempProfile(array([0.0, 0.5, 1.0]),
        ...                       array([500.0, 1500.0, 2000.0])
        """
        return _cantera.stflow_setFixedTempProfile(self._hndl, pos, temp)

    def solveSpeciesEqs(self, flag = 1):
        """Enable or disable solving the species equations. If invoked
        with no arguments or with a non-zero argument, the species
        equations will be solved. If invoked with a zero argument,
        they will not be, and instead the species profiles will be
        held at their initial values. Default: species equations
        enabled."""
        return _cantera.stflow_solveSpeciesEqs(self._hndl, _onoff[flag])

    def solveEnergyEqn(self, flag = 1):
        """Enable or disable solving the energy equation. If invoked
        with no arguments or with a non-zero argument, the energy
        equations will be solved. If invoked with a zero argument,
        it will not be, and instead the temperature profiles will be
        held to the one specified by the call to :meth:`.setFixedTempProfile`.
        Default: energy equation enabled."""
        return _cantera.stflow_solveEnergyEqn(self._hndl, _onoff[flag])

    def set(self, **opt):
        """Set parameters.
        In addition to the parameters that may be set by Domain1D.set,
        this method can be used to set the pressure and energy flag

        >>> d.set(pressure=OneAtm, energy='on')
        """
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

    This class is largely a shadow class for C++ kernel class Sim1D.
    """
    def __init__(self, domains = None):
        self._hndl = 0
        nd = len(domains)
        hndls = zeros(nd,'i')
        for n in range(nd):
            hndls[n] = domains[n].domain_hndl()
        self._hndl = _cantera.sim1D_new(hndls)
        self._domains = domains
        self._initialized = False

    def __del__(self):
        _cantera.sim1D_del(self._hndl)

    def setValue(self, dom, comp, localPoint, value):
        """Set the value of one component in one domain at one point
        to 'value'.

        :param dom:
            domain object
        :param comp:
            component number
        :param localPoint:
            grid point number within domain *dom* starting with zero on the left
        :param value:
            numerical value

        >>> s.set(d, 3, 5, 6.7)
        """
        idom = dom.domain_hndl()
        _cantera.sim1D_setValue(self._hndl, idom,
                                comp, localPoint, value)

    def setProfile(self, dom, comp, pos, v):
        """Set an initial estimate for a profile of one component in
        one domain.

        :param dom:
            domain object
        :param comp:
            component name
        :param pos:
            sequence of relative positions, from 0 on the left to 1 on the right
        :param v:
            sequence of values at the relative positions specified in 'pos'

        >>> s.setProfile(d, 'T', [0.0, 0.2, 1.0], [400.0, 800.0, 1500.0])
        """

        idom = dom.index()
        icomp = dom.componentIndex(comp)
        _cantera.sim1D_setProfile(self._hndl, idom, icomp,
                                  asarray(pos), asarray(v))

    def setFlatProfile(self, dom, comp, v):
        """Set a flat profile for one component in one domain.

        :param dom:
            domain object
        :param comp:
            component name
        :param v:
            value

        >>> s.setFlatProfile(d, 'u', -3.0)
        """
        idom = dom.index()
        icomp = dom.componentIndex(comp)
        _cantera.sim1D_setFlatProfile(self._hndl, idom, icomp, v)

    def showSolution(self, fname='-'):
        """Show the current solution. If called with no argument,
        the solution is printed to the screen. If a filename is
        supplied, it is written to the file.

        >>> s.showSolution()
        >>> s.showSolution('soln.txt')
        """
        if not self._initialized:
            self.init()
        _cantera.sim1D_showSolution(self._hndl, fname)

    def setTimeStep(self, stepsize, nsteps):
        """Set the sequence of time steps to try when Newton fails.

        :param stepsize:
            initial time step size [s]
        :param nsteps:
            sequence of integer step numbers

        >>> s.setTimeStep(1.0e-5, [1, 2, 5, 10])
        """
        # 3/20/09
        # The use of asarray seems to set the nsteps array to be of
        # type double. This needs to be checked out further.
        # Probably a function of python version and Numerics version
        _cantera.sim1D_setTimeStep(self._hndl, stepsize, asarray(nsteps))

    def getInitialSoln(self):
        """Load the initial solution from each domain into the global
        solution vector."""
        _cantera.sim1D_getInitialSoln(self._hndl)

    def solve(self, loglevel=1, refine_grid=1):
        """Solve the problem.

        :param loglevel:
            integer flag controlling the amount of diagnostic output. Zero
            suppresses all output, and 5 produces very verbose output. Default: 1
        :param refine_grid:
            if non-zero, enable grid refinement."""

        return _cantera.sim1D_solve(self._hndl, loglevel, refine_grid)

    def refine(self, loglevel=1):
        """Refine the grid, adding points where solution is not
        adequately resolved."""
        return _cantera.sim1D_refine(self._hndl, loglevel)

    def setRefineCriteria(self, domain = None, ratio = 10.0, slope = 0.8,
                          curve = 0.8, prune = 0.05):
        """Set the criteria used to refine one domain.

        :param domain:
            domain object
        :param ratio:
            additional points will be added if the ratio of the spacing
            on either side of a grid point exceeds this value
        :param slope:
            maximum difference in value between two adjacent points, scaled by
            the maximum difference in the profile (0.0 < slope < 1.0). Adds
            points in regions of high slope.
        :param curve:
            maximum difference in slope between two adjacent intervals, scaled
            by the maximum difference in the profile (0.0 < curve < 1.0). Adds
            points in regions of high curvature.
        :param prune:
            if the slope or curve criteria are satisfied to the level of
            'prune', the grid point is assumed not to be needed and is removed.
            Set prune significantly smaller than 'slope' and 'curve'. Set to
            zero to disable pruning the grid.

        >>> s.setRefineCriteria(d, ratio=5.0, slope=0.2, curve=0.3,
        ...                     prune=0.03)
        """
        idom = domain.index()
        return _cantera.sim1D_setRefineCriteria(self._hndl,
                                                idom, ratio, slope, curve, prune)

    def setGridMin(self, domain, gridmin):
        """
        Set the minimum allowable grid spacing in a domain.

        :param domain:
            domain object
        :param gridmin:
            The minimum allowable grid spacing [m] for this domain
        """
        idom = domain.index()
        return _cantera.sim1D_setGridMin(self._hndl, idom, gridmin)

    def save(self, file = 'soln.xml', id = 'solution', desc = 'none'):
        """Save the solution in XML format.

        >>> s.save(file='save.xml', id='energy_off',
        ...        desc='solution with energy eqn. disabled')

        """
        return _cantera.sim1D_save(self._hndl, file, id, desc)

    def restore(self, file = 'soln.xml', id = 'solution'):
        """Set the solution vector to a previously-saved solution.

        :param file:
            solution file
        :param id:
            solution name within the file

        >>> s.restore(file = 'save.xml', id = 'energy_off')
        """
        self._initialized = True
        return _cantera.sim1D_restore(self._hndl, file, id)

    def showStats(self, printTime = 1):
        """Show the statistics for the last solution.
           If invoked with no arguments or with a non-zero argument, the
           timing statistics will be printed. If invoked with a zero argument,
           the timing will not be printed.
           Default: print timing enabled.
        """
        return _cantera.sim1D_writeStats(self._hndl, _onoff[printTime])

    def domainIndex(self, name):
        """Integer index of the domain with name 'name'"""
        return _cantera.sim1D_domainIndex(self._hndl, name)

    def value(self, domain, component, localPoint):
        """Solution value at one point.

        :param domain:
            domain object
        :param component:
            component name
        :param localPoint:
            grid point number in the domain, starting with zero at the left

        >>> t = s.value(flow, 'T', 6)
        """
        icomp = domain.componentIndex(component)
        idom = domain.index()
        return _cantera.sim1D_value(self._hndl, idom, icomp, localPoint)

    def profile(self, domain, component):
        """Spatial profile of one component in one domain.

        >>> print s.profile(flow, 'T')
        """
        np = domain.nPoints()
        x = zeros(np,'d')
        for n in range(np):
            x[n] = self.value(domain, component, n)
        return x

    def workValue(self, dom, icomp, localPoint):
        """Internal work array value at one point. After calling eval,
        this array contains the values of the residual function.

        :param domain:
            domain object
        :param component:
            component name
        :param localPoint:
            grid point number in the domain, starting with zero at the left

        >>> t = s.value(flow, 'T', 6)
        """
        idom = dom.index()
        return _cantera.sim1D_workValue(self._hndl, idom, icomp, localPoint)

    def eval(self, rdt, count=1):
        """Evaluate the residual function. If count = 0, do is 'silently',
        without adding to the function evaluation counter"""
        return _cantera.sim1D_eval(self._hndl, rdt, count)

    def setMaxJacAge(self, ss_age, ts_age):
        """Set the maximum number of times the Jacobian will be used
        before it must be re-evaluated.

        :param ss_age:
            age criterion during steady-state mode
        :param ts_age:
            age criterion during time-stepping mode
        """
        return _cantera.sim1D_setMaxJacAge(self._hndl, ss_age, ts_age)

    def timeStepFactor(self, tfactor):
        """Set the factor by which the time step will be increased
        after a successful step, or decreased after an unsuccessful one.

        >>> s.timeStepFactor(3.0)
        """
        return _cantera.sim1D_timeStepFactor(self._hndl, tfactor)

    def setTimeStepLimits(self, tsmin, tsmax):
        """Set the maximum and minimum time steps."""
        return _cantera.sim1D_setTimeStepLimits(self._hndl, tsmin, tsmax)

    def setFixedTemperature(self, temp):
        """This is a temporary fix."""
        _cantera.sim1D_setFixedTemperature(self._hndl, temp)

def clearDomains():
    """Clear all domains."""
    _cantera.domain_clear()

def clearSim1D():
    """Clear all stacks."""
    _cantera.sim1D_clear()
