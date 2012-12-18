cdef class Domain1D:
    cdef CxxDomain1D* domain
    def __cinit__(self, *args, **kwargs):
        self.domain = NULL

    def __init__(self, *args, **kwargs):
        if self.domain is NULL:
            raise TypeError("Can't instantiate abstract class Domain1D.")

    property index:
        """
        Index of this domain in a stack. Returns -1 if this domain is not part
        of a stack.
        """
        def __get__(self):
            return self.domain.domainIndex()

    property nComponents:
        """Number of solution components at each grid point."""
        def __get__(self):
            return self.domain.nComponents()

    property nPoints:
        """Number of grid points belonging to this domain."""
        def __get__(self):
            return self.domain.nPoints()

    def componentName(self, int n):
        """Name of the nth component."""
        return pystr(self.domain.componentName(n))

    property componentNames:
        """List of the names of all components of this domain."""
        def __get__(self):
            return [self.componentName(n) for n in range(self.nComponents)]

    def componentIndex(self, str name):
        """Index of the component with name 'name'"""
        return self.domain.componentIndex(stringify(name))

    def setBounds(self, *, default=None, Y=None, **kwargs):
        """
        Set the lower and upper bounds on the solution.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (lower_bound, upper_bound) tuples as
        the values.  The keyword *default* may be used to specify default
        bounds for all unspecified components. The keyword *Y* can be used to
        stand for all species mass fractions in flow domains.

        >>> d.setBounds(default=(0, 1), Y=(-1.0e-5, 2.0))
        """
        if default is not None:
            for n in range(self.nComponents):
                self.domain.setBounds(n, default[0], default[1])

        if Y is not None:
            for n in range(4, self.nComponents):
                self.domain.setBounds(n, Y[0], Y[1])

        for name,(lower,upper) in kwargs.items():
            self.domain.setBounds(self.componentName(name), lower, upper)

    def setSteadyTolerances(self, *, default=None, Y=None, **kwargs):
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (rtol, atol) tuples as the values.
        The keyword *default* may be used to specify default bounds for all
        unspecified components. The keyword *Y* can be used to stand for all
        species mass fractions in flow domains.
        """
        self._setTolerances(0, default, Y, kwargs)

    def setTransientTolerances(self, *, default=None, Y=None, **kwargs):
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (rtol, atol) tuples as the values.
        The keyword *default* may be used to specify default bounds for all
        unspecified components. The keyword *Y* can be used to stand for all
        species mass fractions in flow domains.
        """
        self._setTolerances(1, default, Y, kwargs)

    def _setTolerances(self, isTransient, default, Y, components):
        if default is not None:
            for n in range(self.nComponents):
                self.domain.setTolerances(n, default[0], default[1],
                                          isTransient)

        if Y is not None:
            for n in range(4, self.nComponents):
                self.domain.setTolerances(n, Y[0], Y[1], isTransient)

        for name,(lower,upper) in components.items():
            self.domain.setTolerances(self.componentName(name),
                                      lower, upper, isTransient)

    def bounds(self, component):
        """
        Return the (lower, upper) bounds for a solution component.

        >>> d.bounds('T')
        (200.0, 5000.0)
        """
        n = self.componentIndex(component)
        return self.domain.lowerBound(n), self.domain.upperBound(n)

    def tolerances(self, component):
        """
        Return the (relative, absolute) error tolerances for a solution
        component.

        >>> rtol, atol = d.tolerances('u')
        """
        k = self.componentIndex(component)
        return self.domain.rtol(k), self.domain.atol(k)

    property grid:
        """ The grid for this domain """
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] grid = np.empty(self.nPoints)
            cdef int i
            for i in range(self.nPoints):
                grid[i] = self.domain.grid(i)
            return grid

        def __set__(self, grid):
            cdef np.ndarray[np.double_t, ndim=1] data = \
                np.ascontiguousarray(grid, dtype=np.double)
            self.domain.setupGrid(len(data), &data[0])

    property name:
        """ The name / id of this domain """
        def __get__(self):
            return pystr(self.domain.id())
        def __set__(self, name):
            self.domain.setID(stringify(name))

    property description:
        """ A description of this domain """
        def __get__(self):
            return pystr(self.domain.desc())
        def __set__(self, desc):
            self.domain.setDesc(stringify(desc))


cdef class Boundary1D(Domain1D):
    """ Base class for boundary domains. """
    cdef CxxBdry1D* boundary
    def __cinit__(self, *args, **kwargs):
        self.boundary = NULL

    def __init__(self, *args, **kwargs):
        if self.boundary is NULL:
            raise TypeError("Can't instantiate abstract class Boundary1D.")
        self.domain = <CxxDomain1D*>(self.boundary)
        Domain1D.__init__(self, *args, **kwargs)

    property T:
        """ The temperature [K] at this boundary. """
        def __get__(self):
            return self.boundary.temperature()
        def __set__(self, T):
            self.boundary.setTemperature(T)

    property mdot:
        """ The mass flow rate per unit area [kg/m^2] """
        def __get__(self):
            return self.boundary.mdot()
        def __set__(self, mdot):
            self.boundary.setMdot(mdot)

    property X:
        """ Species mole fractions at this boundary. """
        def __set__(self, X):
            cdef np.ndarray[np.double_t, ndim=1] data
            if isinstance(X, str):
                self.boundary.setMoleFractions(stringify(X))
            else:
                data = np.ascontiguousarray(X, dtype=np.double)
                self.boundary.setMoleFractions(&data[0])

    property Y:
        """ Species mass fractions at this boundary. """
        def __get__(self):
            cdef int nsp = self.boundary.nSpecies()
            cdef np.ndarray[np.double_t, ndim=1] Y = np.empty(nsp)
            cdef int k
            for k in range(nsp):
                Y[k] = self.boundary.massFraction(k)
            return Y


cdef class Inlet1D(Boundary1D):
    """
    A one-dimensional inlet. Note that an inlet can only be a terminal
    domain - it must be either the leftmost or rightmost domain in a
    stack.
    """
    cdef CxxInlet1D* inlet
    def __cinit__(self, *args, **kwargs):
        self.inlet = new CxxInlet1D()
        self.boundary = <CxxBdry1D*>(self.inlet)

    def __dealloc__(self):
        del self.inlet

    property spreadRate:
        def __get__(self):
            return self.inlet.spreadRate()
        def __set__(self, s):
            self.inlet.setSpreadRate(s)


cdef class Outlet1D(Boundary1D):
    """
    A one-dimensional outlet. An outlet imposes a zero-gradient boundary
    condition on the flow.
    """
    cdef CxxOutlet1D* outlet
    def __cinit__(self, *args, **kwargs):
        self.outlet = new CxxOutlet1D()
        self.boundary = <CxxBdry1D*>(self.outlet)

    def __dealloc__(self):
        del self.outlet


cdef class OutletReservoir1D(Boundary1D):
    """
    A one-dimensional outlet into a reservoir.
    """
    cdef CxxOutletRes1D* outlet
    def __cinit__(self, *args, **kwargs):
        self.outlet = new CxxOutletRes1D()
        self.boundary = <CxxBdry1D*>(self.outlet)

    def __dealloc__(self):
        del self.outlet


cdef class SymmetryPlane1D(Boundary1D):
    """A symmetry plane."""
    cdef CxxSymm1D* symm
    def __cinit__(self, *args, **kwargs):
        self.symm = new CxxSymm1D()
        self.boundary = <CxxBdry1D*>(self.symm)

    def __dealloc__(self):
        del self.symm


cdef class Surface1D(Boundary1D):
    """A solid surface."""
    cdef CxxSurf1D* surf
    def __cinit__(self, *args, **kwargs):
        self.surf = new CxxSurf1D()
        self.boundary = <CxxBdry1D*>(self.surf)

    def __dealloc__(self):
        del self.surf


cdef class ReactingSurface1D(Boundary1D):
    """A reacting solid surface."""
    cdef CxxReactingSurf1D* surf
    def __cinit__(self, *args, **kwargs):
        self.surf = new CxxReactingSurf1D()
        self.boundary = <CxxBdry1D*>(self.surf)

    def __dealloc__(self):
        del self.surf

    def setKinetics(self, Kinetics kin):
        """Set the kinetics manager (surface reaction mechanism object)."""
        if kin.kinetics.type() not in (kinetics_type_interface,
                                       kinetics_type_edge):
            raise TypeError('Kinetics object must be derived from '
                            'InterfaceKinetics.')
        self.surf.setKineticsMgr(<CxxInterfaceKinetics*>kin.kinetics)

    def enableCoverageEquations(self, on=True):
        """ Turn solving the surface coverage equations on or off. """
        self.surf.enableCoverageEquations(<cbool>on)


cdef class _FlowBase(Domain1D):
    """ Base class for 1D flow domains """
    cdef CxxStFlow* flow
    cdef _SolutionBase gas
    def __cinit__(self, *args, **kwargs):
        self.flow = NULL

    def __init__(self, _SolutionBase thermo, *args, **kwargs):
        self.domain = <CxxDomain1D*>(self.flow)
        super().__init__(*args, **kwargs)
        self.gas = thermo
        self.flow.setKinetics(deref(self.gas.kinetics))
        self.flow.setTransport(deref(self.gas.transport))
        self.flow.setPressure(self.gas.P)
        self.flow.solveEnergyEqn()

    def setPressure(self, P):
        self.flow.setPressure(P)

    def setTransport(self, _SolutionBase phase):
        self.gas = phase
        self.flow.setTransport(deref(self.gas.transport))

    property soretEnabled:
        """
        Determines whether or not to include diffusive mass fluxes due to the
        Soret effect. Enabling this option works only when using the
        multicomponent transport model.
        """
        def __get__(self):
            return self.flow.withSoret()
        def __set__(self, enable):
            self.flow.enableSoret(<cbool>enable)

    property energyEnabled:
        """ Determines whether or not to solve the energy equation."""
        def __get__(self):
            return self.flow.doEnergy(0)
        def __set__(self, enable):
            if enable:
                self.flow.solveEnergyEqn()
            else:
                self.flow.fixTemperature()

    def setFixedTempProfile(self, pos, T):
        """Set the fixed temperature profile. This profile is used
        whenever the energy equation is disabled.

        :param pos:
            arrray of relative positions from 0 to 1
        :param temp:
            array of temperature values

        >>> d.setFixedTempProfile(array([0.0, 0.5, 1.0]),
        ...                       array([500.0, 1500.0, 2000.0])
        """
        cdef vector[double] x, y
        for p in pos:
            x.push_back(p)
        for t in T:
            y.push_back(t)
        self.flow.setFixedTempProfile(x, y)

    def __dealloc__(self):
        del self.flow


cdef CxxIdealGasPhase* getIdealGasPhase(ThermoPhase phase) except *:
    if phase.thermo.eosType() != thermo_type_ideal_gas:
        raise TypeError('ThermoPhase object is not an IdealGasPhase')
    return <CxxIdealGasPhase*>(phase.thermo)


cdef class StagnationFlow(_FlowBase):
    def __cinit__(self, _SolutionBase thermo, *args, **kwargs):
        gas = getIdealGasPhase(thermo)
        self.flow = new CxxStFlow(gas, thermo.nSpecies(), 2)


cdef class FreeFlow(_FlowBase):
    def __cinit__(self, _SolutionBase thermo, *args, **kwargs):
        gas = getIdealGasPhase(thermo)
        self.flow = <CxxStFlow*>(new CxxFreeFlame(gas, thermo.nSpecies, 2))


cdef class AxisymmetricStagnationFlow(_FlowBase):
    """
    An axisymmetric flow domain.

    In an axisymmetric flow domain, the equations solved are the similarity
    equations for the flow in a finite-height gap of infinite radial extent.
    The solution variables are:

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

    It may be shown that if the boundary conditions on these variables are
    independent of radius, then a similarity solution to the exact governing
    equations exists in which these variables are all independent of radius.
    This solution holds only in in low-Mach-number limit, in which case
    (dP/dz) = 0, and lambda is a constant. (Lambda is treated as a spatially-
    varying solution variable for numerical reasons, but in the final solution
    it is always independent of z.) As implemented here, the governing
    equations assume an ideal gas mixture.  Arbitrary chemistry is allowed, as
    well as arbitrary variation of the transport properties.
    """
    def __cinit__(self, _SolutionBase thermo, *args, **kwargs):
        gas = getIdealGasPhase(thermo)
        self.flow = <CxxStFlow*>(new CxxAxiStagnFlow(gas, thermo.nSpecies, 2))


cdef class Sim1D:
    """
    Class Sim1D is a container for one-dimensional domains. It also holds the
    multi-domain solution vector, and controls the process of finding the
    solution.

    Domains are ordered left-to-right, with domain number 0 at the left.
    """
    cdef CxxSim1D* sim
    cdef readonly object domains
    cdef object _initialized

    def __cinit__(self, *args, **kwargs):
        self.sim = NULL

    def __init__(self, domains, *args, **kwargs):
        cdef vector[CxxDomain1D*] D
        cdef Domain1D d
        for d in domains:
            D.push_back(d.domain)

        self.sim = new CxxSim1D(D)
        self.domains = tuple(domains)

        self._initialized = False

    def domainIndex(self, dom):
        """
        Get the index of a domain, specified either by name or as a Domain1D
        object.
        """
        if isinstance(dom, Domain1D):
            idom = self.domains.index(dom)
        elif isinstance(dom, int):
            idom = dom
        else:
            idom = None
            for i,d in enumerate(self.domains):
                if d.name == dom:
                    idom = i
                    dom = d
            if idom is None:
                raise KeyError('Domain named "{}" not found.'.format(dom))

        assert 0 <= idom < len(self.domains)
        return idom

    def _get_indices(self, dom, comp):
        idom = self.domainIndex(dom)
        dom = self.domains[idom]
        if isinstance(comp, str):
            kcomp = dom.componentIndex(comp)
        else:
            kcomp = comp

        assert 0 <= kcomp < dom.nComponents

        return idom, kcomp

    def value(self, domain, component, point):
        """
        Solution value at one point

        :param domain:
            Domain1D object, name, or index
        :param component:
            component name or index
        :param point:
            grid point number within *domain* starting with 0 on the left

        >>> t = s.value('flow', 'T', 6)
        """
        dom, comp = self._get_indices(domain, component)
        return self.sim.value(dom, comp, point)

    def setValue(self, domain, component, point, value):
        """
        Set the value of one component in one domain at one point to 'value'.

        :param domain:
            Domain1D object, name, or index
        :param component:
            component name or index
        :param point:
            grid point number within *domain* starting with 0 on the left
        :param value:
            numerical value

        >>> s.set(d, 3, 5, 6.7)
        >>> s.set(1, 0, 5, 6.7)
        >>> s.set('flow', 'T', 5, 500)
        """
        dom, comp = self._get_indices(domain, component)
        self.sim.setValue(dom, comp, point, value)

    def workValue(self, domain, component, point):
        """
        Internal work array value at one point. After calling eval, this array
        contains the values of the residual function.

        :param domain:
            Domain1D object, name, or index
        :param component:
            component name or index
        :param point:
            grid point number in the domain, starting with zero at the left

        >>> t = s.value(flow, 'T', 6)
        """
        dom, comp = self._get_indices(domain, component)
        return self.sim.workValue(dom, comp, point)

    def profile(self, domain, component):
        """
        Spatial profile of one component in one domain.

        :param domain:
            Domain1D object, name, or index
        :param component:
            component name or index

        >>> T = s.profile(flow, 'T')
        """
        idom, kcomp = self._get_indices(domain, component)
        dom = self.domains[idom]
        cdef int j
        cdef np.ndarray[np.double_t, ndim=1] data = np.empty(dom.nPoints)
        for j in range(dom.nPoints):
            data[j] = self.sim.value(idom, kcomp, j)
        return data

    def setProfile(self, domain, component, positions, values):
        """
        Set an initial estimate for a profile of one component in one domain.

        :param domain:
            Domain1D object, name, or index
        :param component:
            component name or index
        :param positions:
            sequence of relative positions, from 0 on the left to 1 on the right
        :param values:
            sequence of values at the relative positions specified in *positions*

        >>> s.setProfile(d, 'T', [0.0, 0.2, 1.0], [400.0, 800.0, 1500.0])
        """
        dom, comp = self._get_indices(domain, component)

        cdef vector[double] pos_vec, val_vec
        for p in positions:
            pos_vec.push_back(p)
        for v in values:
            val_vec.push_back(v)

        self.sim.setProfile(dom, comp, pos_vec, val_vec)

    def setFlatProfile(self, domain, component, value):
        """Set a flat profile for one component in one domain.

        :param domain:
            Domain1D object, name, or index
        :param component:
            component name or index
        :param v:
            value

        >>> s.setFlatProfile(d, 'u', -3.0)
        """
        dom, comp = self._get_indices(domain, component)
        self.sim.setFlatProfile(dom, comp, value)

    def showSolution(self):
        """ print the current solution. """
        if not self._initialized:
            self.setInitialGuess()
        self.sim.showSolution()

    def setTimeStep(self, stepsize, nSteps):
        """Set the sequence of time steps to try when Newton fails.

        :param stepsize:
            initial time step size [s]
        :param nSteps:
            sequence of integer step numbers

        >>> s.setTimeStep(1.0e-5, [1, 2, 5, 10])
        """
        cdef vector[int] data
        for n in nSteps:
            data.push_back(n)
        self.sim.setTimeStep(stepsize, data.size(), &data[0])

    def setInitialGuess(self):
        """
        Set the initial guess for the solution. Derived classes extend this
        function to set approximations for the temperature and composition
        profiles.
        """
        self._getInitialSolution()
        self._initialized = True

    def _getInitialSolution(self):
        """
        Load the initial solution from each domain into the global solution
        vector.
        """
        self.sim.getInitialSoln()

    def solve(self, loglevel=1, refine_grid=True):
        """
        Solve the problem.

        :param loglevel:
            integer flag controlling the amount of diagnostic output. Zero
            suppresses all output, and 5 produces very verbose output.
        :param refine_grid:
            if True, enable grid refinement.
        """
        if not self._initialized:
            self.setInitialGuess()
        self.sim.solve(loglevel, <cbool>refine_grid)

    def refine(self, loglevel=1):
        """
        Refine the grid, adding points where solution is not adequately
        resolved.
        """
        self.sim.refine(loglevel)

    def setRefineCriteria(self, domain, ratio=10.0, slope=0.8, curve=0.8,
                          prune=0.05):
        """
        Set the criteria used to refine one domain.

        :param domain:
            domain object, index, or name
        :param ratio:
            additional points will be added if the ratio of the spacing on
            either side of a grid point exceeds this value
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

        >>> s.setRefineCriteria(d, ratio=5.0, slope=0.2, curve=0.3, prune=0.03)
        """
        idom = self.domainIndex(domain)
        self.sim.setRefineCriteria(idom, ratio, slope, curve, prune)

    def setMaxJacAge(self, ss_age, ts_age):
        """
        Set the maximum number of times the Jacobian will be used before it
        must be re-evaluated.

        :param ss_age:
            age criterion during steady-state mode
        :param ts_age:
            age criterion during time-stepping mode
        """
        self.sim.setJacAge(ss_age, ts_age)

    def setTimeStepFactor(self, tfactor):
        """
        Set the factor by which the time step will be increased after a
        successful step, or decreased after an unsuccessful one.
        """
        self.sim.setTimeStepFactor(tfactor)

    def setMinTimeStep(self, tsmin):
        """ Set the minimum time step. """
        self.sim.setMinTimeStep(tsmin)

    def setMaxTimeStep(self, tsmax):
        """ Set the maximum time step. """
        self.sim.setMaxTimeStep(tsmax)

    def setFixedTemperature(self, T):
        """
        Set the temperature used to fix the spatial location of a freely
        propagating flame.
        """
        self.sim.setFixedTemperature(T)

    def save(self, filename='soln.xml', name='solution', description='none'):
        """
        Save the solution in XML format.

        >>> s.save(file='save.xml', name='energy_off',
        ...        description='solution with energy eqn. disabled')

        """
        self.sim.save(stringify(filename), stringify(name),
                      stringify(description))

    def restore(self, filename='soln.xml', name='solution'):
        """Set the solution vector to a previously-saved solution.

        :param filename:
            solution file
        :param name:
            solution name within the file

        >>> s.restore(filename='save.xml', id='energy_off')
        """
        self.sim.restore(stringify(filename), stringify(name))
        self._initialized = True

    def showStats(self, printTime=True):
        """
        Show the statistics for the last solution.

        If invoked with no arguments or with a non-zero argument, the timing
        statistics will be printed. Otherwise, the timing will not be printed.
        """
        self.sim.writeStats(printTime)

    def __dealloc__(self):
        del self.sim
