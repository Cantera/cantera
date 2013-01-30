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

    property n_components:
        """Number of solution components at each grid point."""
        def __get__(self):
            return self.domain.nComponents()

    property n_points:
        """Number of grid points belonging to this domain."""
        def __get__(self):
            return self.domain.nPoints()

    def component_name(self, int n):
        """Name of the nth component."""
        return pystr(self.domain.componentName(n))

    property component_names:
        """List of the names of all components of this domain."""
        def __get__(self):
            return [self.component_name(n) for n in range(self.n_components)]

    def component_index(self, str name):
        """Index of the component with name 'name'"""
        return self.domain.componentIndex(stringify(name))

    def set_bounds(self, *, default=None, Y=None, **kwargs):
        """
        Set the lower and upper bounds on the solution.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (lower_bound, upper_bound) tuples as
        the values.  The keyword *default* may be used to specify default
        bounds for all unspecified components. The keyword *Y* can be used to
        stand for all species mass fractions in flow domains.

        >>> d.set_bounds(default=(0, 1), Y=(-1.0e-5, 2.0))
        """
        if default is not None:
            for n in range(self.n_components):
                self.domain.setBounds(n, default[0], default[1])

        if Y is not None:
            for n in range(4, self.n_components):
                self.domain.setBounds(n, Y[0], Y[1])

        for name,(lower,upper) in kwargs.items():
            self.domain.setBounds(self.component_name(name), lower, upper)

    def set_steady_tolerances(self, *, default=None, Y=None, **kwargs):
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (rtol, atol) tuples as the values.
        The keyword *default* may be used to specify default bounds for all
        unspecified components. The keyword *Y* can be used to stand for all
        species mass fractions in flow domains.
        """
        self._set_tolerances(0, default, Y, kwargs)

    def set_transient_tolerances(self, *, default=None, Y=None, **kwargs):
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (rtol, atol) tuples as the values.
        The keyword *default* may be used to specify default bounds for all
        unspecified components. The keyword *Y* can be used to stand for all
        species mass fractions in flow domains.
        """
        self._set_tolerances(1, default, Y, kwargs)

    def _set_tolerances(self, is_transient, default, Y, components):
        if default is not None:
            for n in range(self.n_components):
                self.domain.setTolerances(n, default[0], default[1],
                                          is_transient)

        if Y is not None:
            for n in range(4, self.n_components):
                self.domain.setTolerances(n, Y[0], Y[1], is_transient)

        for name,(lower,upper) in components.items():
            self.domain.setTolerances(self.component_name(name),
                                      lower, upper, is_transient)

    def bounds(self, component):
        """
        Return the (lower, upper) bounds for a solution component.

        >>> d.bounds('T')
        (200.0, 5000.0)
        """
        n = self.component_index(component)
        return self.domain.lowerBound(n), self.domain.upperBound(n)

    def tolerances(self, component):
        """
        Return the (relative, absolute) error tolerances for a solution
        component.

        >>> rtol, atol = d.tolerances('u')
        """
        k = self.component_index(component)
        return self.domain.rtol(k), self.domain.atol(k)

    property grid:
        """ The grid for this domain """
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] grid = np.empty(self.n_points)
            cdef int i
            for i in range(self.n_points):
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
            if isinstance(X, (str, unicode)):
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

    property spread_rate:
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

    def set_kinetics(self, Kinetics kin):
        """Set the kinetics manager (surface reaction mechanism object)."""
        if kin.kinetics.type() not in (kinetics_type_interface,
                                       kinetics_type_edge):
            raise TypeError('Kinetics object must be derived from '
                            'InterfaceKinetics.')
        self.surf.setKineticsMgr(<CxxInterfaceKinetics*>kin.kinetics)

    def enable_coverage_equations(self, on=True):
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
        self.P = self.gas.P
        self.flow.solveEnergyEqn()

    property P:
        """ Pressure [Pa] """
        def __get__(self):
            return self.flow.pressure()
        def __set__(self, P):
            self.flow.setPressure(P)

    def set_transport(self, _SolutionBase phase):
        self.gas = phase
        self.flow.setTransport(deref(self.gas.transport))

    property soret_enabled:
        """
        Determines whether or not to include diffusive mass fluxes due to the
        Soret effect. Enabling this option works only when using the
        multicomponent transport model.
        """
        def __get__(self):
            return self.flow.withSoret()
        def __set__(self, enable):
            self.flow.enableSoret(<cbool>enable)

    property energy_enabled:
        """ Determines whether or not to solve the energy equation."""
        def __get__(self):
            return self.flow.doEnergy(0)
        def __set__(self, enable):
            if enable:
                self.flow.solveEnergyEqn()
            else:
                self.flow.fixTemperature()

    def set_fixed_temp_profile(self, pos, T):
        """Set the fixed temperature profile. This profile is used
        whenever the energy equation is disabled.

        :param pos:
            arrray of relative positions from 0 to 1
        :param temp:
            array of temperature values

        >>> d.set_fixed_temp_profile(array([0.0, 0.5, 1.0]),
        ...                          array([500.0, 1500.0, 2000.0])
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


cdef class FreeFlow(_FlowBase):
    def __cinit__(self, _SolutionBase thermo, *args, **kwargs):
        gas = getIdealGasPhase(thermo)
        self.flow = <CxxStFlow*>(new CxxFreeFlame(gas, thermo.n_species, 2))


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
        self.flow = <CxxStFlow*>(new CxxAxiStagnFlow(gas, thermo.n_species, 2))


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

    def domain_index(self, dom):
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
        idom = self.domain_index(dom)
        dom = self.domains[idom]
        if isinstance(comp, (str, unicode)):
            kcomp = dom.component_index(comp)
        else:
            kcomp = comp

        assert 0 <= kcomp < dom.n_components

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

    def set_value(self, domain, component, point, value):
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

    def work_value(self, domain, component, point):
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
        cdef np.ndarray[np.double_t, ndim=1] data = np.empty(dom.n_points)
        for j in range(dom.n_points):
            data[j] = self.sim.value(idom, kcomp, j)
        return data

    def set_profile(self, domain, component, positions, values):
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

        >>> s.set_profile(d, 'T', [0.0, 0.2, 1.0], [400.0, 800.0, 1500.0])
        """
        dom, comp = self._get_indices(domain, component)

        cdef vector[double] pos_vec, val_vec
        for p in positions:
            pos_vec.push_back(p)
        for v in values:
            val_vec.push_back(v)

        self.sim.setProfile(dom, comp, pos_vec, val_vec)

    def set_flat_profile(self, domain, component, value):
        """Set a flat profile for one component in one domain.

        :param domain:
            Domain1D object, name, or index
        :param component:
            component name or index
        :param v:
            value

        >>> s.set_flat_profile(d, 'u', -3.0)
        """
        dom, comp = self._get_indices(domain, component)
        self.sim.setFlatProfile(dom, comp, value)

    def show_solution(self):
        """ print the current solution. """
        if not self._initialized:
            self.set_initial_guess()
        self.sim.showSolution()

    def set_time_step(self, stepsize, n_steps):
        """Set the sequence of time steps to try when Newton fails.

        :param stepsize:
            initial time step size [s]
        :param n_steps:
            sequence of integer step numbers

        >>> s.set_time_step(1.0e-5, [1, 2, 5, 10])
        """
        cdef vector[int] data
        for n in n_steps:
            data.push_back(n)
        self.sim.setTimeStep(stepsize, data.size(), &data[0])

    def set_initial_guess(self):
        """
        Set the initial guess for the solution. Derived classes extend this
        function to set approximations for the temperature and composition
        profiles.
        """
        self._get_initial_solution()
        self._initialized = True

    def _get_initial_solution(self):
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
            self.set_initial_guess()
        self.sim.solve(loglevel, <cbool>refine_grid)

    def refine(self, loglevel=1):
        """
        Refine the grid, adding points where solution is not adequately
        resolved.
        """
        self.sim.refine(loglevel)

    def set_refine_criteria(self, domain, ratio=10.0, slope=0.8, curve=0.8,
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

        >>> s.set_refine_criteria(d, ratio=5.0, slope=0.2, curve=0.3, prune=0.03)
        """
        idom = self.domain_index(domain)
        self.sim.setRefineCriteria(idom, ratio, slope, curve, prune)

    def set_max_jac_age(self, ss_age, ts_age):
        """
        Set the maximum number of times the Jacobian will be used before it
        must be re-evaluated.

        :param ss_age:
            age criterion during steady-state mode
        :param ts_age:
            age criterion during time-stepping mode
        """
        self.sim.setJacAge(ss_age, ts_age)

    def set_time_step_factor(self, tfactor):
        """
        Set the factor by which the time step will be increased after a
        successful step, or decreased after an unsuccessful one.
        """
        self.sim.setTimeStepFactor(tfactor)

    def set_min_time_step(self, tsmin):
        """ Set the minimum time step. """
        self.sim.setMinTimeStep(tsmin)

    def set_max_time_step(self, tsmax):
        """ Set the maximum time step. """
        self.sim.setMaxTimeStep(tsmax)

    def set_fixed_temperature(self, T):
        """
        Set the temperature used to fix the spatial location of a freely
        propagating flame.
        """
        self.sim.setFixedTemperature(T)

    def save(self, filename='soln.xml', name='solution', description='none',
             loglevel=1):
        """
        Save the solution in XML format.

        >>> s.save(file='save.xml', name='energy_off',
        ...        description='solution with energy eqn. disabled')

        """
        self.sim.save(stringify(filename), stringify(name),
                      stringify(description), loglevel)

    def restore(self, filename='soln.xml', name='solution', loglevel=2):
        """Set the solution vector to a previously-saved solution.

        :param filename:
            solution file
        :param name:
            solution name within the file
        :param loglevel:
            Amount of logging information to display while restoring,
            from 0 (disabled) to 2 (most verbose).

        >>> s.restore(filename='save.xml', id='energy_off')
        """
        self.sim.restore(stringify(filename), stringify(name), loglevel)
        self._initialized = True

    def show_stats(self, print_time=True):
        """
        Show the statistics for the last solution.

        If invoked with no arguments or with a non-zero argument, the timing
        statistics will be printed. Otherwise, the timing will not be printed.
        """
        self.sim.writeStats(print_time)

    def __dealloc__(self):
        del self.sim


class FlameBase(Sim1D):
    """ Base class for flames with a single flow domain """

    def __init__(self, domains, gas, grid):
        """
        :param gas:
            object to use to evaluate all gas properties and reaction rates
        :param grid:
            array of initial grid points
        """
        self.flame.grid = grid
        super().__init__(domains)
        self.gas = gas
        self.flame.P = gas.P

    def set_refine_criteria(self, ratio=10.0, slope=0.8, curve=0.8, prune=0.0):
        super().set_refine_criteria(self.flame, ratio, slope, curve, prune)

    def set_profile(self, component, locations, values):
        super().set_profile(self.flame, component, locations, values)

    @property
    def transport_model(self):
        return self.gas.transport_model

    @transport_model.setter
    def transport_model(self, model):
        self.gas.transport_model = model
        self.flame.set_transport(self.gas)

    @property
    def energy_enabled(self):
        return self.flame.energy_enabled

    @energy_enabled.setter
    def energy_enabled(self, enable):
        self.flame.energy_enabled = enable

    @property
    def soret_enabled(self):
        return self.flame.soret_enabled

    @soret_enabled.setter
    def soret_enabled(self, enable):
        self.flame.soret_enabled = enable

    @property
    def grid(self):
        """ Array of grid point positions along the flame. """
        return self.flame.grid

    @property
    def P(self):
        return self.flame.P

    @P.setter
    def P(self, P):
        self.flame.P = P

    @property
    def T(self):
        """ Array containing the temperature [K] at each grid point. """
        return self.profile(self.flame, 'T')

    @property
    def u(self):
        """
        Array containing the velocity [m/s] normal to the flame at each point.
        """
        return self.profile(self.flame, 'u')

    @property
    def V(self):
        """
        Array containing the tangential velocity gradient [1/s] at each point.
        """
        return self.profile(self.flame, 'V')

    def solution(self, component, point=None):
        if point is None:
            return self.profile(self.flame, component)
        else:
            return self.value(self.flame, component, point)

    def set_gas_state(self, point):
        k0 = self.flame.component_index(self.gas.species_name(0))
        Y = [self.solution(k, point)
             for k in range(k0, k0 + self.gas.n_species)]
        self.gas.TPY = self.value(self.flame, 'T', point), self.P, Y

def _trim(docstring):
    """Remove block indentation from a docstring."""
    if not docstring:
        return ''
    lines = docstring.splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = 999
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < 999:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())

    # Return a single string, with trailing and leading blank lines stripped
    return '\n'.join(trimmed).strip('\n')

def _array_property(attr, size=None):
    """
    Generate a property that retrieves values at each point in the flame. The
    'size' argument is the attribute name of the gas object used to set the
    leading dimension of the resulting array.
    """
    def getter(self):
        if size is None:
            # 1D array for scalar property
            vals = np.empty(self.flame.n_points)
        else:
            # 2D array
            vals = np.empty((getattr(self.gas, size), self.flame.n_points))

        for i in range(self.flame.n_points):
            self.set_gas_state(i)
            vals[...,i] = getattr(self.gas, attr)

        return vals

    if size is None:
        extradoc = "\nReturns an array of length `n_points`."
    else:
        extradoc = "\nReturns an array of size `%s` x `n_points`." % size

    doc = _trim(getattr(Solution, attr).__doc__) + extradoc
    return property(getter, doc=doc)

# Add scalar properties to FlameBase
for attr in ['density', 'density_mass', 'density_mole', 'volume_mass',
             'volume_mole', 'int_energy_mole', 'int_energy_mass', 'h',
             'enthalpy_mole', 'enthalpy_mass', 's', 'entropy_mole',
             'entropy_mass', 'g', 'gibbs_mole', 'gibbs_mass', 'cv',
             'cv_mole', 'cv_mass', 'cp', 'cp_mole', 'cp_mass',
             'isothermal_compressibility', 'thermal_expansion_coeff',
             'viscosity', 'thermal_conductivity']:
    setattr(FlameBase, attr, _array_property(attr))
FlameBase.volume = _array_property('v') # avoid confusion with velocity gradient 'V'
FlameBase.int_energy = _array_property('u') # avoid collision with velocity 'u'

# Add properties with values for each species
for attr in ['X', 'Y', 'concentrations', 'partial_molar_enthalpies',
             'partial_molar_entropies', 'partial_molar_int_energies',
             'chemical_potentials', 'electrochemical_potentials', 'partial_molar_cp',
             'partial_molar_volumes', 'standard_enthalpies_RT',
             'standard_entropies_R', 'standard_int_energies_RT',
             'standard_gibbs_RT', 'standard_cp_R', 'creation_rates',
             'destruction_rates', 'net_production_rates', 'mix_diff_coeffs',
             'mix_diff_coeffs_mass', 'mix_diff_coeffs_mole', 'thermal_diff_coeffs']:
    setattr(FlameBase, attr, _array_property(attr, 'n_species'))

# Add properties with values for each reaction
for attr in ['forward_rates_of_progress', 'reverse_rates_of_progress', 'net_rates_of_progress',
             'equilibrium_constants', 'forward_rate_constants', 'reverse_rate_constants',
             'delta_enthalpy', 'delta_gibbs', 'delta_entropy',
             'delta_standard_enthalpy', 'delta_standard_gibbs',
             'delta_standard_entropy']:
    setattr(FlameBase, attr, _array_property(attr, 'n_reactions'))


class FreeFlame(FlameBase):
    """A freely-propagating flat flame."""

    def __init__(self, gas, grid):
        """
        A domain of type FreeFlow named 'flame' will be created to represent
        the flame. The three domains comprising the stack are stored as
        ``self.inlet``, ``self.flame``, and ``self.outlet``.
        """
        self.inlet = Inlet1D()
        self.inlet.name = 'reactants'
        self.outlet = Outlet1D()
        self.outlet.name = 'products'
        self.flame = FreeFlow(gas)
        self.flame.name = 'flame'

        super().__init__((self.inlet, self.flame, self.outlet), gas, grid)

    def set_initial_guess(self):
        """
        Set the initial guess for the solution. The adiabatic flame
        temperature and equilibrium composition are computed for the inlet gas
        composition. The temperature profile rises linearly over 20% of the
        domain width to Tad, then is flat. The mass fraction profiles are set
        similarly.
        """
        super().set_initial_guess()
        self.gas.TPY = self.inlet.T, self.P, self.inlet.Y
        Y0 = self.inlet.Y
        u0 = self.inlet.mdot/self.gas.density
        T0 = self.inlet.T

        # get adiabatic flame temperature and composition
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y
        u1 = self.inlet.mdot/self.gas.density

        locs = [0.0, 0.3, 0.5, 1.0]
        self.set_profile('u', locs, [u0, u0, u1, u1])
        self.set_profile('T', locs, [T0, T0, Teq, Teq])
        self.set_fixed_temperature(0.5 * (T0 + Teq))
        for n in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(n),
                             locs, [Y0[n], Y0[n], Yeq[n], Yeq[n]])


class BurnerFlame(FlameBase):
    """A burner-stabilized flat flame."""

    def __init__(self, gas, grid):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            Array of initial grid points

        A domain of class `AxisymmetricStagnationFlow` named ``flame`` will
        be created to represent the flame. The three domains comprising the
        stack are stored as ``self.burner``, ``self.flame``, and
        ``self.outlet``.
        """
        self.burner = Inlet1D()
        self.burner.name = 'burner'
        self.burner.T = gas.T
        self.outlet = Outlet1D()
        self.outlet.name = 'outlet'
        self.flame = AxisymmetricStagnationFlow(gas)
        self.flame.name = 'flame'

        super().__init__((self.burner, self.flame, self.outlet), gas, grid)

    def set_initial_guess(self):
        """
        Set the initial guess for the solution. The adiabatic flame
        temperature and equilibrium composition are computed for the burner
        gas composition. The temperature profile rises linearly in the first
        20% of the flame to Tad, then is flat. The mass fraction profiles are
        set similarly.
        """
        super().set_initial_guess()

        self.gas.TPY = self.burner.T, self.P, self.burner.Y
        Y0 = self.burner.Y
        u0 = self.burner.mdot/self.gas.density
        T0 = self.burner.T

        # get adiabatic flame temperature and composition
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y
        u1 = self.burner.mdot/self.gas.density

        locs = [0.0, 0.2, 1.0]
        self.set_profile('u', locs, [u0, u1, u1])
        self.set_profile('T', locs, [T0, Teq, Teq])
        for n in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(n),
                             locs, [Y0[n], Yeq[n], Yeq[n]])


class CounterflowDiffusionFlame(FlameBase):
    """ A counterflow diffusion flame """

    def __init__(self, gas, grid):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            Array of initial grid points

        A domain of class `AxisymmetricStagnationFlow` named ``flame`` will
        be created to represent the flame. The three domains comprising the
        stack are stored as ``self.fuel_inlet``, ``self.flame``, and
        ``self.oxidizer_inlet``.
        """
        self.fuel_inlet = Inlet1D()
        self.fuel_inlet.name = 'fuel_inlet'
        self.fuel_inlet.T = gas.T

        self.oxidizer_inlet = Inlet1D()
        self.oxidizer_inlet.name = 'oxidizer_inlet'
        self.oxidizer_inlet.T = gas.T

        self.flame = AxisymmetricStagnationFlow(gas)
        self.flame.name = 'flame'

        super().__init__((self.fuel_inlet, self.flame, self.oxidizer_inlet),
                         gas, grid)

    def set_initial_guess(self, fuel, oxidizer='O2', stoich=None):
        """
        Set the initial guess for the solution. The fuel species must be
        specified:

        >>> f.set_initial_guess(fuel='CH4')

        The oxidizer and corresponding stoichiometry must be specified if it
        is not 'O2'. The initial guess is generated by assuming infinitely-
        fast chemistry.
        """

        super().set_initial_guess()

        if stoich is None:
            if oxidizer == 'O2':
                stoich = 0.0
                if 'H' in self.gas.element_names:
                    stoich += 0.25 * self.gas.n_atoms(fuel, 'H')
                if 'C' in self.gas.element_names:
                    stoich += self.gas.n_atoms(fuel, 'C')
            else:
                raise Exception('oxidizer/fuel stoichiometric ratio must be '
                                'specified since the oxidizer is not O2')

        kFuel = self.gas.species_index(fuel)
        kOx = self.gas.species_index(oxidizer)

        s = stoich * self.gas.molecular_weights[kOx] / self.gas.molecular_weights[kFuel]
        phi = s * self.fuel_inlet.Y[kFuel] / self.oxidizer_inlet.Y[kOx]
        zst = 1.0 / (1.0 + phi)

        Yin_f = self.fuel_inlet.Y
        Yin_o = self.oxidizer_inlet.Y
        Yst = zst * Yin_f + (1.0 - zst) * Yin_o

        self.gas.TPY = self.fuel_inlet.T, self.P, Yin_f
        mdotf = self.fuel_inlet.mdot
        u0f = mdotf / self.gas.density
        T0f = self.fuel_inlet.T

        self.gas.TPY = self.oxidizer_inlet.T, self.P, Yin_o
        mdoto = self.oxidizer_inlet.mdot
        u0o = mdoto/self.gas.density
        T0o = self.oxidizer_inlet.T

        # get adiabatic flame temperature and composition
        Tbar = 0.5 * (T0f + T0o)
        self.gas.TPY = Tbar, self.P, Yst
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y

        # estimate strain rate
        zz = self.flame.grid
        dz = zz[-1] - zz[0]
        a = (u0o + u0f)/dz
        f = np.sqrt(a / (2.0 * self.gas.mix_diff_coeffs[kOx]))

        x0 = mdotf * dz / (mdotf + mdoto)
        nz = len(zz)

        Y = np.zeros((nz, self.gas.n_species))
        T = np.zeros(nz)
        for j in range(nz):
            x = zz[j]
            zeta = f * (x - x0)
            zmix = 0.5 * (1.0 - math.erf(zeta))
            if zmix > zst:
                Y[j] = Yeq + (Yin_f - Yeq) * (zmix - zst) / (1.0 - zst)
                T[j] = Teq + (T0f - Teq) * (zmix - zst) / (1.0 - zst)
            else:
                Y[j] = Yin_o + zmix * (Yeq - Yin_o) / zst
                T[j] = T0o + (Teq - T0o) * zmix / zst

        T[0] = T0f
        T[-1] = T0o
        zrel = zz/dz

        self.set_profile('u', [0.0, 1.0], [u0f, -u0o])
        self.set_profile('V', [0.0, x0/dz, 1.0], [0.0, a, 0.0])
        self.set_profile('T', zrel, T)
        for k,spec in enumerate(self.gas.species_names):
            self.set_profile(spec, zrel, Y[:,k])
