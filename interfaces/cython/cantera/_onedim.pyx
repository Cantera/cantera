# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from .interrupts import no_op
import warnings
import numpy as np

from ._utils cimport stringify, pystr, anymap_to_py
from ._utils import CanteraError
from cython.operator import dereference as deref

cdef class Domain1D:
    _domain_type = "none"
    def __cinit__(self, _SolutionBase phase not None, *args, **kwargs):
        self.domain = NULL

    def __init__(self, phase, *args, **kwargs):
        if self.domain is NULL:
            raise TypeError("Can't instantiate abstract class Domain1D.")

        self.gas = phase
        self.set_default_tolerances()

    property phase:
        """
        Phase describing the domain (that is, a gas phase or surface phase).
        """
        def __get__(self):
            return self.gas

    property index:
        """
        Index of this domain in a stack. Returns -1 if this domain is not part
        of a stack.
        """
        def __get__(self):
            return self.domain.domainIndex()

    property domain_type:
        """
        String indicating the domain implemented.
        """
        def __get__(self):
            return pystr(self.domain.domainType())

    property n_components:
        """Number of solution components at each grid point."""
        def __get__(self):
            return self.domain.nComponents()

    property n_points:
        """Number of grid points belonging to this domain."""
        def __get__(self):
            return self.domain.nPoints()

    def _to_array(self, SolutionArrayBase dest, cbool normalize):
        """
        Retrieve domain data as a `SolutionArray` object. Service method used by
        `FlameBase.to_array`, which adds information not available in Cython.

        .. versionadded:: 3.0
        """
        dest._base = self.domain.toArray(normalize)
        dest.base = dest._base.get()
        return dest

    def _from_array(self, SolutionArrayBase arr):
        """
        Restore domain data from a `SolutionArray` object. Service method used by
        `FlameBase.from_array`.

        .. versionadded:: 3.0
        """
        self.domain.fromArray(arr._base)

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

    def _has_component(self, str name):
        """Check whether `Domain1D` has component"""
        return self.domain.hasComponent(stringify(name))

    @property
    def grid(self):
        """The grid for this domain."""
        cdef vector[double] grid_vec = self.domain.grid()
        return np.asarray(grid_vec)

    @grid.setter
    def grid(self, grid):
        cdef vector[double] grid_vec
        for g in grid:
            grid_vec.push_back(g)
        self.domain.setupGrid(grid_vec)

    def get_values(self, str component):
        """
        Retrieve spatial profile of a component.

        :param component:
            component name

        >>> T = d.get_values('T')

        .. versionadded:: 3.2
        """
        cdef vector[double] values
        values.resize(self.n_points)
        self.domain.getValues(stringify(component), values)
        return np.asarray(values)

    def set_values(self, str component, values):
        """
        Specify spatial profile of a component.

        :param component:
            component name
        :param values:
            array containing values

        >>> d.set_values('T', T)

        .. versionadded:: 3.2
        """
        cdef vector[double] values_vec
        for v in values:
            values_vec.push_back(v)
        self.domain.setValues(stringify(component), values)

    def set_profile(self, component, positions, values):
        """
        Set an initial estimate for a profile of one component in one domain.

        :param component:
            component name
        :param positions:
            sequence of relative positions, from 0 on the left to 1 on the right
        :param values:
            sequence of values at the relative positions specified in ``positions``

        >>> d.set_profile('T', [0.0, 0.2, 1.0], [400.0, 800.0, 1500.0])

        .. versionadded:: 3.2
        """
        cdef vector[double] pos_vec, val_vec
        for p in positions:
            pos_vec.push_back(p)
        for v in values:
            val_vec.push_back(v)

        self.domain.setProfile(stringify(component), pos_vec, val_vec)

    def set_flat_profile(self, component, value):
        """
        Set a flat profile for a component.

        :param component:
            component name
        :param v:
            value

        >>> d.set_flat_profile('u', -3.0)

        .. versionadded:: 3.2
        """
        self.domain.setFlatProfile(stringify(component), value)

    def set_bounds(self, *, default=None, Y=None, **kwargs):
        """
        Set the lower and upper bounds on the solution.

        The argument list should consist of keyword/value pairs, with
        component names as keywords and (lower bound, upper bound) tuples as
        the values.  The keyword ``default`` may be used to specify default
        bounds for all unspecified components. The keyword ``Y`` can be used to
        stand for all species mass fractions in flow domains.

        >>> d.set_bounds(default=(0, 1), Y=(-1.0e-5, 2.0))
        """
        if default is not None:
            for n in range(self.n_components):
                self.domain.setBounds(n, default[0], default[1])

        if Y is not None:
            k0 = self.component_index(self.gas.species_name(0))
            for n in range(k0, k0 + self.gas.n_species):
                self.domain.setBounds(n, Y[0], Y[1])

        for name,(lower,upper) in kwargs.items():
            self.domain.setBounds(self.component_index(name), lower, upper)

    def set_steady_tolerances(self, *, default=None, Y=None, abs=None, rel=None,
                              **kwargs):
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with component names as
        keywords and (relative tolerance, absolute tolerance) tuples as the values.
        The keyword ``default`` may be used to specify default bounds for all
        unspecified components. The keyword ``Y`` can be used to stand for all
        species mass fractions in flow domains. Alternatively, the keywords
        ``abs`` and ``rel`` can be used to specify arrays for the absolute and
        relative tolerances for each solution component.
        """
        self.have_user_tolerances = True
        if default is not None:
            self.domain.setSteadyTolerances(default[0], default[1])

        if abs is not None or rel is not None:
            if rel is None:
                rel = self.steady_reltol()
            if abs is None:
                abs = self.steady_abstol()
            assert len(abs) == len(rel) == self.n_components
            for n, (r, a) in enumerate(zip(rel, abs)):
                self.domain.setSteadyTolerances(r, a, n)

        if Y is not None:
            k0 = self.component_index(self.gas.species_name(0))
            for n in range(k0, k0 + self.gas.n_species):
                self.domain.setSteadyTolerances(Y[0], Y[1], n)

        for name, (rtol, atol) in kwargs.items():
            self.domain.setSteadyTolerances(rtol, atol,
                                            self.component_index(name))

    def set_transient_tolerances(self, *, default=None, Y=None, abs=None,
                                 rel=None, **kwargs):
        """
        Set the error tolerances for the steady-state problem.

        The argument list should consist of keyword/value pairs, with component names as
        keywords and (relative tolerance, absolute tolerance) tuples as the values.
        The keyword ``default`` may be used to specify default bounds for all
        unspecified components. The keyword ``Y`` can be used to stand for all
        species mass fractions in flow domains. Alternatively, the keywords
        ``abs`` and ``rel`` can be used to specify arrays for the absolute and
        relative tolerances for each solution component.
        """
        self.have_user_tolerances = True
        if default is not None:
            self.domain.setTransientTolerances(default[0], default[1])

        if abs is not None or rel is not None:
            if rel is None:
                rel = self.transient_reltol()
            if abs is None:
                abs = self.transient_abstol()
            assert len(abs) == len(rel) == self.n_components
            for n, (r, a) in enumerate(zip(rel, abs)):
                self.domain.setTransientTolerances(r, a, n)

        if Y is not None:
            k0 = self.component_index(self.gas.species_name(0))
            for n in range(k0, k0 + self.gas.n_species):
                self.domain.setTransientTolerances(Y[0], Y[1], n)

        for name, (rtol, atol) in kwargs.items():
            self.domain.setTransientTolerances(rtol, atol,
                                               self.component_index(name))

    def set_default_tolerances(self):
        """
        Set all tolerances to their default values
        """
        self.set_steady_tolerances(default=(1e-4, 1e-9))
        self.set_transient_tolerances(default=(1e-4, 1e-11))
        self.have_user_tolerances = False

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

    def steady_reltol(self, component=None):
        """
        Return the relative error tolerance for the steady state problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.steady_rtol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.steady_rtol(self.component_index(component))

    def steady_abstol(self, component=None):
        """
        Return the absolute error tolerance for the steady state problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.steady_atol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.steady_atol(self.component_index(component))

    def transient_reltol(self, component=None):
        """
        Return the relative error tolerance for the transient problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.transient_rtol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.transient_rtol(self.component_index(component))

    def transient_abstol(self, component=None):
        """
        Return the absolute error tolerance for the transient problem for a
        specified solution component, or all components if none is specified.
        """
        if component is None:
            return np.array([self.domain.transient_atol(n)
                             for n in range(self.n_components)])
        else:
            return self.domain.transient_atol(self.component_index(component))

    property name:
        """ The name / id of this domain """
        def __get__(self):
            return pystr(self.domain.id())
        def __set__(self, name):
            self.domain.setID(stringify(name))

    def __reduce__(self):
        raise NotImplementedError('Domain1D object is not picklable')

    def __copy__(self):
        raise NotImplementedError('Domain1D object is not copyable')

    property settings:
        """
        Return comprehensive dictionary describing type, name, and simulation settings
        that are specific to domain types.

        .. versionchanged:: 3.0

            Added missing domain-specific simulation settings and updated structure.
        """
        def __get__(self):
            cdef shared_ptr[CxxSolutionArray] arr
            arr = self.domain.toArray(False)
            return anymap_to_py(arr.get().meta())


cdef class Boundary1D(Domain1D):
    """
    Base class for boundary domains.

    :param phase:
        The (gas) phase corresponding to the adjacent flow domain
    """
    def __cinit__(self, _SolutionBase phase, *args, name="", **kwargs):
        if self._domain_type in {"None"}:
            self.boundary = NULL
        else:
            self._domain = CxxNewDomain1D(
                stringify(self._domain_type), phase._base, stringify(name))
            self.domain = self._domain.get()
            self.boundary = <CxxBoundary1D*>self.domain

    def __init__(self, phase, name=None):
        if self.boundary is NULL:
            raise TypeError("Can't instantiate abstract class Boundary1D.")
        Domain1D.__init__(self, phase, name=name)

    property T:
        """ The temperature [K] at this boundary. """
        def __get__(self):
            return self.boundary.temperature()
        def __set__(self, T):
            self.boundary.setTemperature(T)

    property mdot:
        """The mass flow rate per unit area [kg/s/m²]"""
        def __get__(self):
            return self.boundary.mdot()
        def __set__(self, mdot):
            self.boundary.setMdot(mdot)

    property X:
        """
        Species mole fractions at this boundary. May be set as either a string
        or as an array.
        """
        def __get__(self):
            self.gas.TPY = self.gas.T, self.gas.P, self.Y
            return self.gas.X

        def __set__(self, X):
            self.gas.TPX = None, None, X
            cdef np.ndarray[np.double_t, ndim=1] data = self.gas.X
            self.boundary.setMoleFractions(&data[0])

    property Y:
        """
        Species mass fractions at this boundary. May be set as either a string
        or as an array.
        """
        def __get__(self):
            cdef int nsp = self.boundary.nSpecies()
            cdef np.ndarray[np.double_t, ndim=1] Y = np.empty(nsp)
            cdef int k
            for k in range(nsp):
                Y[k] = self.boundary.massFraction(k)
            return Y

        def __set__(self, Y):
            self.gas.TPY = self.gas.T, self.gas.P, Y
            self.X = self.gas.X

    property spread_rate:
        """
        Get/set the tangential velocity gradient [1/s] at this boundary.
        """
        def __get__(self):
            return self.boundary.spreadRate()
        def __set__(self, s):
            self.boundary.setSpreadRate(s)


cdef class Inlet1D(Boundary1D):
    """
    A one-dimensional inlet. Note that an inlet can only be a terminal
    domain - it must be either the leftmost or rightmost domain in a stack.
    """
    _domain_type = "inlet"


cdef class Outlet1D(Boundary1D):
    """
    A one-dimensional outlet. An outlet imposes a zero-gradient boundary
    condition on the flow.
    """
    _domain_type = "outlet"


cdef class OutletReservoir1D(Boundary1D):
    """
    A one-dimensional outlet into a reservoir.
    """
    _domain_type = "outlet-reservoir"


cdef class SymmetryPlane1D(Boundary1D):
    """A symmetry plane."""
    _domain_type = "symmetry-plane"


cdef class Surface1D(Boundary1D):
    """A solid surface."""
    _domain_type = "surface"


cdef class ReactingSurface1D(Boundary1D):
    """A reacting solid surface.

    :param phase:
        The (surface) phase corresponding to the boundary

    .. versionchanged:: 3.0

        Starting in Cantera 3.0, parameter `phase` should reference surface instead of
        gas phase.
    """
    _domain_type = "reacting-surface"

    def __init__(self, _SolutionBase phase, name=None):
        gas = None
        for val in phase._adjacent.values():
            if val.phase_of_matter == "gas":
                gas = val
                break
        if gas is None:
            raise CanteraError("ReactingSurface1D needs an adjacent gas phase")
        super().__init__(gas, name=name)

        self.surface = phase

    property phase:
        """
        Get the `Interface` object representing species and reactions on the surface
        """
        def __get__(self):
            return self.surface

    property coverage_enabled:
        """Controls whether or not to solve the surface coverage equations."""
        def __set__(self, value):
            (<CxxReactingSurf1D*>self.domain).enableCoverageEquations(<cbool>value)
        def __get__(self):
            return (<CxxReactingSurf1D*>self.domain).coverageEnabled()


cdef class FlowBase(Domain1D):
    """ Base class for 1D flow domains """
    def __cinit__(self, _SolutionBase phase, *args, name="", **kwargs):
        self._domain = CxxNewDomain1D(
            stringify(self._domain_type), phase._base, stringify(name))
        self.domain = self._domain.get()
        self.flow = <CxxFlow1D*>self.domain

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.P = self.gas.P
        self.flow.solveEnergyEqn()

    def __getattr__(self, name):
        # used to access fields by alias rather than canonical names;
        # also provides access to species that are valid Python identifiers
        # (replicates SolutionArray behavior)
        component_name = name
        if (not self._has_component(name) and
            self._has_component(name.replace("_", "-"))):
            component_name = name.replace("_", "-")

        if self._has_component(component_name):
            return self.get_values(component_name)

        raise AttributeError(
            f"{self.__class__.__name__!r} object has no attribute {name!r}")

    @property
    def P(self):
        """Pressure [Pa]"""
        return self.flow.pressure()

    @P.setter
    def P(self, P):
        self.flow.setPressure(P)

    @property
    def T(self):
        """
        Array containing the temperature [K] at each grid point.

        .. versionadded:: 3.2
        """
        return self.get_values("T")

    @property
    def velocity(self):
        """
        Array containing the velocity [m/s] normal to the flame at each point.

        .. versionadded:: 3.2
        """
        return self.get_values("velocity")

    @property
    def spread_rate(self):
        """
        Array containing the tangential velocity gradient [1/s] (that is, radial
        velocity divided by radius) at each point. Note: This value is named
        ``spreadRate`` in the C++ code and is only defined for axisymmetric flows.

        .. versionadded:: 3.2
        """
        return self.get_values("spread_rate")

    @property
    def radial_pressure_gradient(self):
        """
        Array containing the radial pressure gradient (1/r)(dP/dr) [N/m⁴] at
        each point. Note: This value is named ``Lambda`` in the C++ code and is only
        defined for axisymmetric flows.

        .. versionadded:: 3.2
        """
        return self.get_values("Lambda")

    @property
    def electric_field(self):
        """
        Array containing the electric field strength at each point.
        Note: This value is named ``eField`` in the C++ code and is only defined if
        the transport model is ``ionized-gas``.

        .. versionadded:: 3.2
        """
        return self.get_values("eField")

    @property
    def oxidizer_velocity(self):
        """
        Array containing the oxidizer velocity (right boundary velocity) [m/s] at
        each point.
        Note: This value is named ``Uo`` in the C++ code and is only defined when using
        two-point control.

        .. versionchanged:: 3.2
        """
        return self.get_values("Uo")

    property transport_model:
        """
        Get/set the transport model used for calculating transport properties.

        .. versionadded:: 3.0
        """
        def __get__(self):
            return pystr(self.flow.transportModel())
        def __set__(self, model):
            self.flow.setTransportModel(stringify(model))
            # ensure that transport remains accessible
            self.gas.transport = self.gas.base.transport().get()

    def set_default_tolerances(self):
        """
        Set all tolerances to their default values
        """
        super().set_default_tolerances()
        if self.transport_model != "ionized-gas":
            return

        chargetol = {}
        for S in self.gas.species():
            if S.composition == {'E': 1.0}:
                chargetol[S.name] = (1e-5, 1e-20)
            elif S.charge != 0:
                chargetol[S.name] = (1e-5, 1e-16)
        self.set_steady_tolerances(**chargetol)
        self.set_transient_tolerances(**chargetol)
        self.have_user_tolerances = False

    property soret_enabled:
        """
        Determines whether or not to include diffusive mass fluxes due to the
        Soret effect. Enabling this option only works for multicomponent and
        mixture-averaged diffusion models.
        """
        def __get__(self):
            return self.flow.withSoret()
        def __set__(self, enable):
            self.flow.enableSoret(<cbool>enable)

    property flux_gradient_basis:
        """
        Get/Set whether or not species diffusive fluxes are computed with
        respect to their mass fraction gradients ('mass')
        or mole fraction gradients ('molar', default) when
        using the mixture-averaged transport model.
        """
        def __get__(self):
            if self.flow.fluxGradientBasis() == ThermoBasis.molar:
                return 'molar'
            else:
                return 'mass'
        def __set__(self, basis):
            if basis == 'molar':
                self.flow.setFluxGradientBasis(ThermoBasis.molar)
            elif basis == 'mass':
                self.flow.setFluxGradientBasis(ThermoBasis.mass)
            else:
                raise ValueError("Valid choices are 'mass' or 'molar'."
                                 " Got {!r}.".format(basis))

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
            array of relative positions from 0 to 1
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

    def get_settings3(self):
        """
        Temporary method returning new behavior of settings getter.

        .. versionadded:: 3.0
        """
        return self.settings

    property boundary_emissivities:
        """ Set/get boundary emissivities. """
        def __get__(self):
            return self.flow.leftEmissivity(), self.flow.rightEmissivity()
        def __set__(self, tuple epsilon):
            if len(epsilon) != 2:
                raise ValueError('Setting the boundary emissivities requires a '
                                 'tuple of length 2.')
            self.flow.setBoundaryEmissivities(epsilon[0], epsilon[1])

    property radiation_enabled:
        """ Determines whether or not to include radiative heat transfer """
        def __get__(self):
            return self.flow.radiationEnabled()
        def __set__(self, do_radiation):
            self.flow.enableRadiation(<cbool>do_radiation)

    property radiative_heat_loss:
        """
        Return radiative heat loss (only non-zero if radiation is enabled).
        """
        def __get__(self):
            cdef int j
            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.n_points)
            for j in range(self.n_points):
                data[j] = self.flow.radiativeHeatLoss(j)
            return data

    def set_free_flow(self):
        """
        Set flow configuration for freely-propagating flames, using an internal
        point with a fixed temperature as the condition to determine the inlet
        mass flux.
        """
        self.flow.setFreeFlow()

    def set_axisymmetric_flow(self):
        """
        Set flow configuration for axisymmetric counterflow or burner-stabilized
        flames, using specified inlet mass fluxes.
        """
        self.flow.setAxisymmetricFlow()

    @property
    def type(self):
        """
        Return the type of flow domain being represented.

        Examples:
        - ``free-flow``/``free-ion-flow``,
        - ``axisymmetric-flow``/``axisymmetric-ion-flow``,
        - ``unstrained-flow``/``unstrained-ion-flow``
        """
        return pystr(self.flow.type())

    @property
    def solving_stage(self):
        """
        Solving stage mode for handling ionized species (only relevant if transport
        model is ``ionized-gas``):

        - ``stage == 1``: the fluxes of charged species are set to zero
        - ``stage == 2``: the electric field equation is solved, and the drift flux for
          ionized species is evaluated

        .. deprecated:: 3.2

            To be removed after Cantera 3.2. Use the `electric_field_enabled` property
            instead.
        """
        return self.flow.getSolvingStage()

    @solving_stage.setter
    def solving_stage(self, stage):
        self.flow.setSolvingStage(stage)

    @property
    def electric_field_enabled(self):
        """
        Determines whether or not to solve the electric field equation (only relevant
        if transport model is ``ionized-gas``).
        """
        return self.flow.doElectricField()

    @electric_field_enabled.setter
    def electric_field_enabled(self, enable):
        if enable:
            self.flow.solveElectricField()
        else:
            self.flow.fixElectricField()

    property left_control_point_temperature:
        """ Get/Set the left control point temperature [K] """
        def __get__(self):
            return self.flow.leftControlPointTemperature()
        def __set__(self, T):
            self.flow.setLeftControlPointTemperature(T)

    property right_control_point_temperature:
        """ Get/Set the right control point temperature [K] """
        def __get__(self):
            return self.flow.rightControlPointTemperature()
        def __set__(self, T):
            self.flow.setRightControlPointTemperature(T)

    property left_control_point_coordinate:
        """ Get the left control point coordinate [m] """
        def __get__(self):
            return self.flow.leftControlPointCoordinate()

    property right_control_point_coordinate:
        """ Get the right control point coordinate [m] """
        def __get__(self):
            return self.flow.rightControlPointCoordinate()

    property two_point_control_enabled:
        """ Get/Set the state of the two-point flame control """
        def __get__(self):
            return self.flow.twoPointControlEnabled()
        def __set__(self, enable):
            self.flow.enableTwoPointControl(<cbool>enable)


cdef class FreeFlow(FlowBase):
    r"""A free flow domain. The equations solved are standard equations for adiabatic
    one-dimensional flow. The solution variables are:

    *velocity*
        axial velocity
    *T*
        temperature
    *Y_k*
        species mass fractions
    """
    _domain_type = "free-flow"


cdef class UnstrainedFlow(FlowBase):
    r"""An unstrained flow domain. The equations solved are standard equations for
    adiabatic one-dimensional flow. The solution variables are:

    *velocity*
        axial velocity
    *T*
        temperature
    *Y_k*
        species mass fractions
    """
    _domain_type = "unstrained-flow"


cdef class AxisymmetricFlow(FlowBase):
    r"""
    An axisymmetric flow domain. The equations solved are the similarity equations for
    the flow in a finite-height gap of infinite radial extent. The solution variables
    are:

    *velocity*
        axial velocity
    *spread_rate*
        radial velocity divided by radius
    *T*
        temperature
    *lambda*
        :math:`(1/r)(dP/dr)`
    *Y_k*
        species mass fractions

    It may be shown that if the boundary conditions on these variables are independent
    of radius, then a similarity solution to the exact governing equations exists in
    which these variables are all independent of radius. This solution holds only in
    the low-Mach-number limit, in which case :math:`(dP/dz) = 0`, and :math:`\Lambda` is
    a constant. (Lambda is treated as a spatially-varying solution variable for
    numerical reasons, but in the final solution it is always independent of :math:`z`.)
    As implemented here, the governing equations assume an ideal gas mixture. Arbitrary
    chemistry is allowed, as well as arbitrary variation of the transport properties.
    """
    _domain_type = "axisymmetric-flow"


cdef class Sim1D:
    """
    Class Sim1D is a container for one-dimensional domains. It also holds the
    multi-domain solution vector, and controls the process of finding the
    solution.

    Domains are ordered left-to-right, with domain number 0 at the left.
    """

    def __init__(self, domains, *args, **kwargs):
        cdef vector[shared_ptr[CxxDomain1D]] cxx_domains
        cdef Domain1D d
        for d in domains:
            cxx_domains.push_back(d._domain)

        self._sim = CxxNewSim1D(cxx_domains)
        self.sim = self._sim.get()
        self.domains = tuple(domains)
        self.set_interrupt(no_op)
        self._initialized = False
        self._initial_guess_args = ()
        self._initial_guess_kwargs = {}

    def phase(self, domain=None):
        """
        Return phase describing a domain (that is, a gas phase or surface phase).

        :param domain: Index of domain within `Sim1D.domains` list; the default
            is to return the phase of the parent `Sim1D` object.
        """
        if domain is None:
            return self.gas

        dom = self.domains[self.domain_index(domain)]
        return dom.phase

    def set_interrupt(self, f):
        """
        Set an interrupt function to be called each time that :ct:`OneDim::eval` is
        called. The signature of ``f`` is ``float f(float)``. The default
        interrupt function is used to trap `KeyboardInterrupt` exceptions so
        that ``ctrl-c`` can be used to break out of the C++ solver loop.
        """
        if f is None:
            self.sim.setInterrupt(NULL)
            self._interrupt = None
            return

        if not isinstance(f, Func1):
            f = Func1(f)
        self._interrupt = f
        self.sim.setInterrupt(self._interrupt.func)

    def set_time_step_callback(self, f):
        """
        Set a callback function to be called after each successful timestep.
        The signature of ``f`` is ``float f(float)``. The argument passed to ``f`` is
        the size of the timestep. The output is ignored.
        """
        if f is None:
            self.sim.setTimeStepCallback(NULL)
            self._time_step_callback = None
            return

        if not isinstance(f, Func1):
            f = Func1(f)
        self._time_step_callback = f
        self.sim.setTimeStepCallback(self._time_step_callback.func)

    def set_steady_callback(self, f):
        """
        Set a callback function to be called after each successful steady-state
        solve, before regridding. The signature of ``f`` is ``float f(float)``. The
        argument passed to ``f`` is 0.0 and the output is ignored.
        """
        if f is None:
            self.sim.setSteadyCallback(NULL)
            self._steady_callback = None
            return

        if not isinstance(f, Func1):
            f = Func1(f)
        self._steady_callback = f
        self.sim.setSteadyCallback(self._steady_callback.func)

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
                raise KeyError('Domain named "{0}" not found.'.format(dom))

        assert 0 <= idom < len(self.domains)
        return idom

    def _get_indices(self, dom, comp):
        idom = self.domain_index(dom)
        dom = self.domains[idom]
        if isinstance(comp, (str, bytes)):
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
            grid point number within ``domain`` starting with 0 on the left

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
            grid point number within ``domain`` starting with 0 on the left
        :param value:
            numerical value

        >>> s.set(d, 3, 5, 6.7)
        >>> s.set(1, 0, 5, 6.7)
        >>> s.set('flow', 'T', 5, 500)
        """
        dom, comp = self._get_indices(domain, component)
        self.sim.setValue(dom, comp, point, value)

    def eval(self, rdt=0.0):
        """
        Evaluate the governing equations using the current solution estimate,
        storing the residual in the array which is accessible with the
        `work_value` function.

        :param rdt:
           Reciprocal of the time-step
        """
        self.sim.eval(rdt)

    def work_value(self, domain, component, point):
        """
        Internal work array value at one point. After calling `eval`, this array
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

        .. deprecated:: 3.2

            To be removed after Cantera 3.2. Replaceable by `Domain1D.get_values`.
        """
        warnings.warn("To be removed after Cantera 3.2. Replaceable by "
                      "'Domain1D.get_values'.", DeprecationWarning)
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
            sequence of values at the relative positions specified in ``positions``

        >>> s.set_profile(d, 'T', [0.0, 0.2, 1.0], [400.0, 800.0, 1500.0])

        .. deprecated:: 3.2

            To be removed after Cantera 3.2. Replaceable by Domain1D.set_profile.
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

        .. deprecated:: 3.2

            To be removed after Cantera 3.2. Replaceable by Domain1D.set_flat_profile.
        """
        dom, comp = self._get_indices(domain, component)
        self.sim.setFlatProfile(dom, comp, value)

    def show(self):
        """ print the current solution. """
        if not self._initialized:
            self.set_initial_guess()
        self.sim.show()

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

    property max_time_step_count:
        """
        Get/Set the maximum number of time steps allowed before reaching the
        steady-state solution
        """
        def __get__(self):
            return self.sim.maxTimeStepCount()
        def __set__(self, nmax):
            self.sim.setMaxTimeStepCount(nmax)

    def set_jacobian_perturbation(self, relative, absolute, threshold):
        """
        Configure perturbations used to evaluate finite difference Jacobian

        :param relative:
            Relative perturbation (multiplied by the absolute value of each component).
            Default ``1.0e-5``.
        :param absolute:
            Absolute perturbation (independent of component value). Default ``1.0e-10``.
        :param threshold:
            Threshold below which to exclude elements from the Jacobian. Default ``0.0``.
        """
        self.sim.setJacobianPerturbation(relative, absolute, threshold)

    @property
    def linear_solver(self):
        """
        Get/Set the the linear solver used to hold the Jacobian matrix and solve linear
        systems as part of each Newton iteration. The default is a banded, direct
        solver. See :ref:`sec-python-jacobians` for available solvers.
        """
        return SystemJacobian.wrap(self.sim.linearSolver())

    @linear_solver.setter
    def linear_solver(self, SystemJacobian precon):
        self.sim.setLinearSolver(precon._base)

    def set_initial_guess(self, *args, **kwargs):
        """
        Store arguments for initial guess and prepare storage for solution.
        """
        self._initial_guess_args = args
        self._initial_guess_kwargs = kwargs
        self._get_initial_solution()
        self._initialized = True

    def _get_initial_solution(self):
        """
        Load the initial solution from each domain into the global solution
        vector.
        """
        self.sim.resize()
        self.sim.getInitialSoln()

    def extinct(self):
        """
        Method overloaded for some flame types to indicate if the flame has been
        extinguished. Base class method always returns 'False'
        """
        return False

    def solve(self, loglevel=1, refine_grid=True, auto=False):
        """
        Solve the problem.

        :param loglevel:
            integer flag controlling the amount of diagnostic output. Zero
            suppresses all output, and 5 produces very verbose output.
        :param refine_grid:
            if True, enable grid refinement.
        :param auto: if True, sequentially execute the different solution stages
            and attempt to automatically recover from errors. Attempts to first
            solve on the initial grid with energy enabled. If that does not
            succeed, a fixed-temperature solution will be tried followed by
            enabling the energy equation, and then with grid refinement enabled.
            If non-default tolerances have been specified or multicomponent
            transport is enabled, an additional solution using these options
            will be calculated.
        """

        if not auto:
            if not self._initialized:
                self.set_initial_guess()
            self.sim.solve(loglevel, <cbool>refine_grid)
            return

        def set_transport(multi):
            for dom in self.domains:
                if isinstance(dom, FlowBase):
                    dom.transport_model = multi

        # Do initial solution steps with default tolerances
        have_user_tolerances = any(dom.have_user_tolerances for dom in self.domains)
        if have_user_tolerances:
            # Save the user-specified tolerances
            atol_ss_final = [dom.steady_abstol() for dom in self.domains]
            rtol_ss_final = [dom.steady_reltol() for dom in self.domains]
            atol_ts_final = [dom.transient_abstol() for dom in self.domains]
            rtol_ts_final = [dom.transient_reltol() for dom in self.domains]

        def restore_tolerances():
            if have_user_tolerances:
                for i in range(len(self.domains)):
                    self.domains[i].set_steady_tolerances(abs=atol_ss_final[i],
                                                        rel=rtol_ss_final[i])
                    self.domains[i].set_transient_tolerances(abs=atol_ts_final[i],
                                                            rel=rtol_ts_final[i])

        for dom in self.domains:
            dom.set_default_tolerances()

        # Do initial steps without Soret diffusion
        soret_doms = [dom for dom in self.domains if getattr(dom, 'soret_enabled', False)]

        def set_soret(soret):
            for dom in soret_doms:
                dom.soret_enabled = soret

        set_soret(False)

        # Do initial solution steps without multicomponent transport
        transport = self.transport_model
        solve_multi = self.transport_model == 'multicomponent'
        if solve_multi:
            set_transport('mixture-averaged')

        def log(msg, *args):
            if loglevel:
                print('\n{:*^78s}'.format(' ' + msg.format(*args) + ' '))

        flow_domains = [D for D in self.domains if isinstance(D, FlowBase)]
        zmin = [D.grid[0] for D in flow_domains]
        zmax = [D.grid[-1] for D in flow_domains]

        # 'data' entry is used for restart
        data = self._initial_guess_kwargs.get('data')
        if data is not None:
           nPoints = [len(flow_domains[0].grid)]
        else:
           nPoints = [len(flow_domains[0].grid), 12, 24, 48]

        for N in nPoints:
            for i,D in enumerate(flow_domains):
                if N > self.get_max_grid_points(D):
                    raise CanteraError('Maximum number of grid points exceeded')

                if N != len(D.grid):
                    D.grid = np.linspace(zmin[i], zmax[i], N)

            if data is None:
                self.set_initial_guess(*self._initial_guess_args,
                                       **self._initial_guess_kwargs)

            # Try solving with energy enabled, which usually works
            log('Solving on {} point grid with energy equation enabled', N)
            self.energy_enabled = True
            try:
                self.sim.solve(loglevel, <cbool>False)
                solved = True
            except CanteraError as e:
                log(str(e))
                solved = False
            except Exception as e:
                # restore settings before re-raising exception
                set_transport(transport)
                set_soret(True)
                restore_tolerances()
                raise e

            # If initial solve using energy equation fails, fall back on the
            # traditional fixed temperature solve followed by solving the energy
            # equation
            if not solved or self.extinct():
                if self.extinct():
                    self.set_initial_guess(*self._initial_guess_args,
                                           **self._initial_guess_kwargs)
                log('Initial solve failed; Retrying with energy equation disabled')
                self.energy_enabled = False
                try:
                    self.sim.solve(loglevel, <cbool>False)
                    solved = True
                except CanteraError as e:
                    log(str(e))
                    solved = False
                except Exception as e:
                    # restore settings before re-raising exception
                    set_transport(transport)
                    set_soret(True)
                    restore_tolerances()
                    raise e

                if solved:
                    log('Solving on {} point grid with energy equation re-enabled', N)
                    self.energy_enabled = True
                    try:
                        self.sim.solve(loglevel, <cbool>False)
                        solved = True
                    except CanteraError as e:
                        log(str(e))
                        solved = False
                    except Exception as e:
                        # restore settings before re-raising exception
                        set_transport(transport)
                        set_soret(True)
                        restore_tolerances()
                        raise e

            if solved and not self.extinct() and refine_grid:
                # Found a non-extinct solution on the fixed grid
                log('Solving with grid refinement enabled')
                try:
                    self.sim.solve(loglevel, <cbool>True)
                    solved = True
                except CanteraError as e:
                    log(str(e))
                    solved = False
                except Exception as e:
                    # restore settings before re-raising exception
                    set_transport(transport)
                    set_soret(True)
                    restore_tolerances()
                    raise e

                if solved and not self.extinct():
                    # Found a non-extinct solution on the refined grid
                    break

            if self.extinct():
                log('Flame is extinct on {} point grid', N)

            if not refine_grid:
                break

        if not solved:
            raise CanteraError('Could not find a solution for the 1D problem')

        if solve_multi:
            log('Solving with multicomponent transport')
            set_transport('multicomponent')

        if soret_doms:
            log('Solving with Soret diffusion')
            set_soret(True)

        if have_user_tolerances:
            log('Solving with user-specified tolerances')
            restore_tolerances()

        # Final call with expensive options enabled
        if have_user_tolerances or solve_multi or soret_doms:
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

    def get_refine_criteria(self, domain):
        """
        Get a dictionary of the criteria used to refine one domain. The items in
        the dictionary are the ``ratio``, ``slope``, ``curve``, and ``prune``,
        as defined in `~Sim1D.set_refine_criteria`.

        :param domain:
            domain object, index, or name

        >>> s.set_refine_criteria(d, ratio=5.0, slope=0.2, curve=0.3, prune=0.03)
        >>> s.get_refine_criteria(d)
        {'ratio': 5.0, 'slope': 0.2, 'curve': 0.3, 'prune': 0.03}
        """
        idom = self.domain_index(domain)
        c = self.sim.getRefineCriteria(idom)
        return {'ratio': c[0], 'slope': c[1], 'curve': c[2], 'prune': c[3]}

    def set_grid_min(self, dz, domain=None):
        """
        Set the minimum grid spacing on ``domain``. If ``domain`` is `None`, then
        set the grid spacing for all domains.
        """
        if domain is None:
            idom = -1
        else:
            idom = self.domain_index(domain)
        self.sim.setGridMin(idom, dz)

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

    property fixed_temperature:
        """
        Set the temperature used to fix the spatial location of a freely
        propagating flame.
        """
        def __get__(self):
            return self.sim.fixedTemperature()
        def __set__(self, T):
            self.sim.setFixedTemperature(T)

    property fixed_temperature_location:
        """
        Return the location of the point where temperature is fixed for a freely
        propagating flame.
        """
        def __get__(self):
            return self.sim.fixedTemperatureLocation()

    def set_left_control_point(self, T):
        """
        Set the left control point using the specified temperature. This user-provided
        temperature will be used to locate the closest grid point to that temperature,
        which will serve to locate the left control point's coordinate. Starting from
        the left boundary, the first grid point that is equal to or exceeds the
        specified temperature will be used to locate the left control point's
        coordinate.
        """
        self.sim.setLeftControlPoint(T)

    def set_right_control_point(self, T):
        """
        Set the right control point using a specified temperature. This user-provided
        temperature will be used to locate the closest grid point to that temperature,
        which will serve to locate the right control point's coordinate.Starting from
        the right boundary, the first grid point that is equal to or exceeds the
        specified temperature will be used to locate the right control point's
        coordinate.
        """
        self.sim.setRightControlPoint(T)

    def save(self, filename='soln.yaml', name='solution', description=None,
             loglevel=None, *, overwrite=False, compression=0, basis=None):
        """
        Save current simulation data to a data file (CSV, YAML or HDF).

        In order to save the content of the current object, individual domains are
        converted to `SolutionArray` objects and saved using the `~SolutionArray.save`
        method. For HDF and YAML output, all domains are written to a single container
        file with shared header information. Simulation settings of individual domains
        are preserved as meta data of the corresponding `SolutionArray` objects.
        For CSV files, only state and auxiliary data of the main 1D domain are saved.

        The complete state of the current object can be restored from HDF and YAML
        container files using the `restore` method, while individual domains can be
        loaded using `SolutionArray.restore` for further analysis. While CSV files do
        not contain complete information, they can be used for setting initial states
        of individual simulation objects (example: `~FreeFlame.set_initial_guess`).

        :param filename:
            Name of output file (CSV, YAML or HDF)
        :param name:
            Identifier of storage location within the container file; this node/group
            contains header information and multiple subgroups holding domain-specific
            `SolutionArray` data (YAML/HDF only).
        :param description:
            Custom comment describing the dataset to be stored (YAML/HDF only).
        :param overwrite:
            Force overwrite if file and/or data entry exists; optional (default=`False`)
        :param compression:
            Compression level (0-9); optional (default=0; HDF only)
        :param basis:
            Output mass (``Y``/``mass``) or mole (``Y``/``mass``) fractions;
            if not specified (`None`), the native basis of the underlying `ThermoPhase`
            manager is used.

        >>> s.save(filename='save.yaml', name='energy_off',
        ...        description='solution with energy eqn. disabled')

        .. versionchanged:: 3.0
            Argument loglevel is no longer supported
        """
        if loglevel is not None:
            warnings.warn("Sim1D.save: Argument 'loglevel' is deprecated and will be "
                "ignored.", DeprecationWarning)
        self.sim.save(stringify(str(filename)), stringify(name),
                      stringify(description), overwrite, compression, stringify(basis))

    def restore(self, filename='soln.yaml', name='solution', loglevel=None):
        """Retrieve data and settings from a previously saved simulation.

        This method restores a simulation object from YAML or HDF data previously saved
        using the `save` method.

        :param filename:
            Name of container file (YAML or HDF)
        :param name:
            Identifier of location within the container file; this node/group contains
            header information and subgroups with domain-specific `SolutionArray` data
        :param loglevel:
            Amount of logging information to display while restoring,
            from 0 (disabled) to 2 (most verbose).
        :return:
            Dictionary containing header information

        >>> s.restore(filename='save.yaml', name='energy_off')

        .. versionchanged:: 3.0
            Implemented return value for meta data; loglevel is no longer supported
        """
        if loglevel is not None:
            warnings.warn("Sim1D.restore: Argument 'loglevel' is deprecated and will be"
                 " ignored.", DeprecationWarning)
        cdef CxxAnyMap header
        header = self.sim.restore(stringify(str(filename)), stringify(name))
        self._initialized = True
        return anymap_to_py(header)

    def restore_time_stepping_solution(self):
        """
        Set the current solution vector to the last successful time-stepping
        solution. This can be used to examine the solver progress after a failed
        integration.
        """
        self.sim.restoreTimeSteppingSolution()

    def restore_steady_solution(self):
        """
        Set the current solution vector to the last successful steady-state
        solution. This can be used to examine the solver progress after a
        failure during grid refinement.
        """
        self.sim.restoreSteadySolution()

    def show_stats(self, print_time=True):
        """
        Show the statistics for the last solution.

        If invoked with no arguments or with a non-zero argument, the timing
        statistics will be printed. Otherwise, the timing will not be printed.
        """
        self.sim.writeStats(print_time)

    def clear_stats(self):
        """
        Clear solver statistics.
        """
        self.sim.clearStats()

    def solve_adjoint(self, perturb, n_params, dgdx, g=None, dp=1e-5):
        """
        Find the sensitivities of an objective function using an adjoint method.

        For an objective function :math:`g(x, p)` where :math:`x` is the state
        vector of the system and :math:`p` is a vector of parameters, this
        computes the vector of sensitivities :math:`dg/dp`. This assumes that
        the system of equations has already been solved to find :math:`x`.

        :param perturb:
            A function with the signature ``perturb(sim, i, dp)`` which
            perturbs parameter ``i`` by a relative factor of ``dp``. To
            perturb a reaction rate constant, this function could be defined
            as::

                def perturb(sim, i, dp):
                    sim.gas.set_multiplier(1+dp, i)

            Calling ``perturb(sim, i, 0)`` should restore that parameter to its
            default value.
        :param n_params:
            The length of the vector of sensitivity parameters
        :param dgdx:
            The vector of partial derivatives of the function :math:`g(x, p)`
            with respect to the system state :math:`x`.
        :param g:
            A function with the signature ``value = g(sim)`` which computes the
            value of :math:`g(x,p)` at the current system state. This is used to
            compute :math:`\partial g/\partial p`. If this is identically zero
            (that is, :math:`g` is independent of :math:`p`) then this argument may
            be omitted.
        :param dp:
            A relative value by which to perturb each parameter
        """
        n_vars = self.sim.size()
        cdef np.ndarray[np.double_t, ndim=1] L = np.empty(n_vars)
        cdef np.ndarray[np.double_t, ndim=1] gg = \
                np.ascontiguousarray(dgdx, dtype=np.double)

        self.sim.solveAdjoint(&gg[0], &L[0])

        cdef np.ndarray[np.double_t, ndim=1] dgdp = np.empty(n_params)
        cdef np.ndarray[np.double_t, ndim=2] dfdp = np.empty((n_vars, n_params))
        cdef np.ndarray[np.double_t, ndim=1] fplus = np.empty(n_vars)
        cdef np.ndarray[np.double_t, ndim=1] fminus = np.empty(n_vars)
        gplus = gminus = 0

        for i in range(n_params):
            perturb(self, i, dp)
            if g:
                gplus = g(self)
            self.sim.getResidual(0, &fplus[0])

            perturb(self, i, -dp)
            if g:
                gminus = g(self)
            self.sim.getResidual(0, &fminus[0])

            perturb(self, i, 0)
            dgdp[i] = (gplus - gminus)/(2*dp)
            dfdp[:,i] = (fplus - fminus) / (2*dp)

        return dgdp - np.dot(L, dfdp)

    property grid_size_stats:
        """Return total grid size in each call to solve()"""
        def __get__(self):
            return self.sim.gridSizeStats()

    property jacobian_time_stats:
        """Return CPU time spent evaluating Jacobians in each call to solve()"""
        def __get__(self):
            return self.sim.jacobianTimeStats()

    property jacobian_count_stats:
        """Return number of Jacobian evaluations made in each call to solve()"""
        def __get__(self):
            return self.sim.jacobianCountStats()

    property eval_time_stats:
        """
        Return CPU time spent on non-Jacobian function evaluations in each call
        to solve()
        """
        def __get__(self):
            return self.sim.evalTimeStats()

    property eval_count_stats:
        """
        Return number of non-Jacobian function evaluations made in each call to
        solve()
        """
        def __get__(self):
            return self.sim.evalCountStats()

    property time_step_stats:
        """Return number of time steps taken in each call to solve()"""
        def __get__(self):
            return self.sim.timeStepStats()

    def set_max_grid_points(self, domain, npmax):
        """ Set the maximum number of grid points in the specified domain. """
        idom = self.domain_index(domain)
        self.sim.setMaxGridPoints(idom, npmax)

    def get_max_grid_points(self, domain):
        """ Get the maximum number of grid points in the specified domain. """
        idom = self.domain_index(domain)
        return self.sim.maxGridPoints(idom)
