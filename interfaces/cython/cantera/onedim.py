# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from math import erf as _erf
from pathlib import Path as _Path
import warnings
import numpy as np

from ._onedim import (
    AxisymmetricFlow, FreeFlow, Inlet1D, Outlet1D, ReactingSurface1D, Sim1D,
    Surface1D, SymmetryPlane1D, UnstrainedFlow,
)
from ._utils import CanteraError, __git_commit__, __version__, hdf_support
from .composite import Solution, SolutionArray


class FlameBase(Sim1D):
    """ Base class for flames with a single flow domain """
    __slots__ = ('gas',)

    def __init__(self, domains, gas, grid=None):
        """
        :param gas:
            object to use to evaluate all gas properties and reaction rates
        :param grid:
            array of initial grid points
        """
        if grid is None:
            grid = np.linspace(0.0, 0.1, 6)
        self.flame.grid = grid
        super().__init__(domains)

        #: The `Solution` object representing the species and reactions in the flame
        self.gas = gas
        self.flame.P = gas.P

    def set_refine_criteria(self, ratio=10.0, slope=0.8, curve=0.8, prune=0.0):
        """
        Set the criteria used for grid refinement.

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

        >>> f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0)
        """
        super().set_refine_criteria(self.flame, ratio, slope, curve, prune)

    def get_refine_criteria(self):
        """
        Get a dictionary of the criteria used for grid refinement. The items in
        the dictionary are the ``ratio``, ``slope``, ``curve``, and ``prune``,
        as defined in `~FlameBase.set_refine_criteria`.

        >>> f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0)
        >>> f.get_refine_criteria()
        {'ratio': 3.0, 'slope': 0.1, 'curve': 0.2, 'prune': 0.0}
        """
        return super().get_refine_criteria(self.flame)

    def set_initial_guess(self, *args, data=None, group=None, **kwargs):
        """
        Set the initial guess for the solution, and load restart data if
        provided. Derived classes extend this function to set approximations
        for the temperature and composition profiles.

        :param data:
            Restart data, which are typically based on an earlier simulation
            result. Restart data may be specified using a `SolutionArray`,
            `pandas.DataFrame`, or previously saved CSV, YAML or HDF container files.
            Note that restart data do not overwrite boundary conditions.
            HDF input requires Cantera compiled with HDF support. DataFrame input
            requires a working installation of *pandas*, which can be installed using
            pip or conda (package: ``pandas``).
        :param group:
            Group identifier within a HDF container file (only used in
            combination with HDF restart data).
        """
        super().set_initial_guess(*args, data=data, group=group, **kwargs)
        if data is None:
            return

        # load restart data into SolutionArray
        if isinstance(data, SolutionArray):
            # already a solution array
            arr = data

        elif isinstance(data, (str, _Path)):
            data = str(data)
            arr = SolutionArray(self.gas)
            if any(data.endswith(suffix) for suffix in [".hdf5", ".h5", ".hdf"]):
                # data source identifies a HDF file
                if "native" in hdf_support():
                    arr.restore(data, name=group, sub=self.domains[1].name)
                else:
                    arr.read_hdf(data, group=group, subgroup=self.domains[1].name)
            elif data.endswith(".yaml") or data.endswith(".yml"):
                # data source identifies a YAML file
                arr.restore(data, name=group, sub=self.domains[1].name)
            elif data.endswith('.csv'):
                # data source identifies a CSV file
                arr.read_csv(data)
            else:
                raise ValueError(f"'{data}' does not identify CSV, YAML or HDF file.")
        else:
            # data source is presumably a pandas DataFrame
            arr = SolutionArray(self.gas)
            arr.from_pandas(data)

        # get left and right boundaries
        left = self.domains[0]
        right = self.domains[2]

        if isinstance(left, Inlet1D) and isinstance(right, Inlet1D):
            # find stagnation plane
            i = np.flatnonzero(arr.velocity > 0)[-1]

            # adjust temperatures
            grid = arr.grid
            xi = (grid - grid[0]) / (grid[-1] - grid[0])
            T = arr.T
            T += (left.T - T[0]) * (1 - xi) + (right.T - T[-1]) * xi
            arr.TP = T, self.P

            # adjust velocities
            u = arr.velocity

            self.gas.TPY = left.T, self.P, left.Y
            arr[:i].velocity = u[:i] * left.mdot / self.gas.density / u[0]

            self.gas.TPY = right.T, self.P, right.Y
            arr[i:].velocity = - u[i:] * right.mdot / self.gas.density / u[-1]

        elif isinstance(left, Inlet1D):
            # adjust temperatures
            T = arr.T
            arr.TP = T + left.T - T[0], self.P

            # adjust velocities
            if not self.flame.domain_type.startswith("free"):
                self.gas.TPY = left.T, self.P, left.Y
                u0 = left.mdot / self.gas.density
                arr.velocity = u0 * arr.velocity / arr.velocity[0]

        self.from_array(arr)

    @property
    def max_grid_points(self):
        """
        Get/Set the maximum number of grid points used in the solution of
        this flame.
        """
        return super().get_max_grid_points(self.flame)

    @max_grid_points.setter
    def max_grid_points(self, npmax):
        super().set_max_grid_points(self.flame, npmax)

    @property
    def transport_model(self):
        """
        Get/Set the transport model used by the `Solution` object used for this
        simulation.
        """
        return self.flame.transport_model

    @transport_model.setter
    def transport_model(self, model):
        self.flame.transport_model = model

    @property
    def energy_enabled(self):
        """ Get/Set whether or not to solve the energy equation."""
        return self.flame.energy_enabled

    @energy_enabled.setter
    def energy_enabled(self, enable):
        self.flame.energy_enabled = enable

    @property
    def soret_enabled(self):
        """
        Get/Set whether or not to include diffusive mass fluxes due to the
        Soret effect. Enabling this option works only when using the
        multicomponent transport model.
        """
        return self.flame.soret_enabled

    @soret_enabled.setter
    def soret_enabled(self, enable):
        self.flame.soret_enabled = enable

    @property
    def flux_gradient_basis(self):
        """
        Get/Set whether or not species diffusive fluxes are computed with
        respect to their mass fraction gradients ('mass')
        or mole fraction gradients ('molar', default) when
        using the mixture-averaged transport model.
        """
        return self.flame.flux_gradient_basis

    @flux_gradient_basis.setter
    def flux_gradient_basis(self, basis):
        self.flame.flux_gradient_basis = basis

    @property
    def radiation_enabled(self):
        """
        Get/Set whether or not to include radiative heat transfer
        """
        return self.flame.radiation_enabled

    @radiation_enabled.setter
    def radiation_enabled(self, enable):
        self.flame.radiation_enabled = enable

    @property
    def boundary_emissivities(self):
        """ Set/get boundary emissivities. """
        return self.flame.boundary_emissivities

    @boundary_emissivities.setter
    def boundary_emissivities(self, epsilon):
        if len(epsilon) != 2:
            raise ValueError("Boundary emissivities must both be set at the same time.")
        self.flame.boundary_emissivities = epsilon[0], epsilon[1]

    @property
    def grid(self):
        """ Array of grid point positions along the flame. """
        return self.flame.grid

    @property
    def P(self):
        """ Get/Set the pressure of the flame [Pa] """
        return self.flame.P

    @P.setter
    def P(self, P):
        self.flame.P = P

    @property
    def T(self):
        """ Array containing the temperature [K] at each grid point. """
        return self.flame.values("T")

    @property
    def velocity(self):
        """
        Array containing the velocity [m/s] normal to the flame at each point.
        """
        return self.flame.values("velocity")

    @property
    def spread_rate(self):
        """
        Array containing the tangential velocity gradient [1/s] (that is, radial
        velocity divided by radius) at each point.
        """
        return self.flame.values("spreadRate")

    @property
    def L(self):
        """
        Array containing the radial pressure gradient (1/r)(dP/dr) [N/m⁴] at
        each point. Note: This value is named ``Lambda`` in the C++ code.
        """
        return self.flame.values("Lambda")

    @property
    def E(self):
        """
        Array containing the electric field strength at each point.
        Note: This value is named ``eField`` in the C++ code and is only defined if
        the transport model is ``ionized-gas``.
        """
        return self.flame.values("eField")

    @property
    def Uo(self):
        """
        Array containing the oxidizer velocity (right boundary velocity) [m/s] at
        each point.
        Note: This value is named ``Uo`` in the C++ code and is only defined when using
        two-point control.
        """
        return self.flame.values("Uo")

    @property
    def left_control_point_temperature(self):
        """ Get/Set the left control point temperature [K] """
        return self.flame.left_control_point_temperature

    @left_control_point_temperature.setter
    def left_control_point_temperature(self, T):
        self.flame.left_control_point_temperature = T

    @property
    def left_control_point_coordinate(self):
        """ Get the left control point coordinate [m] """
        return self.flame.left_control_point_coordinate

    @property
    def right_control_point_temperature(self):
        """ Get/Set the right control point temperature [K] """
        return self.flame.right_control_point_temperature

    @right_control_point_temperature.setter
    def right_control_point_temperature(self, T):
        self.flame.right_control_point_temperature = T

    @property
    def right_control_point_coordinate(self):
        """ Get the right control point coordinate [m] """
        return self.flame.right_control_point_coordinate

    def elemental_mass_fraction(self, m):
        r"""
        Get the elemental mass fraction :math:`Z_{\mathrm{mass},m}` of element
        :math:`m` at each grid point, which is defined as:

        .. math:: Z_{\mathrm{mass},m} = \sum_k \frac{a_{m,k} M_m}{M_k} Y_k

        with :math:`a_{m,k}` being the number of atoms of element :math:`m` in
        species :math:`k`, :math:`M_m` the atomic weight of element :math:`m`,
        :math:`M_k` the molecular weight of species :math:`k`, and :math:`Y_k`
        the mass fraction of species :math:`k`.

        :param m:
            Base element, may be specified by name or by index.

        >>> phase.elemental_mass_fraction('H')
        [1.0, ..., 0.0]
        """
        vals = np.empty(self.flame.n_points)
        for i in range(self.flame.n_points):
            self.flame.update_state(i)
            vals[i] = self.gas.elemental_mass_fraction(m)
        return vals

    def elemental_mole_fraction(self, m):
        r"""
        Get the elemental mole fraction :math:`Z_{\mathrm{mole},m}` of element
        :math:`m` at each grid point, which is defined as:

        .. math:: Z_{\mathrm{mole},m} = \sum_k \frac{a_{m,k}}{\sum_j a_{j,k}} X_k

        with :math:`a_{m,k}` being the number of atoms of element :math:`m` in
        species :math:`k` and :math:`X_k` the mole fraction of species
        :math:`k`.

        :param m:
            Base element, may be specified by name or by index.

        >>> phase.elemental_mole_fraction('H')
        [1.0, ..., 0.0]
        """
        vals = np.empty(self.flame.n_points)
        for i in range(self.flame.n_points):
            self.flame.update_state(i)
            vals[i] = self.gas.elemental_mole_fraction(m)
        return vals

    def to_array(self, domain=None, normalize=False):
        """
        Retrieve domain data as a `SolutionArray` object.

        :param domain:
            Domain to be converted; by default, the method retrieves the flow domain
        :param normalize:
            Boolean flag indicating whether mass/mole fractions should be normalized

        .. versionadded:: 3.0
        """
        if domain is None:
            domain = self.flame
        else:
            domain = self.domains[self.domain_index(domain)]
        dest = SolutionArray(domain.phase, init=False)
        dest = domain._to_array(dest, normalize)
        dest.shape = dest._api_shape()
        return dest

    def from_array(self, arr, domain=None):
        """
        Restore the solution vector from a `SolutionArray` object.

        :param arr:
            `SolutionArray` containing data to be restored.
        :param domain:
            Domain to be converted; by default, the method retrieves the flow domain

        .. versionadded:: 3.0
        """
        if domain is None:
            domain = self.flame
        else:
            domain = self.domains[self.domain_index(domain)]
        domain._from_array(arr)

    def to_pandas(self, species='X', normalize=True):
        """
        Return the solution vector as a `pandas.DataFrame`.

        :param species:
            Attribute to use obtaining species profiles, for example ``X`` for
            mole fractions or ``Y`` for mass fractions.
        :param normalize:
            Boolean flag to indicate whether the mole/mass fractions should
            be normalized (default is ``True``)

        This method uses `to_array` and requires a working *pandas*
        installation. Use pip or conda to install ``pandas`` to enable this
        method.
        """
        cols = ('extra', 'T', 'D', species)
        return self.to_array(normalize=normalize).to_pandas(cols=cols)

    @property
    def electric_field_enabled(self):
        """
        Get/Set whether or not to solve the Poisson's equation.

        Used in conjunction with the ``ion-transport`` transport model to model
        diffusion of ionized species. See classes :ct:`IonFlow` and
        :ct:`IonGasTransport`.
        """
        return self.flame.electric_field_enabled

    @electric_field_enabled.setter
    def electric_field_enabled(self, enable):
        self.flame.electric_field_enabled = enable

    @property
    def two_point_control_enabled(self):
        """
        Get/Set whether or not to active two point flame control.
        """
        return self.flame.two_point_control_enabled

    @two_point_control_enabled.setter
    def two_point_control_enabled(self, enable):
        self.flame.two_point_control_enabled = enable


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
            self.flame.update_state(i)
            vals[...,i] = getattr(self.gas, attr)

        return vals

    if size is None:
        extradoc = "\nReturns an array of length `n_points`."
    else:
        extradoc = "\nReturns an array of size `%s` x `n_points`." % size

    basedoc = getattr(Solution, attr).__doc__


    doc = _trim(getattr(Solution, attr).__doc__) +'\n' + extradoc
    return property(getter, doc=doc)

# Add scalar properties to FlameBase
for _attr in ['density', 'density_mass', 'density_mole', 'volume_mass',
              'volume_mole', 'int_energy_mole', 'int_energy_mass', 'h',
              'enthalpy_mole', 'enthalpy_mass', 's', 'entropy_mole',
              'entropy_mass', 'g', 'gibbs_mole', 'gibbs_mass', 'cv',
              'cv_mole', 'cv_mass', 'cp', 'cp_mole', 'cp_mass',
              'isothermal_compressibility', 'thermal_expansion_coeff',
              'sound_speed', 'viscosity', 'thermal_conductivity',
              'heat_release_rate', 'mean_molecular_weight']:
    setattr(FlameBase, _attr, _array_property(_attr))
FlameBase.volume = _array_property('v') # avoid confusion with velocity gradient 'V'
FlameBase.int_energy = _array_property('u') # avoid collision with velocity 'u'

# Add properties with values for each species
for _attr in ['X', 'Y', 'concentrations', 'partial_molar_enthalpies',
              'partial_molar_entropies', 'partial_molar_int_energies',
              'chemical_potentials', 'electrochemical_potentials', 'partial_molar_cp',
              'partial_molar_volumes', 'standard_enthalpies_RT',
              'standard_entropies_R', 'standard_int_energies_RT',
              'standard_gibbs_RT', 'standard_cp_R', 'creation_rates',
              'destruction_rates', 'net_production_rates', 'creation_rates_ddC',
              'creation_rates_ddP', 'creation_rates_ddT', 'destruction_rates_ddC',
              'destruction_rates_ddP', 'destruction_rates_ddT',
              'net_production_rates_ddC', 'net_production_rates_ddP',
              'net_production_rates_ddT', 'mix_diff_coeffs', 'mix_diff_coeffs_mass',
              'mix_diff_coeffs_mole', 'thermal_diff_coeffs', 'activities',
              'activity_coefficients', 'mobilities', 'species_viscosities']:
    setattr(FlameBase, _attr, _array_property(_attr, 'n_species'))

# Remove misleading examples and references to setters that don't exist
FlameBase.X.__doc__ = "Array of mole fractions of size `n_species` x `n_points`"
FlameBase.Y.__doc__ = "Array of mass fractions of size `n_species` x `n_points`"
FlameBase.concentrations.__doc__ = ("Array of species concentrations [kmol/m³]"
                                    " of size `n_species` x `n_points`")

# Add properties with values for each reaction
for _attr in ['forward_rates_of_progress', 'reverse_rates_of_progress', 'net_rates_of_progress',
              'equilibrium_constants', 'forward_rate_constants', 'reverse_rate_constants',
              'delta_enthalpy', 'delta_gibbs', 'delta_entropy',
              'delta_standard_enthalpy', 'delta_standard_gibbs',
              'delta_standard_entropy', 'heat_production_rates',
              'third_body_concentrations', 'forward_rate_constants_ddC',
              'forward_rate_constants_ddP', 'forward_rate_constants_ddT',
              'forward_rates_of_progress_ddC', 'forward_rates_of_progress_ddP',
              'forward_rates_of_progress_ddT', 'net_rates_of_progress_ddC',
              'net_rates_of_progress_ddP', 'net_rates_of_progress_ddT',
              'reverse_rates_of_progress_ddC', 'reverse_rates_of_progress_ddP',
              'reverse_rates_of_progress_ddT']:
    setattr(FlameBase, _attr, _array_property(_attr, 'n_reactions'))


class FreeFlame(FlameBase):
    """A freely-propagating flat flame."""
    __slots__ = ('inlet', 'flame', 'outlet')

    def __init__(self, gas, grid=None, width=None):
        """
        A domain of type `FreeFlow` named 'flame' will be created to represent
        the flame. The three domains comprising the stack are stored as ``self.inlet``,
        ``self.flame``, and ``self.outlet``.

        :param grid:
            A list of points to be used as the initial grid. Not recommended
            unless solving only on a fixed grid; Use the `width` parameter
            instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.
        """

        #: `Inlet1D` at the left of the domain representing premixed reactants
        self.inlet = Inlet1D(name='reactants', phase=gas)

        #: `Outlet1D` at the right of the domain representing the burned products
        self.outlet = Outlet1D(name='products', phase=gas)

        if not hasattr(self, 'flame'):
            # Create flame domain if not already instantiated by a child class
            #: `FreeFlow` domain representing the flame
            self.flame = FreeFlow(gas, name='flame')

        if width is not None:
            if grid is not None:
                raise ValueError("'grid' and 'width' arguments are mutually exclusive")
            grid = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]) * width

        super().__init__((self.inlet, self.flame, self.outlet), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.inlet.T = gas.T
        self.inlet.X = gas.X

    def set_initial_guess(self, locs=[0.0, 0.3, 0.5, 1.0], data=None, group=None):
        """
        Set the initial guess for the solution. By default, the adiabatic flame
        temperature and equilibrium composition are computed for the inlet gas
        composition. Alternatively, a previously calculated result can be
        supplied as an initial guess via 'data' and 'key' inputs (see
        `FlameBase.set_initial_guess`).

        :param locs:
            A list of four locations to define the temperature and mass fraction
            profiles. Profiles rise linearly between the second and third
            location. Locations are given as a fraction of the entire domain
        """
        super().set_initial_guess(data=data, group=group)
        if data is not None:
            # set fixed temperature
            Tmid = .75 * self.T[0] + .25 * self.T[-1]
            i = np.flatnonzero(self.T < Tmid)[-1]
            self.fixed_temperature = self.T[i]

            return

        self.gas.TPY = self.inlet.T, self.P, self.inlet.Y

        if not self.inlet.mdot:
            # nonzero initial guess increases likelihood of convergence
            self.inlet.mdot = 1.0 * self.gas.density

        Y0 = self.inlet.Y
        u0 = self.inlet.mdot / self.gas.density
        T0 = self.inlet.T

        # get adiabatic flame temperature and composition
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y
        u1 = self.inlet.mdot / self.gas.density

        self.flame.set_profile('velocity', locs, [u0, u0, u1, u1])
        self.flame.set_profile('T', locs, [T0, T0, Teq, Teq])

        # Pick the location of the fixed temperature point, using an existing
        # point if a reasonable choice exists
        T = self.T
        Tmid = 0.75 * T0 + 0.25 * Teq
        i = np.flatnonzero(T < Tmid)[-1] # last point less than Tmid
        if Tmid - T[i] < 0.2 * (Tmid - T0):
            self.fixed_temperature = T[i]
        elif T[i+1] - Tmid < 0.2 * (Teq - Tmid):
            self.fixed_temperature = T[i+1]
        else:
            self.fixed_temperature = Tmid

        for n in range(self.gas.n_species):
            self.flame.set_profile(self.gas.species_name(n),
                             locs, [Y0[n], Y0[n], Yeq[n], Yeq[n]])

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
            return super().solve(loglevel, refine_grid, auto)

        # Use a callback function to check that the domain is actually wide
        # enough to contain the flame after each steady-state solve. If the user
        # provided a callback, store this so it can called in addition to our
        # callback, and restored at the end.
        original_callback = self._steady_callback

        class DomainTooNarrow(Exception): pass

        def check_width(t):
            T = self.T
            x = self.grid
            mRef = (T[-1] - T[0]) / (x[-1] - x[0])
            mLeft = (T[1] - T[0]) / (x[1] - x[0]) / mRef
            mRight = (T[-3] - T[-1]) / (x[-3] - x[-1]) / mRef

            # The domain is considered too narrow if gradient at the left or
            # right edge is significant, compared to the average gradient across
            # the domain.
            if mLeft > 0.02 or mRight > 0.02:
                raise DomainTooNarrow()

            if original_callback:
                return original_callback(t)
            else:
                return 0.0

        self.set_steady_callback(check_width)

        for _ in range(12):
            try:
                super().solve(loglevel, refine_grid, auto)
                break
            except DomainTooNarrow:
                self.flame.grid *= 2
                if loglevel > 0:
                    print('Expanding domain to accommodate flame thickness. '
                          'New width: {} m'.format(
                          self.flame.grid[-1] - self.flame.grid[0]))
                if refine_grid:
                    self.refine(loglevel)

        self.set_steady_callback(original_callback)

    def get_flame_speed_reaction_sensitivities(self):
        r"""
        Compute the normalized sensitivities of the laminar flame speed
        :math:`S_u` with respect to the reaction rate constants :math:`k_i`:

        .. math::

            s_i = \frac{k_i}{S_u} \frac{dS_u}{dk_i}
        """

        def g(sim):
            return sim.velocity[0]

        Nvars = sum(D.n_components * D.n_points for D in self.domains)

        # Index of u[0] in the global solution vector
        i_Su = self.inlet.n_components + self.flame.component_index('velocity')

        dgdx = np.zeros(Nvars)
        dgdx[i_Su] = 1

        Su0 = g(self)

        def perturb(sim, i, dp):
            sim.gas.set_multiplier(1+dp, i)

        return self.solve_adjoint(perturb, self.gas.n_reactions, dgdx) / Su0


class BurnerFlame(FlameBase):
    """A burner-stabilized flat flame."""
    __slots__ = ('burner', 'flame', 'outlet')

    def __init__(self, gas, grid=None, width=None):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            A list of points to be used as the initial grid. Not recommended
            unless solving only on a fixed grid; Use the `width` parameter
            instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.

        A domain of class `UnstrainedFlow` named ``flame`` will be created to
        represent the flame. The three domains comprising the stack are stored as
        ``self.burner``, ``self.flame``, and ``self.outlet``.
        """
        #: `Inlet1D` at the left of the domain representing the burner surface through
        #: which reactants flow
        self.burner = Inlet1D(name='burner', phase=gas)

        #: `Outlet1D` at the right of the domain representing the burned gas
        self.outlet = Outlet1D(name='outlet', phase=gas)

        if not hasattr(self, 'flame'):
            # Create flame domain if not already instantiated by a child class
            #: `UnstrainedFlow` domain representing the flame
            self.flame = UnstrainedFlow(gas, name='flame')

        if width is not None:
            if grid is not None:
                raise ValueError("'grid' and 'width' arguments are mutually exclusive")
            grid = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]) * width

        super().__init__((self.burner, self.flame, self.outlet), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.burner.T = gas.T
        self.burner.X = gas.X

    def set_initial_guess(self, data=None, group=None):
        """
        Set the initial guess for the solution. By default, the adiabatic flame
        temperature and equilibrium composition are computed for the burner
        gas composition. The temperature profile rises linearly in the first
        20% of the flame to Tad, then is flat. The mass fraction profiles are
        set similarly. Alternatively, a previously calculated result can be
        supplied as an initial guess  via 'data' and 'key' inputs (see
        `FlameBase.set_initial_guess`).
        """
        super().set_initial_guess(data=data, group=group)
        if data:
            return

        self.gas.TPY = self.burner.T, self.P, self.burner.Y
        Y0 = self.burner.Y
        u0 = self.burner.mdot / self.gas.density
        T0 = self.burner.T

        # get adiabatic flame temperature and composition
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y
        u1 = self.burner.mdot / self.gas.density

        locs = [0.0, 0.2, 1.0]
        self.flame.set_profile('velocity', locs, [u0, u1, u1])
        self.flame.set_profile('T', locs, [T0, Teq, Teq])
        for n in range(self.gas.n_species):
            self.flame.set_profile(self.gas.species_name(n),
                             locs, [Y0[n], Yeq[n], Yeq[n]])

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
        # Use a callback function to check that the flame has not been blown off
        # the burner surface. If the user provided a callback, store this so it
        # can called in addition to our callback, and restored at the end.
        original_callback = self._steady_callback

        class FlameBlowoff(Exception): pass

        if auto:
            def check_blowoff(t):
                T = self.T
                n = max(3, len(self.T) // 5)

                # Near-zero temperature gradient at burner indicates blowoff
                if abs(T[n] - T[0]) / (T[-1] - T[0]) < 1e-6:
                    raise FlameBlowoff()

                if original_callback:
                    return original_callback(t)
                else:
                    return 0.0

            self.set_steady_callback(check_blowoff)

        try:
            return super().solve(loglevel, refine_grid, auto)
        except FlameBlowoff:
            # The eventual solution for a blown off flame is the non-reacting
            # solution, so just set the state to this now
            self.flame.set_flat_profile("T", self.T[0])
            for k,spec in enumerate(self.gas.species_names):
                self.flame.set_flat_profile(spec, self.burner.Y[k])

            self.set_steady_callback(original_callback)
            super().solve(loglevel, False, False)
            if loglevel > 0:
                print('Flame has blown off of burner (non-reacting solution)')

        self.set_steady_callback(original_callback)


class CounterflowDiffusionFlame(FlameBase):
    """ A counterflow diffusion flame """
    __slots__ = ('fuel_inlet', 'flame', 'oxidizer_inlet')

    def __init__(self, gas, grid=None, width=None):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            A list of points to be used as the initial grid. Not recommended
            unless solving only on a fixed grid; Use the `width` parameter
            instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.

        A domain of class `AxisymmetricFlow` named ``flame`` will be created to
        represent the flame. The three domains comprising the stack are stored as
        ``self.fuel_inlet``, ``self.flame``, and ``self.oxidizer_inlet``.
        """

        #: `Inlet1D` at the left of the domain representing the fuel mixture
        self.fuel_inlet = Inlet1D(name='fuel_inlet', phase=gas)
        self.fuel_inlet.T = gas.T

        #: `Inlet1D` at the right of the domain representing the oxidizer mixture
        self.oxidizer_inlet = Inlet1D(name='oxidizer_inlet', phase=gas)
        self.oxidizer_inlet.T = gas.T

        #: `AxisymmetricFlow` domain representing the flame
        self.flame = AxisymmetricFlow(gas, name='flame')

        if width is not None:
            if grid is not None:
                raise ValueError("'grid' and 'width' arguments are mutually exclusive")
            grid = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) * width

        super().__init__((self.fuel_inlet, self.flame, self.oxidizer_inlet), gas, grid)

    def set_initial_guess(self, data=None, group=None):
        """
        Set the initial guess for the solution. By default, the initial guess
        is generated by assuming infinitely-fast chemistry. Alternatively, a
        previously calculated result can be supplied as an initial guess via
        'data' and 'key' inputs (see `FlameBase.set_initial_guess`).
        """
        super().set_initial_guess(data=data, group=group)
        if data:
            return

        # Compute stoichiometric mixture composition
        Yin_f = self.fuel_inlet.Y
        self.gas.TPY = self.fuel_inlet.T, self.P, Yin_f
        mdotf = self.fuel_inlet.mdot
        rho0f = self.gas.density
        u0f = mdotf / rho0f
        T0f = self.fuel_inlet.T

        Yin_o = self.oxidizer_inlet.Y
        self.gas.TPY = self.oxidizer_inlet.T, self.P, Yin_o
        mdoto = self.oxidizer_inlet.mdot
        rho0o = self.gas.density
        u0o = mdoto / rho0o
        T0o = self.oxidizer_inlet.T

        if mdoto == mdotf == 0.0:
            raise CanteraError("Mass flow for fuel and/or oxidizer "
                               "must be positive")

        zst = 1 / (1 + self.gas.stoich_air_fuel_ratio(Yin_f, Yin_o, 'mass'))
        Yst = zst * Yin_f + (1.0 - zst) * Yin_o

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
        kOx = (self.gas.species_index('O2') if 'O2' in self.gas.species_names else
               self.gas.species_index('o2'))
        f = np.sqrt(a / (2.0 * self.gas.mix_diff_coeffs[kOx]))
        L = - 0.5 * (rho0o + rho0f) * a**2

        x0 = np.sqrt(mdotf*u0f) * dz / (np.sqrt(mdotf*u0f) + np.sqrt(mdoto*u0o))
        nz = len(zz)

        Y = np.zeros((nz, self.gas.n_species))
        T = np.zeros(nz)
        for j in range(nz):
            x = zz[j] - zz[0]
            zeta = f * (x - x0)
            zmix = 0.5 * (1.0 - _erf(zeta))
            if zmix > zst:
                Y[j] = Yeq + (Yin_f - Yeq) * (zmix - zst) / (1.0 - zst)
                T[j] = Teq + (T0f - Teq) * (zmix - zst) / (1.0 - zst)
            else:
                Y[j] = Yin_o + zmix * (Yeq - Yin_o) / zst
                T[j] = T0o + (Teq - T0o) * zmix / zst

        T[0] = T0f
        T[-1] = T0o
        zrel = (zz - zz[0])/dz

        self.flame.set_profile('velocity', [0.0, 1.0], [u0f, -u0o])
        self.flame.set_profile('spreadRate', [0.0, x0/dz, 1.0], [0.0, a, 0.0])
        self.flame.set_profile("Lambda", [0.0, 1.0], [L, L])
        self.flame.set_profile('T', zrel, T)
        for k,spec in enumerate(self.gas.species_names):
            self.flame.set_profile(spec, zrel, Y[:,k])

    def extinct(self):
        return max(self.T) - max(self.fuel_inlet.T, self.oxidizer_inlet.T) < 10

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
        if self.flame.transport_model == 'ionized-gas':
            warnings.warn(
                "The 'ionized-gas' transport model is untested for "
                "'CounterflowDiffusionFlame' objects.", UserWarning)

        super().solve(loglevel, refine_grid, auto)
        # Do some checks if loglevel is set
        if loglevel > 0:
            if self.extinct():
                print('WARNING: Flame is extinct.')
            else:
                # Check if the flame is very thick
                # crude width estimate based on temperature
                z_flame = self.grid[self.T > np.max(self.T) / 2]
                flame_width = z_flame[-1] - z_flame[0]
                domain_width = self.grid[-1] - self.grid[0]
                if flame_width / domain_width > 0.4:
                    print('WARNING: The flame is thick compared to the domain '
                          'size. The flame might be affected by the plug-flow '
                          'boundary conditions. Consider increasing the inlet mass '
                          'fluxes or using a larger domain.')

                # Check if the temperature peak is close to a boundary
                z_center = (self.grid[np.argmax(self.T)] - self.grid[0]) / domain_width
                if z_center < 0.25:
                    print('WARNING: The flame temperature peak is close to the '
                          'fuel inlet. Consider increasing the ratio of the '
                          'fuel inlet mass flux to the oxidizer inlet mass flux.')
                if z_center > 0.75:
                    print('WARNING: The flame temperature peak is close to the '
                          'oxidizer inlet. Consider increasing the ratio of the '
                          'oxidizer inlet mass flux to the fuel inlet mass flux.')

    def strain_rate(self, definition, fuel=None, oxidizer='O2', stoich=None):
        r"""
        Return the axial strain rate of the counterflow diffusion flame in 1/s.

        :param definition:
            The definition of the strain rate to be calculated. Options are:
            ``mean``, ``max``, ``stoichiometric``, ``potential_flow_fuel``, and
            ``potential_flow_oxidizer``.
        :param fuel: The fuel species. Used only if ``definition`` is
            ``stoichiometric``.
        :param oxidizer: The oxidizer species, default ``O2``. Used only if
            ``definition`` is ``stoichiometric``.
        :param stoich: The molar stoichiometric oxidizer-to-fuel ratio.
            Can be omitted if the oxidizer is ``O2``. Used only if ``definition``
            is ``stoichiometric``.

        The parameter ``definition`` sets the method to compute the strain rate.
        Possible options are:

        ``mean``:
            The mean axial velocity gradient in the entire domain

           .. math:: a_{mean} = \left| \frac{\Delta u}{\Delta z} \right|

        ``max``:
            The maximum axial velocity gradient

            .. math:: a_{max} = \max \left( \left| \frac{du}{dz} \right| \right)

        ``stoichiometric``:
            The axial velocity gradient at the stoichiometric surface.

            .. math::

                a_{stoichiometric} = \left| \left. \frac{du}{dz}
                \right|_{\phi=1} \right|

            This method uses the additional keyword arguments ``fuel``,
            ``oxidizer``, and ``stoich``.

            >>> f.strain_rate('stoichiometric', fuel='H2', oxidizer='O2',
                              stoich=0.5)

        ``potential_flow_fuel``:
            The corresponding axial strain rate for a potential flow boundary
            condition at the fuel inlet.

            .. math:: a_{f} = \sqrt{-\frac{\Lambda}{\rho_{f}}}

        ``potential_flow_oxidizer``:
            The corresponding axial strain rate for a potential flow boundary
            condition at the oxidizer inlet.

            .. math:: a_{o} = \sqrt{-\frac{\Lambda}{\rho_{o}}}
        """
        if definition == 'mean':
            return - (self.velocity[-1] - self.velocity[0]) / self.grid[-1]

        elif definition == 'max':
            return np.max(np.abs(np.gradient(self.velocity) / np.gradient(self.grid)))

        elif definition == 'stoichiometric':
            if fuel is None:
                raise KeyError('Required argument "fuel" not defined')
            if oxidizer != 'O2' and stoich is None:
                raise KeyError('Required argument "stoich" not defined')

            if stoich is None:
                # oxidizer is O2
                stoich = - 0.5 * self.gas.n_atoms(fuel, 'O')
                if 'H' in self.gas.element_names:
                    stoich += 0.25 * self.gas.n_atoms(fuel, 'H')
                if 'C' in self.gas.element_names:
                    stoich += self.gas.n_atoms(fuel, 'C')

            d_u_d_z = np.gradient(self.velocity) / np.gradient(self.grid)
            phi = (self.X[self.gas.species_index(fuel)] * stoich /
                   np.maximum(self.X[self.gas.species_index(oxidizer)], 1e-20))
            z_stoich = np.interp(-1., -phi, self.grid)
            return np.abs(np.interp(z_stoich, self.grid, d_u_d_z))

        elif definition == 'potential_flow_fuel':
            return np.sqrt(- self.flame.radial_pressure_gradient[0] / self.density[0])

        elif definition == 'potential_flow_oxidizer':
            return np.sqrt(- self.flame.radial_pressure_gradient[0] / self.density[-1])

        else:
            raise ValueError('Definition "' + definition + '" is not available')

    def mixture_fraction(self, m):
        r"""
        Compute the mixture fraction based on element ``m`` or from the
        Bilger mixture fraction by setting ``m="Bilger"``

        The mixture fraction is computed from the elemental mass fraction of
        element ``m``, normalized by its values on the fuel and oxidizer
        inlets:

        .. math:: Z = \frac{Z_{\mathrm{mass},m}(z) -
                            Z_{\mathrm{mass},m}(z_\mathrm{oxidizer})}
                           {Z_{\mathrm{mass},m}(z_\mathrm{fuel}) -
                            Z_{\mathrm{mass},m}(z_\mathrm{oxidizer})}

        or from the Bilger mixture fraction:

        .. math:: Z = \frac{\beta-\beta_{\mathrm{oxidizer}}}
                           {\beta_{\mathrm{fuel}}-\beta_{\mathrm{oxidizer}}}

        with

        .. math:: \beta = 2\frac{Z_C}{M_C}+2\frac{Z_S}{M_S}
                          +\frac{1}{2}\frac{Z_H}{M_H}-\frac{Z_O}{M_O}

        :param m:
            The element based on which the mixture fraction is computed,
            may be specified by name or by index, or "Bilger" for the
            Bilger mixture fraction, which considers the elements C,
            H, S, and O

        >>> f.mixture_fraction('H')
        >>> f.mixture_fraction('Bilger')
        """

        self.flame.update_state(0)
        Yf = self.gas.Y
        self.flame.update_state(self.flame.n_points - 1)
        Yo = self.gas.Y

        vals = np.empty(self.flame.n_points)
        for i in range(self.flame.n_points):
            self.flame.update_state(i)
            vals[i] = self.gas.mixture_fraction(Yf, Yo, 'mass', m)
        return vals

    @property
    def equivalence_ratio(self):
        self.flame.update_state(0)
        Yf = self.gas.Y
        self.flame.update_state(self.flame.n_points - 1)
        Yo = self.gas.Y

        vals = np.empty(self.flame.n_points)
        for i in range(self.flame.n_points):
            self.flame.update_state(i)
            vals[i] = self.gas.equivalence_ratio(Yf, Yo, "mass")
        return vals


class ImpingingJet(FlameBase):
    """An axisymmetric flow impinging on a surface at normal incidence."""
    __slots__ = ('inlet', 'flame', 'surface')

    def __init__(self, gas, grid=None, width=None, surface=None):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            A list of points to be used as the initial grid. Not recommended
            unless solving only on the initial grid; Use the `width` parameter
            instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.
        :param surface:
            A Kinetics object used to compute any surface reactions.

        A domain of class `AxisymmetricFlow` named ``flame`` will be created to
        represent the flame. The three domains comprising the stack are stored as
        ``self.inlet``, ``self.flame``, and ``self.surface``.
        """

        #: `Inlet1D` at the left of the domain representing the incoming reactants
        self.inlet = Inlet1D(name='inlet', phase=gas)

        #: `AxisymmetricFlow` domain representing the flame
        self.flame = AxisymmetricFlow(gas, name='flame')
        self.flame.set_axisymmetric_flow()

        if width is not None:
            if grid is not None:
                raise ValueError("'grid' and 'width' arguments are mutually exclusive")
            grid = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) * width

        if surface is None:
            #: `Surface1D` or `ReactingSurface1D` domain representing the surface the
            #: flow is impinging on
            self.surface = Surface1D(name='surface', phase=gas)
            self.surface.T = gas.T
        else:
            self.surface = ReactingSurface1D(name='surface', phase=surface)
            self.surface.T = surface.T

        super().__init__((self.inlet, self.flame, self.surface), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.inlet.T = gas.T
        self.inlet.X = gas.X

    def set_initial_guess(self, products='inlet', data=None, group=None):
        """
        Set the initial guess for the solution. If products = 'equil', then
        the equilibrium composition at the adiabatic flame temperature will be
        used to form the initial guess. Otherwise the inlet composition will
        be used. Alternatively, a previously calculated result can be supplied
        as an initial guess via 'data' and 'key' inputs (see
        `FlameBase.set_initial_guess`).
        """
        super().set_initial_guess(data=data, group=group, products=products)
        if data:
            return

        Y0 = self.inlet.Y
        T0 = self.inlet.T
        self.gas.TPY = T0, self.flame.P, Y0
        u0 = self.inlet.mdot / self.gas.density

        if products == 'equil':
            self.gas.equilibrate('HP')
            Teq = self.gas.T
            Yeq = self.gas.Y
            locs = np.array([0.0, 0.3, 0.7, 1.0])
            self.flame.set_profile('T', locs, [T0, Teq, Teq, self.surface.T])
            for k in range(self.gas.n_species):
                self.flame.set_profile(self.gas.species_name(k), locs,
                                       [Y0[k], Yeq[k], Yeq[k], Yeq[k]])
        else:
            locs = np.array([0.0, 1.0])
            self.flame.set_profile('T', locs, [T0, self.surface.T])
            for k in range(self.gas.n_species):
                self.flame.set_profile(self.gas.species_name(k), locs, [Y0[k], Y0[k]])

        locs = np.array([0.0, 1.0])
        self.flame.set_profile("velocity", locs, [u0, 0.0])
        self.flame.set_profile("spreadRate", locs, [0.0, 0.0])


class CounterflowPremixedFlame(FlameBase):
    """ A premixed counterflow flame """
    __slots__ = ('reactants', 'flame', 'products')

    def __init__(self, gas, grid=None, width=None):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            Array of initial grid points. Not recommended unless solving only on
            a fixed grid; Use the `width` parameter instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.

        A domain of class `AxisymmetricFlow` named ``flame`` will be created to
        represent the flame. The three domains comprising the stack are stored as
        ``self.reactants``, ``self.flame``, and ``self.products``.
        """

        #: `Inlet1D` at the left of the domain representing premixed reactants
        self.reactants = Inlet1D(name='reactants', phase=gas)
        self.reactants.T = gas.T

        #: `Inlet1D` at the right of the domain representing burned products
        self.products = Inlet1D(name='products', phase=gas)
        self.products.T = gas.T

        #: `AxisymmetricFlow` domain representing the flame
        self.flame = AxisymmetricFlow(gas, name='flame')

        if width is not None:
            if grid is not None:
                raise ValueError("'grid' and 'width' arguments are mutually exclusive")
            # Create grid points aligned with initial guess profile
            grid = np.array([0.0, 0.3, 0.5, 0.7, 1.0]) * width

        super().__init__((self.reactants, self.flame, self.products), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.reactants.X = gas.X

    def set_initial_guess(self, equilibrate=True, data=None, group=None):
        """
        Set the initial guess for the solution.

        If `equilibrate` is True, then the products composition and temperature
        will be set to the equilibrium state of the reactants mixture.
        Alternatively, a previously calculated result can be supplied as an
        initial guess via 'data' and 'key' inputs (see
        `FlameBase.set_initial_guess`).
        """
        super().set_initial_guess(data=data, group=group, equilibrate=equilibrate)
        if data:
            return

        Yu = self.reactants.Y
        Tu = self.reactants.T
        self.gas.TPY = Tu, self.flame.P, Yu
        rhou = self.gas.density
        uu = self.reactants.mdot / rhou

        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y

        if equilibrate:
            Tb = Teq
            Yb = Yeq
            self.products.Y = Yb
            self.products.T = Tb
        else:
            Tb = self.products.T
            Yb = self.products.Y

        self.gas.TPY = Tb, self.flame.P, Yb
        rhob = self.gas.density
        ub = self.products.mdot / rhob

        if uu == ub == 0.0:
            raise CanteraError("Mass flow for reactants and/or products "
                               "must be positive")

        locs = np.array([0.0, 0.4, 0.6, 1.0])
        self.flame.set_profile('T', locs, [Tu, Tu, Teq, Tb])
        for k in range(self.gas.n_species):
            self.flame.set_profile(self.gas.species_name(k), locs,
                                   [Yu[k], Yu[k], Yeq[k], Yb[k]])

        # estimate strain rate
        self.gas.TPY = Teq, self.flame.P, Yeq
        zz = self.flame.grid
        dz = zz[-1] - zz[0]
        a = (uu + ub)/dz
        L = - 0.5 * (rhou + rhob) * a**2
        # estimate stagnation point
        x0 = rhou*uu * dz / (rhou*uu + rhob*ub)

        self.flame.set_profile("velocity", [0.0, 1.0], [uu, -ub])
        self.flame.set_profile("spreadRate", [0.0, x0/dz, 1.0], [0.0, a, 0.0])
        self.flame.set_profile("Lambda", [0.0, 1.0], [L, L])


class CounterflowTwinPremixedFlame(FlameBase):
    """
    A twin premixed counterflow flame. Two opposed jets of the same composition
    shooting into each other.
    """
    __slots__ = ('reactants', 'flame', 'products')

    def __init__(self, gas, grid=None, width=None):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            Array of initial grid points. Not recommended unless solving only on
            a fixed grid; Use the `width` parameter instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.

        A domain of class `AxisymmetricFlow` named ``flame`` will be created to
        represent the flame. The three domains comprising the stack are stored as
        ``self.reactants``, ``self.flame``, and ``self.products``.
        """
        self.reactants = Inlet1D(name='reactants', phase=gas)
        self.reactants.T = gas.T

        #: `AxisymmetricFlow` domain representing the flame
        self.flame = AxisymmetricFlow(gas, name='flame')

        #The right boundary is a symmetry plane
        self.products = SymmetryPlane1D(name='products', phase=gas)

        if width is not None:
            if grid is not None:
                raise ValueError("'grid' and 'width' arguments are mutually exclusive")
            # Create grid points aligned with initial guess profile
            grid = np.array([0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]) * width

        super().__init__((self.reactants, self.flame, self.products), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.reactants.X = gas.X

    def set_initial_guess(self, data=None, group=None):
        """
        Set the initial guess for the solution based on an equilibrium solution.
        Alternatively, a previously calculated result can be supplied as an
        initial guess via 'data' and 'key' inputs (see
        `FlameBase.set_initial_guess`).
        """
        super().set_initial_guess(data=data, group=group)
        if data:
            return

        if self.reactants.mdot == 0:
            raise CanteraError("Reactants velocity must be positive")

        Yu = self.reactants.Y
        Tu = self.reactants.T
        self.gas.TPY = Tu, self.flame.P, Yu
        rhou = self.gas.density
        uu = self.reactants.mdot / rhou

        self.gas.equilibrate('HP')
        Tb = self.gas.T
        Yb = self.gas.Y

        locs = np.array([0.0, 0.4, 0.6, 1.0])
        self.flame.set_profile('T', locs, [Tu, Tu, Tb, Tb])
        for k in range(self.gas.n_species):
            self.flame.set_profile(self.gas.species_name(k), locs,
                                   [Yu[k], Yu[k], Yb[k], Yb[k]])

        # estimate strain rate
        zz = self.flame.grid
        dz = zz[-1] - zz[0]
        a = 2 * uu / dz
        L = - rhou * a**2

        self.flame.set_profile("velocity", [0.0, 1.0], [uu, 0])
        self.flame.set_profile("spreadRate", [0.0, 1.0], [0.0, a])
        self.flame.set_profile("Lambda", [0.0, 1.0], [L, L])
