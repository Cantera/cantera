# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import numpy as np
from ._cantera import *
from .composite import Solution, SolutionArray
import csv as _csv
from math import erf
from email.utils import formatdate

# avoid explicit dependence of cantera on pandas
try:
    import pandas as _pandas
except ImportError as err:
    _pandas = err


class FlameBase(Sim1D):
    """ Base class for flames with a single flow domain """
    __slots__ = ('gas',)
    _extra = () # extra columns used for saving/restoring of simulation data

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

    def set_profile(self, component, locations, values):
        """
        Set an initial estimate for a profile of one component.

        :param component:
            component name or index
        :param positions:
            sequence of relative positions, from 0 on the left to 1 on the right
        :param values:
            sequence of values at the relative positions specified in *positions*

        >>> f.set_profile('T', [0.0, 0.2, 1.0], [400.0, 800.0, 1500.0])
        """
        super().set_profile(self.flame, component, locations, values)

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
        return self.gas.transport_model

    @transport_model.setter
    def transport_model(self, model):
        self.gas.transport_model = model
        self.flame.set_transport(self.gas)

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
    def radiation_enabled(self):
        """
        Get/Set whether or not to include radiative heat transfer
        """
        return self.flame.radiation_enabled

    @radiation_enabled.setter
    def radiation_enabled(self, enable):
        self.flame.radiation_enabled = enable

    def set_boundary_emissivities(self, e_left, e_right):
        self.flame.set_boundary_emissivities(e_left, e_right)

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
        return self.profile(self.flame, 'T')

    @property
    def u(self):
        """
        Array containing the velocity [m/s] normal to the flame at each point.

        .. deprecated:: 2.5

             To be deprecated with version 2.5, and removed thereafter.
             Replaced by property `velocity`.
        """
        warnings.warn("To be removed after Cantera 2.5. "
                      "Replaced by property 'velocity'",
                      DeprecationWarning)
        return self.profile(self.flame, 'velocity')

    @property
    def velocity(self):
        """
        Array containing the velocity [m/s] normal to the flame at each point.
        """
        return self.profile(self.flame, 'velocity')

    @property
    def V(self):
        """
        Array containing the tangential velocity gradient [1/s] at each point.
        """
        return self.profile(self.flame, 'V')

    @property
    def L(self):
        """
        Array containing the radial pressure gradient (1/r)(dP/dr) [N/m^4] at
        each point. Note: This value is named 'lambda' in the C++ code.
        """
        return self.profile(self.flame, 'lambda')

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
            self.set_gas_state(i)
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
            self.set_gas_state(i)
            vals[i] = self.gas.elemental_mole_fraction(m)
        return vals

    def solution(self, component, point=None):
        """
        Get the solution at one point or for the full flame domain (if
        `point=None`) for the specified *component*. The *component* can be
        specified by name or index.
        """
        if point is None:
            return self.profile(self.flame, component)
        else:
            return self.value(self.flame, component, point)

    def set_gas_state(self, point):
        """
        Set the state of the the Solution object used for calculations,
        `self.gas`, to the temperature and composition at the point with index
        *point*.
        """
        k0 = self.flame.component_index(self.gas.species_name(0))
        Y = [self.solution(k, point)
             for k in range(k0, k0 + self.gas.n_species)]
        self.gas.set_unnormalized_mass_fractions(Y)
        self.gas.TP = self.value(self.flame, 'T', point), self.P

    def write_csv(self, filename, species='X', quiet=True):
        """
        Write the velocity, temperature, density, and species profiles
        to a CSV file.

        :param filename:
            Output file name
        :param species:
            Attribute to use obtaining species profiles, e.g. ``X`` for
            mole fractions or ``Y`` for mass fractions.
        """

        z = self.grid
        T = self.T
        u = self.velocity
        V = self.V

        with open(filename, 'w', newline='') as csvfile:
            writer = _csv.writer(csvfile)
            writer.writerow(['z (m)', 'u (m/s)', 'V (1/s)',
                            'T (K)', 'rho (kg/m3)'] + self.gas.species_names)
            for n in range(self.flame.n_points):
                self.set_gas_state(n)
                writer.writerow([z[n], u[n], V[n], T[n], self.gas.density] +
                                list(getattr(self.gas, species)))

        if not quiet:
            print("Solution saved to '{0}'.".format(filename))

    def to_solution_array(self):
        """
        Return the solution vector as a Cantera `SolutionArray` object.

        The `SolutionArray` has the following ``extra`` entries:
         * ``grid``: grid point positions along the flame [m]
         * ``velocity``: normal velocity [m/s]
         * ``gradient``: tangential velocity gradient [1/s] (if applicable)
         * ``lambda``: radial pressure gradient [N/m^4] (if applicable)
         * ``eField``: electric field strength (if applicable)
        """
        # create extra columns
        extra = {}
        for e in self._extra:
            if e == 'grid':
                val = self.grid
            elif e == 'gradient':
                val = self.profile(self.flame, 'V')
            else:
                val = self.profile(self.flame, e)
            extra[e] = np.hstack([np.nan, val, np.nan])
        if self.radiation_enabled:
            qdot = self.flame.radiative_heat_loss()
            extra['qdot'] = np.hstack([np.nan, qdot, np.nan])

        # consider inlet boundaries
        left = self.domains[0]  # left boundary is always an inlet
        if isinstance(self.domains[2], Inlet1D):
            right = self.domains[2]
            n_arr = self.flame.n_points + 2
        else:
            right = None
            n_arr = self.flame.n_points + 1
            for e in extra:
                extra[e] = extra[e][:-1]

        # create solution array object
        arr = SolutionArray(self.gas, n_arr, extra=extra)

        # add left boundary
        self.gas.TPY = left.T, self.P, left.Y
        arr[0].TPY = self.gas.TPY
        arr._extra['grid'][0] = -np.inf
        arr._extra['velocity'][0] = left.mdot/self.gas.density

        # retrieve species concentrations and set states
        for n in range(self.flame.n_points):
            self.set_gas_state(n)
            arr[n + 1].TPY = self.gas.T, self.P, self.gas.Y

        # add right boundary
        if right:
            self.gas.TPY = right.T, self.P, right.Y
            arr[-1].TPY = self.gas.TPY
            arr._extra['grid'][-1] = np.inf
            arr._extra['velocity'][-1] = -right.mdot/self.gas.density

        return arr

    def from_solution_array(self, arr, restore_boundaries=True):
        """
        Restore the solution vector from a Cantera `SolutionArray` object.

        :param arr:
            SolutionArray to be restored
        :param restore_boundaries:
            Boolean flag to indicate whether boundaries should be restored
            (default is ``True``)

        The `SolutionArray` requires the following ``extra`` entries:
         * ``grid``: grid point positions along the flame [m]
         * ``velocity``: normal velocity [m/s]
         * ``gradient``: tangential velocity gradient [1/s] (if applicable)
         * ``lambda``: radial pressure gradient [N/m^4] (if applicable)
         * ``eField``: electric field strength (if applicable)
        """
        # extent (indices) of flame domain
        idx = np.isfinite(arr.grid)
        grid = arr.grid[idx]

        # restore grid
        self.flame.grid = grid
        self._get_initial_solution()
        xi = (grid - grid[0]) / (grid[-1] - grid[0])

        # restore temperature and 'extra' profiles
        self.set_profile('T', xi, arr.T[idx])
        for e in self._extra:
            val = getattr(arr, e)[idx]
            if e in ['grid', 'qdot']:
                pass
            elif e == 'gradient':
                self.set_profile('V', xi, val)
            else:
                self.set_profile(e, xi, val)

        # restore species profiles
        X = arr.X[idx, :]
        for i, spc in enumerate(self.gas.species_names):
            self.set_profile(spc, xi, X[:, i])

        # restore pressure
        self.P = arr.P[0]

        # restore boundaries
        if restore_boundaries:

            # left boundary
            left = self.domains[0]
            left.T = arr[0].T
            left.Y = arr[0].Y
            left.mdot = arr.velocity[0] * arr[0].density

            # right boundary
            if np.isinf(arr.grid[-1]):
                right = self.domains[2]
                right.T = arr[-1].T
                right.Y = arr[-1].Y
                right.mdot = -arr.velocity[-1] * arr[-1].density

    def to_pandas(self, species='X'):
        """
        Return the solution vector as a `pandas.DataFrame`.

        :param species:
            Attribute to use obtaining species profiles, e.g. ``X`` for
            mole fractions or ``Y`` for mass fractions.

        This method uses `to_solution_array` and requires a working pandas
        installation. Use pip or conda to install `pandas` to enable this
        method.
        """
        cols = ('extra', 'T', 'D', species)
        return self.to_solution_array().to_pandas(cols=cols)

    def from_pandas(self, df, restore_boundaries=True):
        """
        Restore the solution vector from a `pandas.DataFrame`.

        :param df:
            `pandas.DataFrame` containing data to be restored
        :param restore_boundaries:
            Boolean flag to indicate whether boundaries should be restored
            (default is ``True``)

        This method is intendend for loading of data that were previously
        exported by `to_pandas`. The method uses `from_solution_array` and
        requires a working pandas installation. The package 'pandas' can be
        installed using pip or conda.
        """
        arr = SolutionArray(self.gas, extra=self._extra)
        arr.from_pandas(df)
        self.from_solution_array(arr, restore_boundaries=restore_boundaries)

    def write_hdf(self, filename, key=None, species='X',
                  mode=None, complevel=None):
        """
        Write the solution vector to a HDF container file. Note that it is
        possible to write multiple data entries to a single HDF container file.
        Simulation settings are stored in tabular form as a separate HDF group
        named ``settings``.

        :param filename:
            HDF container file containing data to be restored
        :param key:
            String identifying the HDF group containing the data. The default
            is the name of the simulated configuration, e.g. ``FreeFlame``.
        :param species:
            Attribute to use obtaining species profiles, e.g. ``X`` for
            mole fractions or ``Y`` for mass fractions.
        :param mode:
            Mode to open file (see `pandas.DataFrame.to_hdf`)
        :param complevel:
            Compression level (see `pandas.DataFrame.to_hdf`)

        The method exports data using `SolutionArray.write_hdf` via
        `to_solution_array` and requires working installations of pandas and
        PyTables. These packages can be installed using pip (`pandas` and
        `tables`) or conda (`pandas` and `pytables`).
        """
        if not key:
            key = type(self).__name__

        # save data
        cols = ('extra', 'T', 'D', species)
        self.to_solution_array().write_hdf(filename, cols=cols, key=key,
                                           mode=mode, complevel=complevel)

        # convert simulation settings to tabular format
        df = _pandas.DataFrame()
        df['key'] = [key]
        df['date'] = formatdate(localtime=True)
        for key, val in self.settings.items():
            df[key] = [val]
        df.set_index('key')

        # store settings to HDF container file as a separate group
        df.to_hdf(filename, key='settings', append=True)

    def read_hdf(self, filename, key=None, restore_boundaries=True):
        """
        Restore the solution vector from a HDF container file.

        :param filename:
            HDF container file containing data to be restored
        :param key:
            String identifying the HDF group containing the data
        :param restore_boundaries:
            Boolean flag to indicate whether boundaries should be restored
            (default is ``True``)

        The method imports data using `SolutionArray.read_hdf` via
        `from_solution_array` and requires working installations of pandas and
        PyTables. These packages can be installed using pip (`pandas` and
        `tables`) or conda (`pandas` and `pytables`).
        """
        arr = SolutionArray(self.gas, extra=self._extra)
        arr.read_hdf(filename, key=key)
        self.from_solution_array(arr, restore_boundaries=restore_boundaries)

    @property
    def settings(self):
        """ Return a dictionary listing simulation settings """
        out = {'type': type(self).__name__}
        out['transport_model'] = self.transport_model
        out['energy_enabled'] = self.energy_enabled
        out['soret_enabled'] = self.soret_enabled
        out['radiation_enabled'] = self.radiation_enabled
        epsilon = self.flame.boundary_emissivities
        out['emissivity_left'] = epsilon[0]
        out['emissivity_right'] = epsilon[1]
        out['fixed_temperature'] = self.fixed_temperature
        out.update(self.get_refine_criteria())
        out['max_time_step_count'] = self.max_time_step_count
        out['max_grid_points'] = self.get_max_grid_points(self.flame)

        # add tolerance settings
        tols = {'steady_abstol': self.flame.steady_abstol(),
                'steady_reltol': self.flame.steady_reltol(),
                'transient_abstol': self.flame.transient_abstol(),
                'transient_reltol': self.flame.transient_reltol()}
        comp = np.array(self.flame.component_names)
        for tname, tol in tols.items():
            # add mode (most frequent tolerance setting)
            values, counts = np.unique(tol, return_counts=True)
            ix = np.argmax(counts)
            out.update({tname: values[ix]})

            # add values deviating from mode
            ix = np.logical_not(np.isclose(tol, values[ix]))
            out.update({'{}_{}'.format(tname, c)
                        for c, t in zip(comp[ix], tol[ix])})

        return out

    def _load_restart_data(self, source, **kwargs):
        """ Load data for restart (called by set_initial_guess) """
        if isinstance(source, SolutionArray):
            # already a solution array
            return source
        elif isinstance(source, str):
            # source identifies a HDF file
            arr = SolutionArray(self.gas, extra=self._extra)
            arr.read_hdf(source, **kwargs)
            return arr
        else:
            # source is a pandas DataFrame
            arr = SolutionArray(self.gas, extra=self._extra)
            arr.from_pandas(source)
            return arr


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

    doc = _trim(getattr(Solution, attr).__doc__) +'\n' + extradoc
    return property(getter, doc=doc)

# Add scalar properties to FlameBase
for _attr in ['density', 'density_mass', 'density_mole', 'volume_mass',
              'volume_mole', 'int_energy_mole', 'int_energy_mass', 'h',
              'enthalpy_mole', 'enthalpy_mass', 's', 'entropy_mole',
              'entropy_mass', 'g', 'gibbs_mole', 'gibbs_mass', 'cv',
              'cv_mole', 'cv_mass', 'cp', 'cp_mole', 'cp_mass',
              'isothermal_compressibility', 'thermal_expansion_coeff',
              'viscosity', 'thermal_conductivity', 'heat_release_rate']:
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
              'destruction_rates', 'net_production_rates', 'mix_diff_coeffs',
              'mix_diff_coeffs_mass', 'mix_diff_coeffs_mole', 'thermal_diff_coeffs']:
    setattr(FlameBase, _attr, _array_property(_attr, 'n_species'))

# Add properties with values for each reaction
for _attr in ['forward_rates_of_progress', 'reverse_rates_of_progress', 'net_rates_of_progress',
              'equilibrium_constants', 'forward_rate_constants', 'reverse_rate_constants',
              'delta_enthalpy', 'delta_gibbs', 'delta_entropy',
              'delta_standard_enthalpy', 'delta_standard_gibbs',
              'delta_standard_entropy', 'heat_production_rates']:
    setattr(FlameBase, _attr, _array_property(_attr, 'n_reactions'))


class FreeFlame(FlameBase):
    """A freely-propagating flat flame."""
    __slots__ = ('inlet', 'outlet', 'flame')
    _extra = ('grid', 'velocity')

    def __init__(self, gas, grid=None, width=None):
        """
        A domain of type IdealGasFlow named 'flame' will be created to represent
        the flame and set to free flow. The three domains comprising the stack
        are stored as ``self.inlet``, ``self.flame``, and ``self.outlet``.

        :param grid:
            A list of points to be used as the initial grid. Not recommended
            unless solving only on a fixed grid; Use the `width` parameter
            instead.
        :param width:
            Defines a grid on the interval [0, width] with internal points
            determined automatically by the solver.
        """
        self.inlet = Inlet1D(name='reactants', phase=gas)
        self.outlet = Outlet1D(name='products', phase=gas)
        if not hasattr(self, 'flame'):
            # Create flame domain if not already instantiated by a child class
            self.flame = IdealGasFlow(gas, name='flame')
            self.flame.set_free_flow()

        if width is not None:
            grid = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]) * width

        super().__init__((self.inlet, self.flame, self.outlet), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.inlet.T = gas.T
        self.inlet.X = gas.X

    def set_initial_guess(self, locs=[0.0, 0.3, 0.5, 1.0], data=None, key=None):
        """
        Set the initial guess for the solution. The adiabatic flame
        temperature and equilibrium composition are computed for the inlet gas
        composition.

        :param locs:
            A list of four locations to define the temperature and mass fraction
            profiles. Profiles rise linearly between the second and third
            location. Locations are given as a fraction of the entire domain
        :param data:
            Restart data, which are typically based on an earlier simulation
            result. Restart data may be specified using a SolutionArray,
            pandas' DataFrame, or a saved HDF container file. Note that restart
            data do not overwrite boundary conditions. DataFrame and HDF input
            require working installations of pandas and PyTables. These packages
            can be installed using pip (`pandas` and `tables`) or conda
            (`pandas` and `pytables`).
        :param key:
            Group identifier within a HDF container file (only used in
            combination with HDF restart data).
        """
        super().set_initial_guess()
        if data:
            data = self._load_restart_data(data, key=key)
            data.TP = data.T + self.inlet.T - data.T[0], self.P
            self.from_solution_array(data, restore_boundaries=False)

            # set fixed temperature
            Tmid = .75 * data.T[0] + .25 * data.T[-1]
            i = np.flatnonzero(data.T < Tmid)[-1]
            self.fixed_temperature = data.T[i]

            return

        self.gas.TPY = self.inlet.T, self.P, self.inlet.Y

        if not self.inlet.mdot:
            # nonzero initial guess increases likelihood of convergence
            self.inlet.mdot = 1.0 * self.gas.density

        Y0 = self.inlet.Y
        u0 = self.inlet.mdot/self.gas.density
        T0 = self.inlet.T

        # get adiabatic flame temperature and composition
        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y
        u1 = self.inlet.mdot/self.gas.density

        self.set_profile('velocity', locs, [u0, u0, u1, u1])
        self.set_profile('T', locs, [T0, T0, Teq, Teq])

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
            self.set_profile(self.gas.species_name(n),
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


class IonFlameBase(FlameBase):

    def write_csv(self, filename, species='X', quiet=True):
        """
        Write the velocity, temperature, density, electric potential,
        electric field strength, and species profiles to a CSV file.
        :param filename:
            Output file name
        :param species:
            Attribute to use obtaining species profiles, e.g. ``X`` for
            mole fractions or ``Y`` for mass fractions.
        """
        z = self.grid
        T = self.T
        u = self.velocity
        V = self.V
        E = self.E

        with open(filename, 'w', newline='') as csvfile:
            writer = _csv.writer(csvfile)
            writer.writerow(['z (m)', 'velocity (m/s)', 'V (1/s)', 'T (K)',
                            'E (V/m)', 'rho (kg/m3)'] + self.gas.species_names)
            for n in range(self.flame.n_points):
                self.set_gas_state(n)
                writer.writerow([z[n], u[n], V[n], T[n], E[n], self.gas.density] +
                                list(getattr(self.gas, species)))

        if not quiet:
            print("Solution saved to '{0}'.".format(filename))

    @property
    def electric_field_enabled(self):
        """ Get/Set whether or not to solve the Poisson's equation."""
        return self.flame.electric_field_enabled

    @electric_field_enabled.setter
    def electric_field_enabled(self, enable):
        self.flame.electric_field_enabled = enable

    @property
    def E(self):
        """
        Array containing the electric field strength at each point.
        """
        return self.profile(self.flame, 'eField')

    def solve(self, loglevel=1, refine_grid=True, auto=False, stage=1, enable_energy=True):
        self.flame.set_solving_stage(stage)
        if stage == 1:
            super().solve(loglevel, refine_grid, auto)
        if stage == 2:
            self.poisson_enabled = True
            super().solve(loglevel, refine_grid, auto)


class IonFreeFlame(IonFlameBase, FreeFlame):
    """A freely-propagating flame with ionized gas."""
    __slots__ = ('inlet', 'outlet', 'flame')
    _extra = ('grid', 'velocity', 'eField')

    def __init__(self, gas, grid=None, width=None):
        if not hasattr(self, 'flame'):
            # Create flame domain if not already instantiated by a child class
            self.flame = IonFlow(gas, name='flame')
            self.flame.set_free_flow()

        super().__init__(gas, grid, width)


class BurnerFlame(FlameBase):
    """A burner-stabilized flat flame."""
    __slots__ = ('burner', 'flame', 'outlet')
    _extra = ('grid', 'velocity')

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

        A domain of class `IdealGasFlow` named ``flame`` will be created to
        represent the flame and set to axisymmetric stagnation flow. The three
        domains comprising the stack are stored as ``self.burner``,
        ``self.flame``, and ``self.outlet``.
        """
        self.burner = Inlet1D(name='burner', phase=gas)
        self.outlet = Outlet1D(name='outlet', phase=gas)
        if not hasattr(self, 'flame'):
            # Create flame domain if not already instantiated by a child class
            self.flame = IdealGasFlow(gas, name='flame')
            self.flame.set_axisymmetric_flow()

        if width is not None:
            grid = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]) * width

        super().__init__((self.burner, self.flame, self.outlet), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.burner.T = gas.T
        self.burner.X = gas.X

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
        self.set_profile('velocity', locs, [u0, u1, u1])
        self.set_profile('T', locs, [T0, Teq, Teq])
        for n in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(n),
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
            self.set_flat_profile(self.flame, 'T', self.T[0])
            for k,spec in enumerate(self.gas.species_names):
                self.set_flat_profile(self.flame, spec, self.burner.Y[k])

            self.set_steady_callback(original_callback)
            super().solve(loglevel, False, False)
            if loglevel > 0:
                print('Flame has blown off of burner (non-reacting solution)')

        self.set_steady_callback(original_callback)


class IonBurnerFlame(IonFlameBase, BurnerFlame):
    """A burner-stabilized flat flame with ionized gas."""
    __slots__ = ('burner', 'flame', 'outlet')
    _extra = ('grid', 'velocity', 'eField')

    def __init__(self, gas, grid=None, width=None):
        if not hasattr(self, 'flame'):
            # Create flame domain if not already instantiated by a child class
            self.flame = IonFlow(gas, name='flame')
            self.flame.set_axisymmetric_flow()

        super().__init__(gas, grid, width)


class CounterflowDiffusionFlame(FlameBase):
    """ A counterflow diffusion flame """
    __slots__ = ('fuel_inlet', 'flame', 'oxidizer_inlet')
    _extra = ('grid', 'velocity', 'gradient', 'lambda')

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

        A domain of class `IdealGasFlow` named ``flame`` will be created to
        represent the flame and set to axisymmetric stagnation flow. The three
        domains comprising the stack are stored as ``self.fuel_inlet``,
        ``self.flame``, and ``self.oxidizer_inlet``.
        """
        self.fuel_inlet = Inlet1D(name='fuel_inlet', phase=gas)
        self.fuel_inlet.T = gas.T

        self.oxidizer_inlet = Inlet1D(name='oxidizer_inlet', phase=gas)
        self.oxidizer_inlet.T = gas.T

        self.flame = IdealGasFlow(gas, name='flame')
        self.flame.set_axisymmetric_flow()

        if width is not None:
            grid = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) * width

        super().__init__((self.fuel_inlet, self.flame, self.oxidizer_inlet), gas, grid)

    def set_initial_guess(self):
        """
        Set the initial guess for the solution. The initial guess is generated
        by assuming infinitely-fast chemistry.
        """

        super().set_initial_guess()

        moles = lambda el: (self.gas.elemental_mass_fraction(el) /
                            self.gas.atomic_weight(el))

        # Compute stoichiometric mixture composition
        Yin_f = self.fuel_inlet.Y
        self.gas.TPY = self.fuel_inlet.T, self.P, Yin_f
        mdotf = self.fuel_inlet.mdot
        u0f = mdotf / self.gas.density
        T0f = self.fuel_inlet.T

        sFuel = moles('O')
        if 'C' in self.gas.element_names:
            sFuel -= 2 * moles('C')
        if 'H' in self.gas.element_names:
            sFuel -= 0.5 * moles('H')

        Yin_o = self.oxidizer_inlet.Y
        self.gas.TPY = self.oxidizer_inlet.T, self.P, Yin_o
        mdoto = self.oxidizer_inlet.mdot
        u0o = mdoto / self.gas.density
        T0o = self.oxidizer_inlet.T

        sOx = moles('O')
        if 'C' in self.gas.element_names:
            sOx -= 2 * moles('C')
        if 'H' in self.gas.element_names:
            sOx -= 0.5 * moles('H')

        zst = 1.0 / (1 - sFuel / sOx)
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

        x0 = np.sqrt(mdotf*u0f) * dz / (np.sqrt(mdotf*u0f) + np.sqrt(mdoto*u0o))
        nz = len(zz)

        Y = np.zeros((nz, self.gas.n_species))
        T = np.zeros(nz)
        for j in range(nz):
            x = zz[j] - zz[0]
            zeta = f * (x - x0)
            zmix = 0.5 * (1.0 - erf(zeta))
            if zmix > zst:
                Y[j] = Yeq + (Yin_f - Yeq) * (zmix - zst) / (1.0 - zst)
                T[j] = Teq + (T0f - Teq) * (zmix - zst) / (1.0 - zst)
            else:
                Y[j] = Yin_o + zmix * (Yeq - Yin_o) / zst
                T[j] = T0o + (Teq - T0o) * zmix / zst

        T[0] = T0f
        T[-1] = T0o
        zrel = (zz - zz[0])/dz

        self.set_profile('velocity', [0.0, 1.0], [u0f, -u0o])
        self.set_profile('V', [0.0, x0/dz, 1.0], [0.0, a, 0.0])
        self.set_profile('T', zrel, T)
        for k,spec in enumerate(self.gas.species_names):
            self.set_profile(spec, zrel, Y[:,k])

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
        super().solve(loglevel, refine_grid, auto)
        # Do some checks if loglevel is set
        if loglevel > 0:
            if self.extinct():
                print('WARNING: Flame is extinct.')

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
        :param fuel: The fuel species. Used only if *definition* is
            ``stoichiometric``.
        :param oxidizer: The oxidizer species, default ``O2``. Used only if
            *definition* is ``stoichiometric``.
        :param stoich: The molar stoichiometric oxidizer-to-fuel ratio.
            Can be omitted if the oxidizer is ``O2``. Used only if *definition*
            is ``stoichiometric``.

        The parameter *definition* sets the method to compute the strain rate.
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

            This method uses the additional keyword arguments *fuel*,
            *oxidizer*, and *stoich*.

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
            return np.sqrt(- self.L[0] / self.density[0])

        elif definition == 'potential_flow_oxidizer':
            return np.sqrt(- self.L[0] / self.density[-1])

        else:
            raise ValueError('Definition "' + definition + '" is not available')

    def mixture_fraction(self, m):
        r"""
        Compute the mixture fraction based on element *m*

        The mixture fraction is computed from the elemental mass fraction of
        element *m*, normalized by its values on the fuel and oxidizer
        inlets:

        .. math:: Z = \frac{Z_{\mathrm{mass},m}(z) -
                            Z_{\mathrm{mass},m}(z_\mathrm{oxidizer})}
                           {Z_{\mathrm{mass},m}(z_\mathrm{fuel}) -
                            Z_{\mathrm{mass},m}(z_\mathrm{oxidizer})}

        :param m:
            The element based on which the mixture fraction is computed,
            may be specified by name or by index

        >>> f.mixture_fraction('H')
        """
        emf = self.elemental_mass_fraction(m)
        return (emf - emf[-1]) / (emf[0] - emf[-1])


class ImpingingJet(FlameBase):
    """An axisymmetric flow impinging on a surface at normal incidence."""
    __slots__ = ('inlet', 'flame', 'surface')
    _extra = ('grid', 'velocity', 'gradient', 'lambda')

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

        A domain of class `IdealGasFlow` named ``flame`` will be created to
        represent the flame and set to axisymmetric stagnation flow. The three
        domains comprising the stack are stored as ``self.inlet``,
        ``self.flame``, and ``self.surface``.
        """
        self.inlet = Inlet1D(name='inlet', phase=gas)
        self.flame = IdealGasFlow(gas, name='flame')
        self.flame.set_axisymmetric_flow()

        if width is not None:
            grid = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]) * width

        if surface is None:
            self.surface = Surface1D(name='surface', phase=gas)
            self.surface.T = gas.T
        else:
            self.surface = ReactingSurface1D(name='surface', phase=gas)
            self.surface.set_kinetics(surface)
            self.surface.T = surface.T

        super().__init__((self.inlet, self.flame, self.surface), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.inlet.T = gas.T
        self.inlet.X = gas.X

    def set_initial_guess(self, products='inlet'):
        """
        Set the initial guess for the solution. If products = 'equil', then
        the equilibrium composition at the adiabatic flame temperature will be
        used to form the initial guess. Otherwise the inlet composition will
        be used.
        """
        super().set_initial_guess(products=products)

        Y0 = self.inlet.Y
        T0 = self.inlet.T
        self.gas.TPY = T0, self.flame.P, Y0
        u0 = self.inlet.mdot / self.gas.density

        if products == 'equil':
            self.gas.equilibrate('HP')
            Teq = self.gas.T
            Yeq = self.gas.Y
            locs = np.array([0.0, 0.3, 0.7, 1.0])
            self.set_profile('T', locs, [T0, Teq, Teq, self.surface.T])
            for k in range(self.gas.n_species):
                self.set_profile(self.gas.species_name(k), locs,
                                 [Y0[k], Yeq[k], Yeq[k], Yeq[k]])
        else:
            locs = np.array([0.0, 1.0])
            self.set_profile('T', locs, [T0, self.surface.T])
            for k in range(self.gas.n_species):
                self.set_profile(self.gas.species_name(k), locs,
                                 [Y0[k], Y0[k]])

        locs = np.array([0.0, 1.0])
        self.set_profile('velocity', locs, [u0, 0.0])
        self.set_profile('V', locs, [0.0, 0.0])


class CounterflowPremixedFlame(FlameBase):
    """ A premixed counterflow flame """
    __slots__ = ('reactants', 'flame', 'products')
    _extra = ('grid', 'velocity', 'gradient', 'lambda')

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

        A domain of class `IdealGasFlow` named ``flame`` will be created to
        represent the flame and set to axisymmetric stagnation flow. The three
        domains comprising the stack are stored as ``self.reactants``,
        ``self.flame``, and ``self.products``.
        """
        self.reactants = Inlet1D(name='reactants', phase=gas)
        self.reactants.T = gas.T

        self.products = Inlet1D(name='products', phase=gas)
        self.products.T = gas.T

        self.flame = IdealGasFlow(gas, name='flame')
        self.flame.set_axisymmetric_flow()

        if width is not None:
            # Create grid points aligned with initial guess profile
            grid = np.array([0.0, 0.3, 0.5, 0.7, 1.0]) * width

        super().__init__((self.reactants, self.flame, self.products), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.reactants.X = gas.X

    def set_initial_guess(self, equilibrate=True):
        """
        Set the initial guess for the solution.

        If `equilibrate` is True, then the products composition and temperature
        will be set to the equilibrium state of the reactants mixture.
        """

        super().set_initial_guess(equilibrate=equilibrate)

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

        locs = np.array([0.0, 0.4, 0.6, 1.0])
        self.set_profile('T', locs, [Tu, Tu, Teq, Tb])
        for k in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(k), locs,
                             [Yu[k], Yu[k], Yeq[k], Yb[k]])

        # estimate strain rate
        self.gas.TPY = Teq, self.flame.P, Yeq
        zz = self.flame.grid
        dz = zz[-1] - zz[0]
        a = (uu + ub)/dz
        # estimate stagnation point
        x0 = rhou*uu * dz / (rhou*uu + rhob*ub)

        self.set_profile('velocity', [0.0, 1.0], [uu, -ub])
        self.set_profile('V', [0.0, x0/dz, 1.0], [0.0, a, 0.0])


class CounterflowTwinPremixedFlame(FlameBase):
    """
    A twin premixed counterflow flame. Two opposed jets of the same composition
    shooting into each other.
    """
    __slots__ = ('reactants', 'flame', 'products')
    _extra = ('grid', 'velocity', 'gradient', 'lambda')

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

        A domain of class `IdealGasFlow` named ``flame`` will be created to
        represent the flame and set to axisymmetric stagnation flow. The three
        domains comprising the stack are stored as ``self.reactants``,
        ``self.flame``, and ``self.products``.
        """
        self.reactants = Inlet1D(name='reactants', phase=gas)
        self.reactants.T = gas.T

        self.flame = IdealGasFlow(gas, name='flame')
        self.flame.set_axisymmetric_flow()

        #The right boundary is a symmetry plane
        self.products = SymmetryPlane1D(name='products', phase=gas)

        if width is not None:
            # Create grid points aligned with initial guess profile
            grid = np.array([0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]) * width

        super().__init__((self.reactants, self.flame, self.products), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain
        self.reactants.X = gas.X

    def set_initial_guess(self):
        """
        Set the initial guess for the solution.
        """
        super().set_initial_guess()

        Yu = self.reactants.Y
        Tu = self.reactants.T
        self.gas.TPY = Tu, self.flame.P, Yu
        uu = self.reactants.mdot / self.gas.density

        self.gas.equilibrate('HP')
        Tb = self.gas.T
        Yb = self.gas.Y

        locs = np.array([0.0, 0.4, 0.6, 1.0])
        self.set_profile('T', locs, [Tu, Tu, Tb, Tb])
        for k in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(k), locs,
                             [Yu[k], Yu[k], Yb[k], Yb[k]])

        # estimate strain rate
        zz = self.flame.grid
        dz = zz[-1] - zz[0]
        a = 2 * uu / dz

        self.set_profile('velocity', [0.0, 1.0], [uu, 0])
        self.set_profile('V', [0.0, 1.0], [0.0, a])
