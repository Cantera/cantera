import numpy as np
from ._cantera import *

try:
    # Python 2.7 or 3.2+
    from math import erf
except ImportError:
    from scipy.special import erf


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
        super(FlameBase, self).__init__(domains)
        self.gas = gas
        self.flame.P = gas.P

    def set_refine_criteria(self, ratio=10.0, slope=0.8, curve=0.8, prune=0.0):
        super(FlameBase, self).set_refine_criteria(self.flame, ratio, slope,
                                                   curve, prune)

    def set_profile(self, component, locations, values):
        super(FlameBase, self).set_profile(self.flame, component, locations,
                                           values)

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

    @property
    def L(self):
        """
        Array containing the radial pressure gradient (1/r)(dP/dr) [N/m^4] at
        each point. Note: This value is named 'lambda' in the C++ code.
        """
        return self.profile(self.flame, 'lambda')

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
        u = self.u
        V = self.V

        csvfile = open(filename, 'w')
        writer = csv.writer(csvfile)
        writer.writerow(['z (m)', 'u (m/s)', 'V (1/s)',
                         'T (K)', 'rho (kg/m3)'] + self.gas.species_names)
        for n in range(self.flame.n_points):
            self.set_gas_state(n)
            writer.writerow([z[n], u[n], V[n], T[n], self.gas.density] +
                            list(getattr(self.gas, species)))
        csvfile.close()
        if not quiet:
            print("Solution saved to '{0}'.".format(filename))


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
    __slots__ = ('inlet', 'outlet', 'flame')

    def __init__(self, gas, grid=None):
        """
        A domain of type FreeFlow named 'flame' will be created to represent
        the flame. The three domains comprising the stack are stored as
        ``self.inlet``, ``self.flame``, and ``self.outlet``.
        """
        self.inlet = Inlet1D(name='reactants', phase=gas)
        self.outlet = Outlet1D(name='products', phase=gas)
        self.flame = FreeFlow(gas, name='flame')

        super(FreeFlame, self).__init__((self.inlet, self.flame, self.outlet),
                                        gas, grid)

    def set_initial_guess(self):
        """
        Set the initial guess for the solution. The adiabatic flame
        temperature and equilibrium composition are computed for the inlet gas
        composition. The temperature profile rises linearly over 20% of the
        domain width to Tad, then is flat. The mass fraction profiles are set
        similarly.
        """
        super(FreeFlame, self).set_initial_guess()
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
    __slots__ = ('burner', 'flame', 'outlet')

    def __init__(self, gas, grid=None):
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
        self.burner = Inlet1D(name='burner', phase=gas)
        self.burner.T = gas.T
        self.outlet = Outlet1D(name='outlet', phase=gas)
        self.flame = AxisymmetricStagnationFlow(gas, name='flame')

        super(BurnerFlame, self).__init__((self.burner, self.flame, self.outlet),
                                          gas, grid)

    def set_initial_guess(self):
        """
        Set the initial guess for the solution. The adiabatic flame
        temperature and equilibrium composition are computed for the burner
        gas composition. The temperature profile rises linearly in the first
        20% of the flame to Tad, then is flat. The mass fraction profiles are
        set similarly.
        """
        super(BurnerFlame, self).set_initial_guess()

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
    __slots__ = ('fuel_inlet', 'flame', 'oxidizer_inlet')

    def __init__(self, gas, grid=None):
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
        self.fuel_inlet = Inlet1D(name='fuel_inlet', phase=gas)
        self.fuel_inlet.T = gas.T

        self.oxidizer_inlet = Inlet1D(name='oxidizer_inlet', phase=gas)
        self.oxidizer_inlet.T = gas.T

        self.flame = AxisymmetricStagnationFlow(gas, name='flame')

        super(CounterflowDiffusionFlame, self).__init__(
                (self.fuel_inlet, self.flame, self.oxidizer_inlet), gas, grid)

    def set_initial_guess(self, fuel, oxidizer='O2', stoich=None):
        """
        Set the initial guess for the solution. The fuel species must be
        specified:

        >>> f.set_initial_guess(fuel='CH4')

        The oxidizer and corresponding stoichiometry must be specified if it
        is not 'O2'. The initial guess is generated by assuming infinitely-
        fast chemistry.
        """

        super(CounterflowDiffusionFlame, self).set_initial_guess()

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
            zmix = 0.5 * (1.0 - erf(zeta))
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
            return - (self.u[-1] - self.u[0]) / self.grid[-1]

        elif definition == 'max':
            return np.max(np.abs(np.gradient(self.u) / np.gradient(self.grid)))

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

            d_u_d_z = np.gradient(self.u) / np.gradient(self.grid)
            phi = (self.X[self.gas.species_index(fuel)]
                   / self.X[self.gas.species_index(oxidizer)] * stoich)
            z_stoich = np.interp(-1., -phi, self.grid)
            return np.abs(np.interp(z_stoich, self.grid, d_u_d_z))

        elif definition == 'potential_flow_fuel':
            return np.sqrt(- self.L[0] / self.density[0])

        elif definition == 'potential_flow_oxidizer':
            return np.sqrt(- self.L[0] / self.density[-1])

        else:
            raise ValueError('Definition "' + definition + '" is not available')


class ImpingingJet(FlameBase):
    """An axisymmetric flow impinging on a surface at normal incidence."""
    __slots__ = ('inlet', 'flame', 'surface')

    def __init__(self, gas, grid=None, surface=None):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            Array of initial grid points
        :param surface:
            A Kinetics object used to compute any surface reactions.

        A domain of class `AxisymmetricStagnationFlow` named ``flame`` will be
        created to represent the flow. The three domains comprising the stack
        are stored as ``self.inlet``, ``self.flame``, and ``self.surface``.
        """
        self.inlet = Inlet1D(name='inlet', phase=gas)
        self.inlet.T = gas.T
        self.flame = AxisymmetricStagnationFlow(gas, name='flame')

        if surface is None:
            self.surface = Surface1D(name='surface', phase=gas)
            self.surface.T = gas.T
        else:
            self.surface = ReactingSurface1D(name='surface', phase=gas)
            self.surface.set_kinetics(surface)
            self.surface.T = surface.T

        super(ImpingingJet, self).__init__(
                (self.inlet, self.flame, self.surface), gas, grid)

    def set_initial_guess(self, products='inlet'):
        """
        Set the initial guess for the solution. If products = 'equil', then
        the equilibrium composition at the adiabatic flame temperature will be
        used to form the initial guess. Otherwise the inlet composition will
        be used.
        """
        super(ImpingingJet, self).set_initial_guess()

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
        self.set_profile('u', locs, [u0, 0.0])
        self.set_profile('V', locs, [0.0, 0.0])


class CounterflowPremixedFlame(FlameBase):
    """ A premixed counterflow flame """
    __slots__ = ('reactants', 'flame', 'products')

    def __init__(self, gas, grid=None):
        """
        :param gas:
            `Solution` (using the IdealGas thermodynamic model) used to
            evaluate all gas properties and reaction rates.
        :param grid:
            Array of initial grid points

        A domain of class `AxisymmetricStagnationFlow` named ``flame`` will
        be created to represent the flame. The three domains comprising the
        stack are stored as ``self.reactants``, ``self.flame``, and
        ``self.products``.
        """
        self.reactants = Inlet1D(name='reactants', phase=gas)
        self.reactants.T = gas.T

        self.products = Inlet1D(name='products', phase=gas)
        self.products.T = gas.T

        self.flame = AxisymmetricStagnationFlow(gas, name='flame')

        super(CounterflowPremixedFlame, self).__init__(
                (self.reactants, self.flame, self.products), gas, grid)

    def set_initial_guess(self, equilibrate=True):
        """
        Set the initial guess for the solution.

        If `equilibrate` is True, then the products composition and temperature
        will be set to the equilibrium state of the reactants mixture.
        """

        super(CounterflowPremixedFlame, self).set_initial_guess()

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

        self.set_profile('u', [0.0, 1.0], [uu, -ub])
        self.set_profile('V', [0.0, x0/dz, 1.0], [0.0, a, 0.0])
