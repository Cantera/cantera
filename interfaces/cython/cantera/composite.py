from ._cantera import *
import numpy as np

class Quantity(object):
    """
    A class representing a specific quantity of a `Solution`. In addition to the
    properties which can be computed for class `Solution`, class `Quantity`
    provides several additional capabilities. A `Quantity` object is created
    from a `Solution` with either the mass or number of moles specified::

        >>> gas = ct.Solution('gri30.xml')
        >>> gas.TPX = 300, 5e5, 'O2:1.0, N2:3.76'
        >>> q1 = ct.Quantity(gas, mass=5) # 5 kg of air

    The state of a `Quantity` can be changed in the same way as a `Solution`::

        >>> q1.TP = 500, 101325

    Quantities have properties which provide access to extensive properties::

        >>> q1.volume
        7.1105094
        >>> q1.enthalpy
        1032237.84

    The size of a `Quantity` can be changed by setting the mass or number of
    moles::

        >>> q1.moles = 3
        >>> q1.mass
        86.552196
        >>> q1.volume
        123.086

    or by multiplication::

        >>> q1 *= 2
        >>> q1.moles
        6.0

    Finally, Quantities can be added, providing an easy way of calculating the
    state resulting from mixing two substances::

        >>> q1.mass = 5
        >>> q2 = ct.Quantity(gas)
        >>> q2.TPX = 300, 101325, 'CH4:1.0'
        >>> q2.mass = 1
        >>> q3 = q1 + q2 # combine at constant UV
        >>> q3.T
        432.31234
        >>> q3.P
        97974.9871
        >>> q3.mole_fraction_dict()
        {'CH4': 0.26452900448117395,
         'N2': 0.5809602821745349,
         'O2': 0.1545107133442912}

    If a different property pair should be held constant when combining, this
    can be specified as follows::

        >>> q1.constant = q2.constant = 'HP'
        >>> q3 = q1 + q2 # combine at constant HP
        >>> q3.T
        436.03320
        >>> q3.P
        101325.0
    """
    def __init__(self, phase, mass=None, moles=None, constant='UV'):
        self.state = phase.TDY
        self._phase = phase

        # A unique key to prevent adding phases with different species
        # definitions
        self._id = hash((phase.name,) + tuple(phase.species_names))

        if mass is not None:
            self.mass = mass
        elif moles is not None:
            self.moles = moles
        else:
            self.mass = 1.0

        assert constant in ('TP','TV','HP','SP','SV','UV')
        self.constant = constant

    @property
    def phase(self):
        """
        Get the underlying `Solution` object, with the state set to match the
        wrapping `Quantity` object.
        """
        self._phase.TDY = self.state
        return self._phase

    @property
    def moles(self):
        """ Get/Set the number of moles [kmol] represented by the `Quantity`. """
        return self.mass / self.phase.mean_molecular_weight

    @moles.setter
    def moles(self, n):
        self.mass = n * self.phase.mean_molecular_weight

    @property
    def volume(self):
        """ Get the total volume [m^3] represented by the `Quantity`. """
        return self.mass * self.phase.volume_mass

    @property
    def int_energy(self):
        """ Get the total internal energy [J] represented by the `Quantity`. """
        return self.mass * self.phase.int_energy_mass

    @property
    def enthalpy(self):
        """ Get the total enthalpy [J] represented by the `Quantity`. """
        return self.mass * self.phase.enthalpy_mass

    @property
    def entropy(self):
        """ Get the total entropy [J/K] represented by the `Quantity`. """
        return self.mass * self.phase.entropy_mass

    @property
    def gibbs(self):
        """
        Get the total Gibbs free energy [J] represented by the `Quantity`.
        """
        return self.mass * self.phase.gibbs_mass

    def equilibrate(self, XY=None, *args, **kwargs):
        """
        Set the state to equilibrium. By default, the property pair
        `self.constant` is held constant. See `ThermoPhase.equilibrate`.
        """
        if XY is None:
            XY = self.constant
        self.phase.equilibrate(XY, *args, **kwargs)
        self.state = self._phase.TDY

    def __imul__(self, other):
        self.mass *= other
        return self

    def __mul__(self, other):
        return Quantity(self.phase, mass=self.mass * other)

    def __rmul__(self, other):
        return Quantity(self.phase, mass=self.mass * other)

    def __iadd__(self, other):
        if (self._id != other._id):
            raise ValueError('Cannot add Quantities with different phase '
                'definitions.')
        assert(self.constant == other.constant)
        a1,b1 = getattr(self.phase, self.constant)
        a2,b2 = getattr(other.phase, self.constant)
        m = self.mass + other.mass
        a = (a1 * self.mass + a2 * other.mass) / m
        b = (b1 * self.mass + b2 * other.mass) / m
        self._phase.Y = (self.Y * self.mass + other.Y * other.mass) / m
        setattr(self._phase, self.constant, (a,b))
        self.state = self._phase.TDY
        self.mass = m
        return self

    def __add__(self, other):
        newquantity = Quantity(self.phase, mass=self.mass, constant=self.constant)
        newquantity += other
        return newquantity

# Synonyms for total properties
Quantity.V = Quantity.volume
Quantity.U = Quantity.int_energy
Quantity.H = Quantity.enthalpy
Quantity.S = Quantity.entropy
Quantity.G = Quantity.gibbs

# Add properties to act as pass-throughs for attributes of class Solution
def _prop(attr):
    def getter(self):
        return getattr(self.phase, attr)

    def setter(self, value):
        setattr(self.phase, attr, value)
        self.state = self._phase.TDY

    return property(getter, setter, doc=getattr(Solution, attr).__doc__)

for _attr in dir(Solution):
    if _attr.startswith('_') or _attr in Quantity.__dict__ or _attr == 'state':
        continue
    else:
        setattr(Quantity, _attr, _prop(_attr))


class SolutionArray(object):
    def __init__(self, phase, shape, states=None):
        self._phase = phase
        if isinstance(shape, int):
            shape = (shape,)
        if states is not None:
            self._shape = states.shape[:-1]
            self._states = states
        else:
            self._shape = shape
            S = np.empty(shape + (2+self._phase.n_species,))
            S[:] = self._phase.state
            self._states = S

        self._indices = list(np.ndindex(self._shape))

    def __getitem__(self, index):
        states = self._states[index]
        shape = states.shape[:-1]
        return SolutionArray(self._phase, shape, states)

    def equilibrate(self, *args, **kwargs):
        """ See `ThermoPhase.equilibrate` """
        for index in self._indices:
            self._phase.state = self._states[index]
            self._phase.equilibrate(*args, **kwargs)
            self._states[index][:] = self._phase.state


def _make_functions():
    # this is wrapped in a function to avoid polluting the module namespace

    scalar = [
        # From ThermoPhase
        'mean_molecular_weight', 'P', 'T', 'density', 'density_mass',
        'density_mole', 'v', 'volume_mass', 'volume_mole', 'u',
        'int_energy_mole', 'int_energy_mass', 'h', 'enthalpy_mole',
        'enthalpy_mass', 's', 'entropy_mole', 'entropy_mass', 'g', 'gibbs_mole',
        'gibbs_mass', 'cv', 'cv_mole', 'cv_mass', 'cp', 'cp_mole', 'cp_mass',
        'critical_temperature', 'critical_pressure', 'critical_density',
        'P_sat', 'T_sat', 'isothermal_compressibility',
        'thermal_expansion_coeff', 'electric_potential',
        # From Transport
        'viscosity', 'electrical_conductivity', 'thermal_conductivity',
    ]
    n_species = [
        # from ThermoPhase
        'Y', 'X', 'concentrations', 'partial_molar_enthalpies',
        'partial_molar_entropies', 'partial_molar_int_energies',
        'chemical_potentials', 'electrochemical_potentials', 'partial_molar_cp',
        'partial_molar_volumes', 'standard_enthalpies_RT',
        'standard_entropies_R', 'standard_int_energies_RT', 'standard_gibbs_RT',
        'standard_cp_R',
        # From Kinetics
        'creation_rates', 'destruction_rates', 'net_production_rates',
        # From Transport
        'mix_diff_coeffs', 'mix_diff_coeffs_mass', 'mix_diff_coeffs_mole',
        'thermal_diff_coeffs'
    ]

    n_species2 = ['multi_diff_coeffs', 'binary_diff_coeffs']

    n_reactions = [
        'forward_rates_of_progress', 'reverse_rates_of_progress',
        'net_rates_of_progress', 'equilibrium_constants',
        'forward_rate_constants', 'reverse_rate_constants',
        'delta_enthalpy', 'delta_gibbs', 'delta_entropy',
        'delta_standard_enthalpy', 'delta_standard_gibbs',
        'delta_standard_entropy'
    ]
    state2 = ['TD', 'TP', 'UV', 'DP', 'HP', 'SP', 'SV']
    state3 = [
        'TDX', 'TDY', 'TPX', 'TPY', 'UVX', 'UVY', 'DPX', 'DPY', 'HPX', 'HPY',
        'SPX', 'SPY', 'SVX', 'SVY'
    ]
    call = ['elemental_mass_fraction', 'elemental_mole_fraction']

    passthrough = [
        # from ThermoPhase
        'name', 'ID', 'basis', 'n_elements', 'element_index',
        'element_name', 'element_names', 'atomic_weight', 'atomic_weights',
        'n_species', 'species_name', 'species_names', 'species_index',
        'species', 'n_atoms', 'molecular_weights', 'min_temp', 'max_temp',
        'reference_pressure',
        # From Kinetics
        'n_total_species', 'n_reactions', 'n_phases', 'reaction_phase_index',
        'kinetics_species_index', 'reaction', 'reactions', 'modify_reaction',
        'is_reversible', 'multiplier', 'set_multiplier', 'reaction_type',
        'reaction_equation', 'reactants', 'products', 'reaction_equations',
        'reactant_stoich_coeff', 'product_stoich_coeff',
        'reactant_stoich_coeffs', 'product_stoich_coeffs',
        # from Transport
        'transport_model',
    ]

    # Factory for creating properties which consist of a tuple of two variables,
    # e.g. 'TP' or 'SV'
    def state2_prop(name):
        def getter(self):
            a = np.empty(self._shape)
            b = np.empty(self._shape)
            for index in self._indices:
                self._phase.state = self._states[index]
                a[index], b[index] = getattr(self._phase, name)
            return a, b

        def setter(self, AB):
            assert len(AB) == 2, "Expected 2 elements, got {}".format(len(AB))
            A, B, _ = np.broadcast_arrays(AB[0], AB[1], self._states[...,0])
            for index in self._indices:
                self._phase.state = self._states[index]
                setattr(self._phase, name, (A[index], B[index]))
                self._states[index][:] = self._phase.state

        return property(getter, setter, doc=getattr(Solution, name).__doc__)

    for name in state2:
        setattr(SolutionArray, name, state2_prop(name))

    # Factory for creating properties which consist of a tuple of three
    # variables, e.g. 'TPY' or 'UVX'
    def state3_prop(name):
        def getter(self):
            a = np.empty(self._shape)
            b = np.empty(self._shape)
            c = np.empty(self._shape + (self._phase.n_species,))
            for index in self._indices:
                self._phase.state = self._states[index]
                a[index], b[index], c[index] = getattr(self._phase, name)
            return a, b, c

        def setter(self, ABC):
            assert len(ABC) == 3, "Expected 3 elements, got {}".format(len(ABC))
            A, B, C, _ = np.broadcast_arrays(ABC[0], ABC[1], ABC[2],
                                             self._states[...,0])
            for index in self._indices:
                self._phase.state = self._states[index]
                setattr(self._phase, name, (A[index], B[index], C[index]))
                self._states[index][:] = self._phase.state

        return property(getter, setter, doc=getattr(Solution, name).__doc__)

    for name in state3:
        setattr(SolutionArray, name, state3_prop(name))

    def scalar_prop(name):
        def getter(self):
            v = np.empty(self._shape)
            for index in self._indices:
                self._phase.state = self._states[index]
                v[index] = getattr(self._phase, name)
            return v
        return property(getter, doc=getattr(Solution, name).__doc__)

    for name in scalar:
        setattr(SolutionArray, name, scalar_prop(name))

    def species_prop(name):
        def getter(self):
            v = np.empty(self._shape + (self._phase.n_species,))
            for index in self._indices:
                self._phase.state = self._states[index]
                v[index] = getattr(self._phase, name)
            return v
        return property(getter, doc=getattr(Solution, name).__doc__)

    for name in n_species:
        setattr(SolutionArray, name, species_prop(name))

    def species2_prop(name):
        def getter(self):
            v = np.empty(self._shape +
                         (self._phase.n_species,self._phase.n_species))
            for index in self._indices:
                self._phase.state = self._states[index]
                v[index] = getattr(self._phase, name)
            return v
        return property(getter, doc=getattr(Solution, name).__doc__)

    for name in n_species2:
        setattr(SolutionArray, name, species2_prop(name))

    def reaction_prop(name):
        def getter(self):
            v = np.empty(self._shape + (self._phase.n_reactions,))
            for index in self._indices:
                self._phase.state = self._states[index]
                v[index] = getattr(self._phase, name)
            return v
        return property(getter, doc=getattr(Solution, name).__doc__)

    for name in n_reactions:
        setattr(SolutionArray, name, reaction_prop(name))

    # Factory for creating wrappers for functions which return a value
    def caller(name):
        def wrapper(self, *args, **kwargs):
            v = np.empty(self._shape)
            for index in self._indices:
                self._phase.state = self._states[index]
                v[index] = getattr(self._phase, name)(*args, **kwargs)
            return v
        return wrapper

    for name in call:
        setattr(SolutionArray, name, caller(name))

    # Factory for creating properties to pass through state-independent
    # functions and properties unmodified. Having a setter is ok even for read-
    # only properties, since the wrapped class will just raise an exception
    def passthrough_prop(name):
        def getter(self):
            return getattr(self._phase, name)

        def setter(self, value):
            setattr(self._phase, name, value)

        return property(getter, setter, doc=getattr(Solution, name).__doc__)

    for name in passthrough:
        setattr(SolutionArray, name, passthrough_prop(name))

_make_functions()
