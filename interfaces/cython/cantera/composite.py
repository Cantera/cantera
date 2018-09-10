# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

from ._cantera import *
import numpy as np
import csv as _csv

class Solution(ThermoPhase, Kinetics, Transport):
    """
    A class for chemically-reacting solutions. Instances can be created to
    represent any type of solution -- a mixture of gases, a liquid solution, or
    a solid solution, for example.

    Class `Solution` derives from classes `ThermoPhase`, `Kinetics`, and
    `Transport`.  It defines no methods of its own, and is provided so that a
    single object can be used to compute thermodynamic, kinetic, and transport
    properties of a solution.

    To skip initialization of the Transport object, pass the keyword argument
    ``transport_model=None`` to the `Solution` constructor.

    The most common way to instantiate `Solution` objects is by using a phase
    definition, species and reactions defined in an input file::

        gas = ct.Solution('gri30.cti')

    If an input file defines multiple phases, the phase *name* (in CTI) or *id*
    (in XML) can be used to specify the desired phase::

        gas = ct.Solution('diamond.cti', 'gas')
        diamond = ct.Solution('diamond.cti', 'diamond')

    `Solution` objects can also be constructed using `Species` and `Reaction`
    objects which can themselves either be imported from input files or defined
    directly in Python::

        spec = ct.Species.listFromFile('gri30.cti')
        rxns = ct.Reaction.listFromFile('gri30.cti')
        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=spec, reactions=rxns)

    where the ``thermo`` and ``kinetics`` keyword arguments are strings
    specifying the thermodynamic and kinetics model, respectively, and
    ``species`` and ``reactions`` keyword arguments are lists of `Species` and
    `Reaction` objects, respectively.

    For non-trivial uses cases of this functionality, see the examples
    `extract_submechanism.py <https://cantera.org/examples/python/extract_submechanism.py.html>`_
    and `mechanism_reduction.py <https://cantera.org/examples/python/mechanism_reduction.py.html>`_.

    In addition, `Solution` objects can be constructed by passing the text of
    the CTI or XML phase definition in directly, using the ``source`` keyword
    argument::

        cti_def = '''
            ideal_gas(name='gas', elements='O H Ar',
                      species='gri30: all',
                      reactions='gri30: all',
                      options=['skip_undeclared_elements', 'skip_undeclared_species', 'skip_undeclared_third_bodies'],
                      initial_state=state(temperature=300, pressure=101325))'''
        gas = ct.Solution(source=cti_def)
    """
    __slots__ = ()


class Interface(InterfacePhase, InterfaceKinetics):
    """
    Two-dimensional interfaces.

    Instances of class `Interface` represent reacting 2D interfaces between bulk
    3D phases. Class `Interface` defines no methods of its own. All of its
    methods derive from either `InterfacePhase` or `InterfaceKinetics`.

    To construct an `Interface` object, adjacent bulk phases which participate
    in reactions need to be created and then passed in as a list in the
    ``phases`` argument to the constructor::

        gas = ct.Solution('diamond.cti', 'gas')
        diamond = ct.Solution('diamond.cti', 'diamond')
        diamond_surf = ct.Interface('diamond.cti', 'diamond_100', [gas, diamond])
    """
    __slots__ = ('_phase_indices',)


class DustyGas(ThermoPhase, Kinetics, DustyGasTransport):
    """
    A composite class which models a gas in a stationary, solid, porous medium.

    The only transport properties computed are the multicomponent diffusion
    coefficients. The model does not compute viscosity or thermal conductivity.

    """
    __slots__ = ()


class Quantity:
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
        return Quantity(self.phase, mass=self.mass * other, constant=self.constant)

    def __rmul__(self, other):
        return Quantity(self.phase, mass=self.mass * other, constant=self.constant)

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


class SolutionArray:
    """
    A class providing a convenient interface for representing many thermodynamic
    states using the same `Solution` object and computing properties for that
    array of states.

    SolutionArray can represent both 1D and multi-dimensional arrays of states,
    with shapes described in the same way as Numpy arrays. All of the states
    can be set in a single call.

        >>> gas = ct.Solution('gri30.cti')
        >>> states = ct.SolutionArray(gas, (6, 10))
        >>> T = np.linspace(300, 1000, 10) # row vector
        >>> P = ct.one_atm * np.linspace(0.1, 5.0, 6)[:,np.newaxis] # column vector
        >>> X = 'CH4:1.0, O2:1.0, N2:3.76'
        >>> states.TPX = T, P, X

    Similar to Numpy arrays, input with fewer non-singleton dimensions than the
    SolutionArray is 'broadcast' to generate input of the appropriate shape. In
    the above example, the single value for the mole fraction input is applied
    to each input, while each row has a constant temperature and each column has
    a constant pressure.

    Computed properties are returned as Numpy arrays with the same shape as the
    array of states, with additional dimensions appended as necessary for non-
    scalar output (e.g. per-species or per-reaction properties)::

        >>> h = states.enthalpy_mass
        >>> h[i,j] # -> enthalpy at P[i] and T[j]
        >>> sk = states.partial_molar_entropies
        >>> sk[i,:,k] # -> entropy of species k at P[i] and each temperature
        >>> ropnet = states.net_rates_of_progress
        >>> ropnet[i,j,n] # -> net reaction rate for reaction n at P[i] and T[j]

    In the case of 1D arrays, additional states can be appended one at a time::

        >>> states = ct.SolutionArray(gas) # creates an empty SolutionArray
        >>> for phi in np.linspace(0.5, 2.0, 20):
        ...     states.append(T=300, P=ct.one_atm,
        ...                   X={'CH4': phi, 'O2': 2, 'N2': 2*3.76})
        >>> # 'states' now contains 20 elements
        >>> states.equilibrate('HP')
        >>> states.T # -> adiabatic flame temperature at various equivalence ratios

    SolutionArray objects can also be 'sliced' like Numpy arrays, which can be
    used both for accessing and setting properties::

        >>> states = ct.SolutionArray(gas, (6, 10))
        >>> states[0].TP = 400, None # set the temperature of the first row to 400 K
        >>> cp = states[:,1].cp_mass # heat capacity of the second column

    If many slices or elements of a property are going to be accessed (i.e.
    within a loop), it is generally more efficient to compute the property array
    once and access this directly, rather than repeatedly slicing the
    SolutionArray object, e.g.::

        >>> mu = states.viscosity
        >>> for i,j in np.ndindex(mu.shape):
        ...     # do something with mu[i,j]

    Information about a subset of species may also be accessed, using
    parentheses to specify the species:

        >>> states('CH4').Y # -> mass fraction of CH4 in each state
        >>> states('H2','O2').partial_molar_cp # -> cp for H2 and O2

    Properties and functions which are not dependent on the thermodynamic state
    act as pass-throughs to the underlying `Solution` object, and are not
    converted to arrays::

        >>> states.element_names
        ['O', 'H', 'C', 'N', 'Ar']
        >>> s.reaction_equation(10)
        'CH4 + O <=> CH3 + OH'

    Data represnted by a SolutionArray can be extracted and saved to a CSV file
    using the `write_csv` method:

        >>> states.write_csv('somefile.csv', cols=('T','P','X','net_rates_of_progress'))

    :param phase: The `Solution` object used to compute the thermodynamic,
        kinetic, and transport properties
    :param shape: A tuple or integer indicating the dimensions of the
        SolutionArray. If the shape is 1D, the array may be extended using the
        `append` method. Otherwise, the shape is fixed.
    :param states: The initial array of states. Used internally to provide
        slicing support.
    """

    _scalar = [
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
    _n_species = [
        # from ThermoPhase
        'Y', 'X', 'concentrations', 'partial_molar_enthalpies',
        'partial_molar_entropies', 'partial_molar_int_energies',
        'chemical_potentials', 'electrochemical_potentials', 'partial_molar_cp',
        'partial_molar_volumes', 'standard_enthalpies_RT',
        'standard_entropies_R', 'standard_int_energies_RT', 'standard_gibbs_RT',
        'standard_cp_R',
        # From Transport
        'mix_diff_coeffs', 'mix_diff_coeffs_mass', 'mix_diff_coeffs_mole',
        'thermal_diff_coeffs'
    ]

    # From Kinetics (differs from Solution.n_species for Interface phases)
    _n_total_species = [
        'creation_rates', 'destruction_rates', 'net_production_rates',
    ]

    _n_species2 = ['multi_diff_coeffs', 'binary_diff_coeffs']

    _n_reactions = [
        'forward_rates_of_progress', 'reverse_rates_of_progress',
        'net_rates_of_progress', 'equilibrium_constants',
        'forward_rate_constants', 'reverse_rate_constants',
        'delta_enthalpy', 'delta_gibbs', 'delta_entropy',
        'delta_standard_enthalpy', 'delta_standard_gibbs',
        'delta_standard_entropy'
    ]
    _state2 = ['TD', 'TP', 'UV', 'DP', 'HP', 'SP', 'SV']
    _call_scalar = ['elemental_mass_fraction', 'elemental_mole_fraction']

    _passthrough = [
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

    _interface_passthrough = ['site_density']
    _interface_n_species = ['coverages']

    def __init__(self, phase, shape=(0,), states=None, extra=None):
        self._phase = phase

        if isinstance(shape, int):
            shape = (shape,)

        if states is not None:
            self._shape = np.shape(states)[:-1]
            self._states = states
        else:
            self._shape = tuple(shape)
            if len(shape) == 1:
                S = [self._phase.state for _ in range(shape[0])]
            else:
                S = np.empty(shape + (2+self._phase.n_species,))
                S[:] = self._phase.state
            self._states = S

        if len(self._shape) == 1:
            self._indices = list(range(self._shape[0]))
            self._output_dummy = self._indices
        else:
            self._indices = list(np.ndindex(self._shape))
            self._output_dummy = self._states[..., 0]

        self._extra_lists = {}
        self._extra_arrays = {}
        if isinstance(extra, dict):
            for name, v in extra.items():
                if not np.shape(v):
                    self._extra_lists[name] = [v]*self._shape[0]
                    self._extra_arrays[name] = np.array(self._extra_lists[name])
                elif len(v) == self._shape[0]:
                    self._extra_lists[name] = list(v)
                else:
                    raise ValueError("Unable to map extra SolutionArray"
                                     "input for named {!r}".format(name))
                self._extra_arrays[name] = np.array(self._extra_lists[name])

        elif extra and self._shape == (0,):
            for name in extra:
                self._extra_lists[name] = []
                self._extra_arrays[name] = np.array(())

        elif extra:
            raise ValueError("Initial values for extra properties must be"
                " supplied in a dict if the SolutionArray is not initially"
                " empty")

    def __getitem__(self, index):
        states = self._states[index]
        shape = states.shape[:-1]
        return SolutionArray(self._phase, shape, states)

    def __getattr__(self, name):
        if name not in self._extra_lists:
            raise AttributeError("'{}' object has no attribute '{}'".format(
                self.__class__.__name__, name))
        L = self._extra_lists[name]
        A = self._extra_arrays[name]
        if len(L) != len(A):
            A = np.array(L)
            self._extra_arrays[name] = A
        return A

    def __call__(self, *species):
        return SolutionArray(self._phase[species], states=self._states,
                             extra=self._extra_lists)

    def append(self, state=None, **kwargs):
        """
        Append an element to the array with the specified state. Elements can
        only be appended in cases where the array of states is one-dimensional.

        The state may be specified in one of three ways:

        - as the array of [temperature, density, mass fractions] which is
          returned by `Solution.state`::

              mystates.append(gas.state)

        - as a tuple of three elements that corresponds to any of the full-state
          setters of `Solution`, e.g. `TPY` or `HPX`::

              mystates.append(TPX=(300, 101325, 'O2:1.0, N2:3.76'))

        - as separate keywords for each of the elements corresponding to one of
          the full-state setters::

              mystates.append(T=300, P=101325, X={'O2':1.0, 'N2':3.76})
        """
        if len(self._shape) != 1:
            raise IndexError("Can only append to 1D SolutionArray")

        for name, value in self._extra_lists.items():
            value.append(kwargs.pop(name))

        if state is not None:
            self._phase.state = state

        elif len(kwargs) == 1:
            attr, value = next(iter(kwargs.items()))
            if frozenset(attr) not in self._phase._full_states:
                raise KeyError("{} does not specify a full thermodynamic state")
            setattr(self._phase, attr, value)

        else:
            try:
                attr = self._phase._full_states[frozenset(kwargs)]
            except KeyError:
                raise KeyError("{} is not a valid combination of properties "
                    "for setting the thermodynamic state".format(tuple(kwargs)))
            setattr(self._phase, attr, [kwargs[a] for a in attr])

        self._states.append(self._phase.state)
        self._indices.append(len(self._indices))
        self._shape = (len(self._indices),)

    def equilibrate(self, *args, **kwargs):
        """ See `ThermoPhase.equilibrate` """
        for index in self._indices:
            self._phase.state = self._states[index]
            self._phase.equilibrate(*args, **kwargs)
            self._states[index][:] = self._phase.state

    def collect_data(self, cols=('extra','T','density','Y'), threshold=0,
                     species='Y'):
        """
        Returns the data specified by *cols* in a single 2D Numpy array, along
        with a list of column labels.

        :param cols: A list of any properties of the solution that are scalars
            or which have a value for each species or reaction. If species names
            are specified, then either the mass or mole fraction of that species
            will be taken, depending on the value of *species*. *cols* may also
            include any arrays which were specified as 'extra' variables when
            defining the SolutionArray object. The special value 'extra' can be
            used to include all 'extra' variables.
        :param threshold: Relative tolerance for including a particular column.
            The tolerance is applied by comparing the maximum absolute value for
            a particular column to the maximum absolute value in all columns for
            the same variable (e.g. mass fraction).
        :param species: Specifies whether to use mass ('Y') or mole ('X')
            fractions for individual species specified in 'cols'
        """
        if len(self._shape) != 1:
            raise TypeError("collect_data only works for 1D SolutionArray")
        data = []
        labels = []

        # Expand cols to include the individual items in 'extra'
        expanded_cols = []
        for c in cols:
            if c == 'extra':
                expanded_cols.extend(self._extra_arrays)
            else:
                expanded_cols.append(c)

        species_names = set(self.species_names)
        for c in expanded_cols:
            single_species = False
            # Determine labels for the items in the current group of columns
            if c in self._extra_arrays:
                collabels = [c]
            elif c in self._scalar:
                collabels = [c]
            elif c in self._n_species:
                collabels = ['{}_{}'.format(c, s) for s in self.species_names]
            elif c in self._n_reactions:
                collabels = ['{} {}'.format(c, r)
                             for r in self.reaction_equations()]
            elif c in species_names:
                single_species = True
                collabels = ['{}_{}'.format(species, c)]
            else:
                raise CanteraError('property "{}" not supported'.format(c))

            # Get the data for the current group of columns
            if single_species:
                d = getattr(self(c), species)
            else:
                d = getattr(self, c)

            if d.ndim == 1:
                d = d[:, np.newaxis]
            elif threshold:
                # Determine threshold value and select columns to keep
                maxval = abs(d).max()
                keep = (abs(d) > threshold * maxval).any(axis=0)
                d = d[:, keep]
                collabels = [label for label, k in zip(collabels, keep) if k]

            data.append(d)
            labels.extend(collabels)

        return np.hstack(data), labels

    def write_csv(self, filename, cols=('extra','T','density','Y'),
                  *args, **kwargs):
        """
        Write a CSV file named *filename* containing the data specified by
        *cols*. The first row of the CSV file will contain column labels.

        Additional arguments are passed on to `collect_data`. This method works
        only with 1D SolutionArray objects.
        """
        data, labels = self.collect_data(cols, *args, **kwargs)
        with open(filename, 'w') as outfile:
            writer = _csv.writer(outfile)
            writer.writerow(labels)
            for row in data:
                writer.writerow(row)


def _make_functions():
    # this is wrapped in a function to avoid polluting the module namespace

    # Factory for creating properties which consist of a tuple of two variables,
    # e.g. 'TP' or 'SV'
    def state2_prop(name, doc_source):
        def getter(self):
            a = np.empty(self._shape)
            b = np.empty(self._shape)
            for index in self._indices:
                self._phase.state = self._states[index]
                a[index], b[index] = getattr(self._phase, name)
            return a, b

        def setter(self, AB):
            assert len(AB) == 2, "Expected 2 elements, got {}".format(len(AB))
            A, B, _ = np.broadcast_arrays(AB[0], AB[1], self._output_dummy)
            for index in self._indices:
                self._phase.state = self._states[index]
                setattr(self._phase, name, (A[index], B[index]))
                self._states[index][:] = self._phase.state

        return property(getter, setter, doc=getattr(doc_source, name).__doc__)

    for name in SolutionArray._state2:
        setattr(SolutionArray, name, state2_prop(name, Solution))

    for name in PureFluid._full_states.values():
        setattr(SolutionArray, name, state2_prop(name, PureFluid))

    # Factory for creating properties which consist of a tuple of three
    # variables, e.g. 'TPY' or 'UVX'
    def state3_prop(name):
        def getter(self):
            a = np.empty(self._shape)
            b = np.empty(self._shape)
            c = np.empty(self._shape + (self._phase.n_selected_species,))
            for index in self._indices:
                self._phase.state = self._states[index]
                a[index], b[index], c[index] = getattr(self._phase, name)
            return a, b, c

        def setter(self, ABC):
            assert len(ABC) == 3, "Expected 3 elements, got {}".format(len(ABC))
            A, B, _ = np.broadcast_arrays(ABC[0], ABC[1], self._output_dummy)
            XY = ABC[2] # composition
            if len(np.shape(XY)) < 2:
                # composition is a single array (or string or dict)
                for index in self._indices:
                    self._phase.state = self._states[index]
                    setattr(self._phase, name, (A[index], B[index], XY))
                    self._states[index][:] = self._phase.state
            else:
                # composition is an array with trailing dimension n_species
                C = np.empty(self._shape + (self._phase.n_selected_species,))
                C[:] = XY
                for index in self._indices:
                    self._phase.state = self._states[index]
                    setattr(self._phase, name, (A[index], B[index], C[index]))
                    self._states[index][:] = self._phase.state

        return property(getter, setter, doc=getattr(Solution, name).__doc__)

    for name in Solution._full_states.values():
        setattr(SolutionArray, name, state3_prop(name))

    # Functions which define empty output arrays of an appropriate size for
    # different properties
    def empty_scalar(self):
        return np.empty(self._shape)

    def empty_species(self):
        return np.empty(self._shape + (self._phase.n_selected_species,))

    def empty_total_species(self):
        return np.empty(self._shape + (self._phase.n_total_species,))

    def empty_species2(self):
        return np.empty(self._shape + (self._phase.n_species,
                                       self._phase.n_species))

    def empty_reactions(self):
        return np.empty(self._shape + (self._phase.n_reactions,))

    # Factory for creating read-only properties
    def make_prop(name, get_container, doc_source):
        def getter(self):
            v = get_container(self)
            for index in self._indices:
                self._phase.state = self._states[index]
                v[index] = getattr(self._phase, name)
            return v
        return property(getter, doc=getattr(doc_source, name).__doc__)

    for name in SolutionArray._scalar:
        setattr(SolutionArray, name, make_prop(name, empty_scalar, Solution))

    for name in SolutionArray._n_species:
        setattr(SolutionArray, name, make_prop(name, empty_species, Solution))

    for name in SolutionArray._interface_n_species:
        setattr(SolutionArray, name, make_prop(name, empty_species, Interface))

    for name in SolutionArray._n_total_species:
        setattr(SolutionArray, name,
                make_prop(name, empty_total_species, Solution))

    for name in SolutionArray._n_species2:
        setattr(SolutionArray, name, make_prop(name, empty_species2, Solution))

    for name in SolutionArray._n_reactions:
        setattr(SolutionArray, name, make_prop(name, empty_reactions, Solution))

    # Factory for creating wrappers for functions which return a value
    def caller(name, get_container):
        def wrapper(self, *args, **kwargs):
            v = get_container(self)
            for index in self._indices:
                self._phase.state = self._states[index]
                v[index] = getattr(self._phase, name)(*args, **kwargs)
            return v
        return wrapper

    for name in SolutionArray._call_scalar:
        setattr(SolutionArray, name, caller(name, empty_scalar))

    # Factory for creating properties to pass through state-independent
    # functions and properties unmodified. Having a setter is ok even for read-
    # only properties, since the wrapped class will just raise an exception
    def passthrough_prop(name, doc_source):
        def getter(self):
            return getattr(self._phase, name)

        def setter(self, value):
            setattr(self._phase, name, value)

        return property(getter, setter, doc=getattr(doc_source, name).__doc__)

    for name in SolutionArray._passthrough:
        setattr(SolutionArray, name, passthrough_prop(name, Solution))

    for name in SolutionArray._interface_passthrough:
        setattr(SolutionArray, name, passthrough_prop(name, Interface))

_make_functions()
