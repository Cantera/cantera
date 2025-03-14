# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.
from __future__ import annotations

from ._cantera import *
import numpy as np
import csv as _csv
import importlib.metadata
import warnings


_pandas = None
def _import_pandas():
    # defer import of pandas
    global _pandas
    if _pandas is not None:
        return
    try:
        importlib.metadata.version('pandas')
    except importlib.metadata.PackageNotFoundError:
        raise ImportError('Method requires a working pandas installation.')
    else:
        import pandas as _pandas


class Solution(Transport, Kinetics, ThermoPhase):
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

        gas = ct.Solution('gri30.yaml')

    If an input file defines multiple phases, the corresponding key in the
    ``phases`` map can be used to specify the desired phase via the ``name`` keyword
    argument of the constructor::

        gas = ct.Solution('diamond.yaml', name='gas')
        diamond = ct.Solution('diamond.yaml', name='diamond')

    The name of the `Solution` object defaults to the *phase* identifier
    specified in the input file. Upon initialization of a `Solution` object,
    a custom name can assigned via::

        gas.name = 'my_custom_name'

    `Solution` objects can also be constructed using `Species` and `Reaction`
    objects which can themselves either be imported from input files or defined
    directly in Python::

        spec = ct.Species.list_from_file("gri30.yaml")
        spec_gas = ct.Solution(thermo='ideal-gas', species=spec)
        rxns = ct.Reaction.list_from_file("gri30.yaml", spec_gas)
        gas = ct.Solution(thermo='ideal-gas', kinetics='gas',
                          species=spec, reactions=rxns, name='my_custom_name')

    where the ``thermo`` and ``kinetics`` keyword arguments are strings
    specifying the thermodynamic and kinetics model, respectively, and
    ``species`` and ``reactions`` keyword arguments are lists of `Species` and
    `Reaction` objects, respectively. Note that importing the reactions from a
    YAML input file requires a `Kinetics` object containing the species, as
    shown.

    Types of underlying models that form the composite `Solution` object are
    queried using the ``thermo_model``, ``kinetics_model`` and
    ``transport_model`` attributes; further, the ``composite`` attribute is a
    shorthand returning a tuple containing the types of the three constitutive
    models.

    For non-trivial uses cases of this functionality, see the examples
    :doc:`extract_submechanism.py </examples/python/kinetics/extract_submechanism>`
    :doc:`mechanism_reduction.py </examples/python/kinetics/mechanism_reduction>`.

    In addition, `Solution` objects can be constructed by passing the text of
    the YAML phase definition in directly, using the ``yaml`` keyword
    argument::

        yaml_def = '''
        phases:
        - name: gas
          thermo: ideal-gas
          kinetics: gas
          elements: [O, H, Ar]
          species:
          - gri30.yaml/species: all
          reactions:
          - gri30.yaml/reactions: declared-species
          skip-undeclared-elements: true
          skip-undeclared-third-bodies: true
          state: {T: 300, P: 1 atm}
        '''
        gas = ct.Solution(yaml=yaml_def)
    """
    __slots__ = ()


class Interface(InterfaceKinetics, InterfacePhase):
    """
    Instances of class `Interface` represent reacting 2D surfaces between bulk 3D
    phases, or 1D edges where multiple surfaces (and bulk phases) meet. Class
    `Interface` defines no methods of its own. All of its methods derive from either
    `InterfacePhase` or `InterfaceKinetics`.

    Constructing an `Interface` object also involves constructing adjacent bulk phases
    that participate in reactions. This is done automatically if the adjacent phases
    are specified as part of the ``adjacent-phases`` entry in the YAML phase definition::

        diamond_surf = ct.Interface("diamond.yaml", name="diamond_100")
        gas = diamond_surf.adjacent["gas"]
        diamond = diamond_surf.adjacent["diamond"]

    This behavior can be overridden by specifying the adjacent phases explicitly,
    either using their name in the input file, or by constructing corresponding
    `Solution` objects::

        gas = ct.Solution("diamond.yaml", name="gas")
        diamond = ct.Solution("diamond.yaml", name="diamond")
        diamond_surf = ct.Interface("diamond.yaml", name="diamond_100",
                                    adjacent=[gas, diamond])
    """
    __slots__ = ('_phase_indices',)


class DustyGas(DustyGasTransport, Kinetics, ThermoPhase):
    """
    A composite class which models a gas in a stationary, solid, porous medium.

    The only transport properties computed are the multicomponent diffusion
    coefficients. The model does not compute viscosity or thermal conductivity.
    """
    __slots__ = ()


# A pure-Python class to store weakrefs to
class _WeakrefProxy:
    pass


class Quantity:
    """
    A class representing a specific quantity of a `Solution`. In addition to the
    properties which can be computed for class `Solution`, class `Quantity`
    provides several additional capabilities. A `Quantity` object is created
    from a `Solution` with either the mass or number of moles specified::

        >>> gas = ct.Solution('gri30.yaml')
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
    __slots__ = ("state", "_phase", "_id", "mass", "constant", "_weakref_proxy")

    def __init__(self, phase, mass=None, moles=None, constant='UV'):
        self.state = phase.state
        self._phase = phase
        self._weakref_proxy = _WeakrefProxy()
        phase._references[self._weakref_proxy] = True

        # A unique key to prevent adding phases with different species
        # definitions
        self._id = hash((phase.name,) + tuple(phase.species_names))

        if mass is not None:
            self.mass = mass
        elif moles is not None:
            self.moles = moles
        else:
            self.mass = 1.0

        if constant not in ('TP','TV','HP','SP','SV','UV'):
            raise ValueError(
                f"Constant {constant} is invalid. "
                "Must be one of 'TP','TV','HP','SP','SV', or 'UV'")
        self.constant = constant

    @property
    def phase(self):
        """
        Get the underlying `Solution` object, with the state set to match the
        wrapping `Quantity` object.
        """
        self._phase.state = self.state
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
        self.state = self._phase.state

    def set_equivalence_ratio(self, phi, fuel, oxidizer, basis="mole", *, diluent=None,
                              fraction=None):
        self._phase.state = self.state
        self._phase.set_equivalence_ratio(phi, fuel, oxidizer, basis, diluent=diluent,
                                          fraction=fraction)
        self.state = self._phase.state
    set_equivalence_ratio.__doc__ = Solution.set_equivalence_ratio.__doc__

    def set_mixture_fraction(self, mixture_fraction, fuel, oxidizer, basis='mole'):
        self._phase.state = self.state
        self._phase.set_mixture_fraction(mixture_fraction, fuel, oxidizer, basis)
        self.state = self._phase.state
    set_mixture_fraction.__doc__ = Solution.set_mixture_fraction.__doc__

    def __imul__(self, other):
        self.mass *= other
        return self

    def __mul__(self, other):
        return Quantity(self.phase, mass=self.mass * other, constant=self.constant)

    def __rmul__(self, other):
        return Quantity(self.phase, mass=self.mass * other, constant=self.constant)

    def __iadd__(self, other):
        if self._id != other._id:
            raise ValueError(
                'Cannot add Quantities with different phase '
                f'definitions. {self._id} != {other._id}')
        if self.constant != other.constant:
            raise ValueError(
                "Cannot add Quantities with different "
                f"constant values. {self.constant} != {other.constant}")

        m = self.mass + other.mass
        Y = (self.Y * self.mass + other.Y * other.mass)
        if self.constant == 'UV':
            U = self.int_energy + other.int_energy
            V = self.volume + other.volume
            if self.basis == 'mass':
                self._phase.UVY = U / m, V / m, Y
            else:
                n = self.moles + other.moles
                self._phase.UVY = U / n, V / n, Y
        else:  # self.constant == 'HP'
            dp_rel = 2 * abs(self.P - other.P) / (self.P + other.P)
            if dp_rel > 1.0e-7:
                raise ValueError(
                    'Cannot add Quantities at constant pressure when '
                    f'pressure is not equal ({self.P} != {other.P})')

            H = self.enthalpy + other.enthalpy
            if self.basis == 'mass':
                self._phase.HPY = H / m, None, Y
            else:
                n = self.moles + other.moles
                self._phase.HPY = H / n, None, Y

        self.state = self._phase.state
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
        self.state = self._phase.state

    return property(getter, setter, doc=getattr(Solution, attr).__doc__)

def _method(orig):
    def f(self, *args, **kwargs):
        return orig(self.phase, *args, **kwargs)
    f.__doc__ = orig.__doc__
    return f

for _attr in dir(Solution):
    if _attr.startswith('_') or _attr in Quantity.__dict__ or _attr == 'state':
        continue
    if _attr.startswith('set_unnormalized'):
        continue
    else:
        _orig = getattr(Solution, _attr)
        if hasattr(_orig, "__call__"):
            setattr(Quantity, _attr, _method(_orig))
        else:
            setattr(Quantity, _attr, _prop(_attr))


class SolutionArray(SolutionArrayBase):
    """
    A class providing a convenient interface for representing many thermodynamic
    states using the same `Solution` object and computing properties for that
    array of states.

    `SolutionArray` can represent both 1D and multi-dimensional arrays of states,
    with shapes described in the same way as *NumPy* arrays. All of the states
    can be set in a single call::

        >>> gas = ct.Solution('gri30.yaml')
        >>> states = ct.SolutionArray(gas, (6, 10))
        >>> T = np.linspace(300, 1000, 10) # row vector
        >>> P = ct.one_atm * np.linspace(0.1, 5.0, 6)[:,np.newaxis] # column vector
        >>> X = 'CH4:1.0, O2:1.0, N2:3.76'
        >>> states.TPX = T, P, X

    Similar to *NumPy* arrays, input with fewer non-singleton dimensions than the
    `SolutionArray` is 'broadcast' to generate input of the appropriate shape. In
    the above example, the single value for the mole fraction input is applied
    to each input, while each row has a constant temperature and each column has
    a constant pressure.

    Computed properties are returned as *NumPy* arrays with the same shape as the
    array of states, with additional dimensions appended as necessary for non-
    scalar output (for example, per-species or per-reaction properties)::

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

    `SolutionArray` objects can also be 'sliced' like *NumPy* arrays, which can be
    used both for accessing and setting properties::

        >>> states = ct.SolutionArray(gas, (6, 10))
        >>> states[0].TP = 400, None # set the temperature of the first row to 400 K
        >>> cp = states[:,1].cp_mass # heat capacity of the second column

    If many slices or elements of a property are going to be accessed (for example,
    within a loop), it is generally more efficient to compute the property array
    once and access this directly, rather than repeatedly slicing the
    `SolutionArray` object, for example::

        >>> mu = states.viscosity
        >>> for i,j in np.ndindex(mu.shape):
        ...     # do something with mu[i,j]

    Information about a subset of species may also be accessed, using
    parentheses to specify the species::

        >>> states('CH4').Y # -> mass fraction of CH4 in each state
        >>> states('H2','O2').partial_molar_cp # -> cp for H2 and O2

    Properties and functions which are not dependent on the thermodynamic state
    act as pass-throughs to the underlying `Solution` object, and are not
    converted to arrays::

        >>> states.element_names
        ['O', 'H', 'C', 'N', 'Ar']
        >>> s.reaction_equation(10)
        'CH4 + O <=> CH3 + OH'

    Data represented by a `SolutionArray` can be extracted to a CSV file using the
    `save` method::

        >>> states.save('somefile.csv', basis="mole')

    As an alternative to the CSV format, `SolutionArray` objects can also be saved in
    YAML or HDF formats, where the keyword argument ``id`` allows for saving and
    accessing of multiple solutions in a single container file::

        >>> states.save('somefile.yaml', id='some_key')

    YAML and HDF files can be read back into `SolutionArray` objects using the
    `restore` method::

        >>> states = ct.SolutionArray(gas)
        >>> states.restore('somefile.yaml', id='some_key')

    As long as stored columns in a CSV file specify a valid thermodynamic state, the
    contents of a `SolutionArray` can be restored using the `read_csv` method, which
    is specific to the Python API::

        >>> states = ct.SolutionArray(gas)
        >>> states.read_csv('somefile.csv')

    Note that `save` and `restore` for HDF requires Cantera to be compiled with HDF
    support, as it depends on external *HighFive* and *HDF5* libraries.

    :param phase: The `Solution` object used to compute the thermodynamic,
        kinetic, and transport properties
    :param shape: A tuple or integer indicating the dimensions of the
        `SolutionArray`. If the shape is 1D, the array may be extended using the
        `append` method. Otherwise, the shape is fixed.
    :param states: The initial array of states. Used internally to provide
        slicing support.
    """

    _scalar = [
        # From ThermoPhase
        'mean_molecular_weight', 'P', 'T', 'Te', 'density', 'density_mass',
        'density_mole', 'v', 'volume_mass', 'volume_mole', 'u',
        'int_energy_mole', 'int_energy_mass', 'h', 'enthalpy_mole',
        'enthalpy_mass', 's', 'entropy_mole', 'entropy_mass', 'g', 'gibbs_mole',
        'gibbs_mass', 'cv', 'cv_mole', 'cv_mass', 'cp', 'cp_mole', 'cp_mass',
        'critical_temperature', 'critical_pressure', 'critical_density',
        'P_sat', 'T_sat', 'isothermal_compressibility',
        'thermal_expansion_coeff', 'sound_speed', 'electric_potential',
        # From Kinetics
        'heat_release_rate',
        # From Transport
        'viscosity', 'electrical_conductivity', 'thermal_conductivity',
    ]
    _strings = ['phase_of_matter']
    _n_species = [
        # from ThermoPhase
        'Y', 'X', 'concentrations', 'partial_molar_enthalpies',
        'partial_molar_entropies', 'partial_molar_int_energies',
        'chemical_potentials', 'electrochemical_potentials', 'partial_molar_cp',
        'partial_molar_volumes', 'standard_enthalpies_RT',
        'standard_entropies_R', 'standard_int_energies_RT', 'standard_gibbs_RT',
        'standard_cp_R', 'activities', 'activity_coefficients',
        # From Transport
        'mix_diff_coeffs', 'mix_diff_coeffs_mass', 'mix_diff_coeffs_mole',
        'thermal_diff_coeffs', 'mobilities', 'species_viscosities',
    ]

    # From Kinetics (differs from Solution.n_species for Interface phases)
    _n_total_species = [
        'creation_rates', 'destruction_rates', 'net_production_rates',
        'creation_rates_ddC', 'creation_rates_ddP', 'creation_rates_ddT',
        'destruction_rates_ddC', 'destruction_rates_ddP', 'destruction_rates_ddT',
        'net_production_rates_ddC', 'net_production_rates_ddP',
        'net_production_rates_ddT'
    ]

    _n_species2 = [
        'multi_diff_coeffs', 'binary_diff_coeffs', 'creation_rates_ddX',
        'destruction_rates_ddX', 'net_production_rates_ddX', 'creation_rates_ddCi',
        'destruction_rates_ddCi', 'net_production_rates_ddCi'
    ]

    _n_reactions = [
        'forward_rates_of_progress', 'reverse_rates_of_progress',
        'net_rates_of_progress', 'equilibrium_constants',
        'forward_rate_constants', 'reverse_rate_constants',
        'delta_enthalpy', 'delta_gibbs', 'delta_entropy',
        'delta_standard_enthalpy', 'delta_standard_gibbs',
        'delta_standard_entropy', 'heat_production_rates',
        'forward_rate_constants_ddC', 'forward_rate_constants_ddP',
        'forward_rate_constants_ddT', 'forward_rates_of_progress_ddC',
        'forward_rates_of_progress_ddP', 'forward_rates_of_progress_ddT',
        'net_rates_of_progress_ddC', 'net_rates_of_progress_ddP',
        'net_rates_of_progress_ddT', 'reverse_rates_of_progress_ddC',
        'reverse_rates_of_progress_ddP', 'reverse_rates_of_progress_ddP',
        'reverse_rates_of_progress_ddT', 'third_body_concentrations',
    ]
    _call_scalar = ['elemental_mass_fraction', 'elemental_mole_fraction']

    _passthrough = [
        # from ThermoPhase
        'name', 'source', 'basis', 'n_elements', 'element_index',
        'element_name', 'element_names', 'atomic_weight', 'atomic_weights',
        'n_species', 'species_name', 'species_names', 'species_index',
        'species', 'n_atoms', 'molecular_weights', 'min_temp', 'max_temp',
        'reference_pressure', 'charges',
        # From Kinetics
        'n_total_species', 'n_reactions', 'n_phases',
        'kinetics_species_index', 'reaction', 'reactions', 'modify_reaction',
        'multiplier', 'set_multiplier', 'reaction_equations', 'reactant_stoich_coeff',
        'product_stoich_coeff', 'reactant_stoich_coeffs', 'product_stoich_coeffs',
        'product_stoich_coeffs_reversible',
        # from Transport
        'transport_model',
    ]

    _interface_passthrough = ['site_density']
    _interface_n_species = ['coverages']

    _purefluid_scalar = ['Q']

    def __init__(self, phase, shape=(0,), states=None, extra=None, meta={}, init=True):
        self._phase = phase
        if not init:
            return

        if states is not None:
            np_states = np.array(states)
            self.shape = np_states.shape[:-1]
            for ix, nd_ix in enumerate(self._indices):
                self._set_state(ix, np_states[nd_ix])
        elif isinstance(shape, int):
            self.shape = (shape,)
        else:
            self.shape = tuple(shape)

        def check_extra(name):
            if not isinstance(name, str):
                raise TypeError(
                    f"Unable to create extra component, passed value '{name!r}' "
                    "is not a string")
            if name in self.__dir__():
                raise ValueError(
                    f"Unable to create extra component '{name}': name is already "
                    "used by SolutionArray objects.")

        if isinstance(extra, dict):
            for name, v in extra.items():
                check_extra(name)
                ndim = self.ndim
                if not np.shape(v):
                    # initialize with scalar
                    self._add_extra(name)
                    self._set_component(name, v)
                elif (self.shape[0] == 1 or np.array(v).shape[:ndim] == self.shape):
                    arr = np.array(v)
                    if arr.dtype == object:
                        raise ValueError(
                            f"Unable to create extra component '{name}': data type "
                            "'object' is not supported.")
                    self._add_extra(name)
                    if len(arr):
                        if arr.ndim >= ndim and arr.shape[:ndim] == self.shape:
                            # direct assignment of multi-column component
                            shape = (self.size,) + arr.shape[ndim:]
                            self._set_component(name, arr.reshape(shape))
                        else:
                            self._set_component(name, arr)
                else:
                    raise ValueError(
                        f"Unable to map extra SolutionArray component named {name!r}")
        elif extra is not None:
            if self.shape != (0,):
                raise ValueError("Initial values for extra components must be "
                                 "supplied in a dictionary if the SolutionArray "
                                 "is not initially empty.")
            if isinstance(extra, np.ndarray):
                extra = extra.flatten()
            elif isinstance(extra, str):
                extra = [extra]

            try:
                iter_extra = iter(extra)
            except TypeError:
                raise ValueError(
                    "Extra components can be created by passing an iterable "
                    "of names for the components. If you want to supply initial "
                    "values for the components, use a dictionary whose keys are "
                    "the names of the components and values are the initial "
                    "values.") from None

            for name in iter_extra:
                check_extra(name)
                self._add_extra(name)

    def __getitem__(self, index):
        selected = np.arange(self.size).reshape(self.shape)[index]
        out = SolutionArray(self._phase, init=False)
        if hasattr(selected, "__len__"):
            self._share(out, selected)
        else:
            self._share(out, [selected])
        out.shape = selected.shape
        return out

    def __getattr__(self, name):
        if self._has_component(name):
            out = self._get_component(name)
            out.setflags(write=False)
            return out.reshape(self.shape + out.shape[1:])
        elif name in self.__dict__:
            super().__getattr__(name)
        else:
            raise AttributeError(
                f"{self.__class__.__name__!r} object has no attribute '{name}'")

    def __setattr__(self, name, value):
        if self._has_extra(name):
            if not self.shape:
                # scalar
                self._set_component(name, [value])
                return
            new = np.array(value)
            if not new.shape:
                # maintain shape of extra entry
                new = np.full(self.__getattr__(name).shape, value)
            elif new.shape[:len(self.shape)] != self.shape:
                raise ValueError(
                    f"Incompatible shapes for extra column '{name}': cannot assign "
                    f"value with shape {new.shape} to SolutionArray with shape "
                    f"{self.shape}")
            self._set_component(name, new)
        else:
            super().__setattr__(name, value)

    def __call__(self, *species):
        out = SolutionArray(self._phase[species], init=False)
        self._share(out, range(self.size))
        out.shape = self.shape
        return out

    def __len__(self) -> int:
        return self.shape[0]

    @property
    def ndim(self) -> int:
        """The number of dimensions in the SolutionArray.

        .. versionadded:: 3.0
        """
        return len(self.shape)

    @property
    def shape(self) -> tuple[int, ...]:
        """The shape of the SolutionArray.

        .. versionadded:: 3.0

        :return: A tuple of integers with the number of elements in each dimension.
        """
        return self._api_shape()

    @shape.setter
    def shape(self, shp):
        self._set_api_shape(shp)
        if len(shp) == 1:
            self._indices = list(range(shp[0]))
            self._output_dummy = self._indices
        else:
            self._indices = list(np.ndindex(shp))
            self._output_dummy = np.empty(shp)

    def append(self, state=None, normalize=True, **kwargs):
        """
        Append an element to the array with the specified state. Elements can
        only be appended in cases where the array of states is one-dimensional.

        The state may be specified in one of three ways:

        - as the array of [temperature, density, mass fractions] which is
          returned by `Solution.state`::

              mystates.append(gas.state)

        - as a tuple of three elements that corresponds to any of the full-state
          setters of `Solution`, such as `TPY` or `HPX`::

              mystates.append(TPX=(300, 101325, 'O2:1.0, N2:3.76'))

        - as separate keywords for each of the elements corresponding to one of
          the full-state setters::

              mystates.append(T=300, P=101325, X={'O2':1.0, 'N2':3.76})

        By default, the mass or mole fractions will be normalized i.e negative values
        are truncated and the mass or mole fractions sum up to 1.0. If this
        is not desired, the ``normalize`` argument can be set to ``False``.
        """
        if len(self.shape) != 1:
            raise IndexError("Can only append to 1D SolutionArray")

        # This check must go before we start appending to any arrays so that
        # array lengths don't get out of sync.
        missing_extra_kwargs = set(self.extra) - set(kwargs.keys())
        if missing_extra_kwargs:
            raise TypeError(
                "Missing keyword arguments for extra values: "
                "'{}'".format(", ".join(missing_extra_kwargs))
            )

        # For the checks of the state below, the kwargs dictionary can
        # only contain keywords that match properties of the state. Here
        # we pop any kwargs that have to do with the extra items so they
        # aren't included in that check. They are put into a temporary
        # storage so that appending can be done at the end of the function
        # all at once.
        extra_temp = {}
        for name in self.extra:
            extra_temp[name] = kwargs.pop(name)

        if state is not None:
            self._phase.state = state

        elif len(kwargs) == 1:
            attr, value = kwargs.popitem()
            if frozenset(attr) not in self._phase._full_states:
                raise KeyError(f"'{attr}' does not specify a full thermodynamic state")
            if normalize or attr.endswith("Q"):
                setattr(self._phase, attr, value)
            else:
                if attr.endswith("X"):
                    self._phase.set_unnormalized_mole_fractions(value[-1])
                elif attr.endswith("Y"):
                    self._phase.set_unnormalized_mass_fractions(value[-1])
                attr = attr[:-1]
                value = value[:-1]
                setattr(self._phase, attr, value)

        else:
            try:
                attr = self._phase._full_states[frozenset(kwargs)]
            except KeyError:
                raise KeyError(
                    "{} is not a valid combination of properties for setting "
                    "the thermodynamic state".format(tuple(kwargs))
                ) from None
            if normalize or attr.endswith("Q"):
                setattr(self._phase, attr, [kwargs[a] for a in attr])
            else:
                if attr.endswith("X"):
                    self._phase.set_unnormalized_mole_fractions(kwargs.pop("X"))
                elif attr.endswith("Y"):
                    self._phase.set_unnormalized_mass_fractions(kwargs.pop("Y"))
                attr = attr[:-1]
                setattr(self._phase, attr, [kwargs[a] for a in attr])

        self._append(self._phase.state, extra_temp)
        self._indices.append(len(self._indices))

    def sort(self, col, reverse=False):
        """
        Sort SolutionArray by column ``col``.

        :param col: Column that is used to sort the SolutionArray.
        :param reverse: If True, the sorted list is reversed (descending order).
        """
        if len(self.shape) != 1:
            raise TypeError("sort only works for 1D SolutionArray objects")

        indices = np.argsort(getattr(self, col))
        if reverse:
            indices = indices[::-1]
        states = [self._get_state(ix) for ix in indices]
        for loc in range(self.size):
            self._set_state(loc, states[loc])

        for k in self.extra:
            v = self._get_component(k)
            self._set_component(k, v[indices])

    def equilibrate(self, *args, **kwargs):
        """ See `ThermoPhase.equilibrate` """
        for loc in range(self.size):
            self._set_loc(loc)
            self._phase.equilibrate(*args, **kwargs)
            self._update_state(loc)

    def restore_data(self, data, normalize=True):
        """
        Restores a `SolutionArray` based on ``data`` specified in an ordered
        dictionary. Thus, this method allows to restore data exported by
        `collect_data`.

        :param data: Dictionary holding data to be restored, where keys
            refer to thermodynamic states (for example, ``T``, ``density``) or extra
            entries, and values contain corresponding data.
        :param normalize: If True, mole or mass fractions are normalized
            so that they sum up to 1.0. If False, mole or mass fractions
            are not normalized.

        The receiving `SolutionArray` either has to be empty or should have
        matching dimensions. Essential state properties and extra entries
        are detected automatically whereas stored information of calculated
        properties is omitted. If the receiving `SolutionArray` has extra
        entries already specified, only those will be imported; if *labels*
        does not contain those entries, an error is raised.
        """

        # check arguments
        if not isinstance(data, dict) or len(data) == 0:
            raise ValueError("'SolutionArray.restore_data' requires a "
                             "non-empty data dictionary")
        labels = list(data.keys())

        shape = data[labels[0]].shape
        if not shape:
            # ensure that data with a single entry have appropriate dimensions
            data = {k: np.array([v]) for k, v in data.items()}
        rows = data[labels[0]].shape[0]

        for col in data.values():
            if not isinstance(col, np.ndarray):
                raise TypeError("'SolutionArray.restore_data' only works for "
                                "dictionaries that contain ndarrays")

            if col.shape[0] != rows:
                raise ValueError("'SolutionArray.restore_data' requires "
                                 "all data entries to have a consistent "
                                 "first dimension")

        if self.shape != (0,) and self.shape != (rows,):
            raise ValueError(
                "incompatible dimensions ({} vs. {}): the receiving "
                "SolutionArray object either needs to be empty "
                "or have a length that matches data rows "
                "to be restored".format(self.shape[0], rows)
            )

        # get full state information (may differ depending on ThermoPhase type)
        states = list(self._phase._full_states.values())

        # add partial and/or potentially non-unique state definitions
        states += list(self._phase._partial_states.values())

        def join(species):
            """ Join tabular species composition data if present """
            prefix = '{}_'.format(species)
            spc = ['{}{}'.format(prefix, s) for s in self.species_names]

            # solution species names also found in labels
            valid_species = {s[2:]: labels.index(s) for s in spc
                             if s in labels}

            # labels that start with prefix (indicating concentration)
            all_species = [l[2:] for l in labels if l.startswith(prefix)]
            if valid_species:

                if len(valid_species) != len(all_species):
                    incompatible = sorted(set(valid_species) ^ set(all_species))
                    raise ValueError('incompatible species information for '
                                    '{}'.format(incompatible))

                state_comp = np.zeros((rows, self.n_species))
                for i, s in enumerate(self.species_names):
                    if s in valid_species:
                        col = data['{}{}'.format(prefix, s)]
                        if col.ndim > 1:
                            raise ValueError("Invalid format: tabular input "
                                             "for 'SolutionArray.restore_data "
                                             "requires 1D ndarrays")
                        state_comp[:, i] = col
                return state_comp
            else:
                return None

        # determine whether complete concentration is available (mass or mole)
        # assumes that 'X' or 'Y' is always in last place
        mode = ''
        state_comp = None
        for species in ['X', 'Y']:
            if not any([species in s for s in states]):
                continue
            elif species in labels:
                state_comp = data[species]
            else:
                state_comp = join(species)

            if state_comp is not None:
                # save species mode and remaining full_state candidates
                mode = species
                states = [v[:-1] for v in states if mode in v]
                break

        if mode == '':
            # concentration/quality specifier ('X' or 'Y') is not used
            states = [st.rstrip('XY') for st in states]

        # determine suitable thermo properties for reconstruction
        basis = 'mass' if self.basis == 'mass' else 'mole'
        prop = {"T": ("T", "temperature"),
                "P": ("P", "pressure"),
                "Q": ("Q", "quality"),
                "D": ("D", "density", f"density_{basis}"),
                "U": ("u", f"int_energy_{basis}"),
                "V": ("v", f"volume_{basis}"),
                "H": ("h", f"enthalpy_{basis}"),
                "S": ("s", f"entropy_{basis}")}
        for st in states:
            # identify property specifiers
            state = [{st[i]: p for p in prop[st[i]] if p in labels}
                     for i in range(len(st))]
            if all(state):
                # all property identifiers match
                mode = st + mode
                state = [list(st.values())[0] for st in state]
                break
        if len(mode) == 1:
            raise ValueError(
                "invalid/incomplete state information (detected "
                "partial information as mode='{}')".format(mode)
            )

        # assemble and restore state information
        state_data = tuple([data[s] for s in state])
        if state_comp is not None:
            state_data += (state_comp,)

        # labels may include calculated properties that must not be restored:
        # compare column labels to names that are reserved for SolutionArray
        # attributes (see `SolutionArray.collect_data`), such as scalar values,
        # arrays with number of species, and arrays with number of reactions.
        exclude = [lab for lab in labels
                   if any([lab in self._scalar,
                           '_'.join(lab.split('_')[:-1]) in self._n_species,
                           lab.split(' ')[0] in self._n_reactions])]
        exclude += ['X', 'Y']
        extra = {lab: data[lab] for lab in labels
                 if lab not in exclude}
        if len(self.extra):
            extra_lists = {k: extra[k] for k in self.extra}
        else:
            extra_lists = extra

        # ensure that SolutionArray accommodates dimensions
        if self.shape == (0,):
            self.shape = (rows,)
        else:
            self.resize(np.prod(self.shape))

        # restore data
        if normalize or mode.endswith("Q"):
            for loc, i in enumerate(self._indices):
                setattr(self._phase, mode, [st[i, ...] for st in state_data])
                self._update_state(loc)
        else:
            for loc, i in enumerate(self._indices):
                if mode.endswith("X"):
                    self._phase.set_unnormalized_mole_fractions(
                        [st[i, ...] for st in state_data][2]
                    )
                    setattr(self._phase, mode[:2],
                            [st[i, ...] for st in state_data[:2]])
                elif mode.endswith("Y"):
                    self._phase.set_unnormalized_mass_fractions(
                        [st[i, ...] for st in state_data][2]
                    )
                    setattr(self._phase, mode[:2],
                            [st[i, ...] for st in state_data[:2]])
                else:
                    setattr(self._phase, mode, [st[i, ...] for st in state_data])
                self._update_state(loc)

        for key, value in extra_lists.items():
            if not self._has_component(key):
                self._add_extra(key)
            self._set_component(key, value)

    def set_equivalence_ratio(self, phi, *args, **kwargs):
        """
        See `ThermoPhase.set_equivalence_ratio`

        Note that ``phi`` either needs to be a scalar value or dimensions have
        to be matched to the `SolutionArray`.
        """

        # If ``phi`` is lower-dimensional than the SolutionArray's shape (for
        # example, a scalar), broadcast it to have the right number of
        # dimensions.
        phi, _ = np.broadcast_arrays(phi, self._output_dummy)

        for loc, index in enumerate(self._indices):
            self._set_loc(loc)
            self._phase.set_equivalence_ratio(phi[index], *args, **kwargs)
            self._update_state(loc)

    def set_mixture_fraction(self, mixture_fraction, *args, **kwargs):
        """
        See `ThermoPhase.set_mixture_fraction`

        Note that ``mixture_fraction`` either needs to be a scalar value or
        dimensions have to be matched to the `SolutionArray`.
        """

        # If ``mixture_fraction`` is lower-dimensional than the SolutionArray's
        # shape (for example, a scalar), broadcast it to have the right number
        # of dimensions.
        mixture_fraction, _ = np.broadcast_arrays(mixture_fraction, self._output_dummy)

        for loc, index in enumerate(self._indices):
            self._set_loc(loc)
            self._phase.set_mixture_fraction(mixture_fraction[index], *args, **kwargs)
            self._update_state(loc)

    def collect_data(self, cols=None, tabular=False, threshold=0, species=None):
        """
        Returns the data specified by ``cols`` in a dictionary, where keys correspond
        to `SolutionArray` attributes to be exported.

        :param cols: A list of any properties of the solution that are scalars
            or which have a value for each species or reaction. If species names
            are specified, then either the mass or mole fraction of that species
            will be taken, depending on the value of ``species``. ``cols`` may also
            include any arrays which were specified as ``extra`` variables when
            defining the `SolutionArray` object. The special value 'extra' can
            be used to include all 'extra' variables.
        :param tabular: Split 2D data into separate 1D columns for each
            species / reaction
        :param threshold: Relative tolerance for including a particular column
            if tabular output is enabled. The tolerance is applied by comparing
            the maximum absolute value for a particular column to the maximum
            absolute value in all columns for the same variable (for example, mass
            fraction).
        :param species: Specifies whether to use mass ('Y') or mole ('X')
            fractions for individual species specified in 'cols'
        """
        if tabular and len(self.shape) != 1:
            raise AttributeError("Tabular output of collect_data only works "
                                 "for 1D SolutionArray")

        # Create default columns (including complete state information)
        if cols is None:
            cols = ('extra',) + self._phase._native_state

        if species is None:
            species = cols[-1]
        else:
            cols = cols[:-1] + (species,)

        # Expand cols to include the individual items in 'extra'
        expanded_cols = []
        for c in cols:
            if c == 'extra':
                expanded_cols.extend(self.extra)
            else:
                expanded_cols.append(c)

        expanded_cols = tuple(['density' if c == 'D' else c
                               for c in expanded_cols])

        def split(c, d):
            """ Split attribute arrays into columns for tabular output """
            # Determine labels for the items in the current group of columns
            if c in self._n_species:
                collabels = ['{}_{}'.format(c, s) for s in self.species_names]
            elif c in self._n_reactions:
                collabels = ['{} {}'.format(c, r)
                             for r in self.reaction_equations()]
            elif c in species_names:
                collabels = ['{}_{}'.format(species, c)]
            elif c in self.extra and d.ndim > 1:
                raise NotImplementedError(
                    "Detected multi-dimensional extra column '{}': "
                    "tabular output is not supported.".format(c))
            else:
                collabels = [c]

            if d.ndim > 1:
                if threshold and d.size:
                    # Determine threshold value and select columns to keep
                    maxval = abs(d).max()
                    keep = (abs(d) > threshold * maxval).any(axis=0)
                    d = d[:, keep]
                    collabels = [label for label, k in zip(collabels, keep) if k]

                return [(cl, d[:, i]) for i, cl in enumerate(collabels)]
            else:
                return [(collabels[0], d)]

        data = []
        attrs = self.__dir__() + self.component_names
        species_names = set(self.species_names)
        for c in expanded_cols:
            if c in species_names:
                d = getattr(self(c), species)
            elif c in attrs:
                d = getattr(self, c)
            else:
                raise AttributeError(f"Component '{c}' not supported")

            if tabular:
                data += split(c, d)

            else:
                data += [(c, d)]

        return dict(data)

    def read_csv(self, filename, normalize=True):
        """
        Read a CSV file named ``filename`` and restore data to the `SolutionArray`
        using `restore_data`. This method allows for recreation of data
        previously exported by `write_csv`.

        The ``normalize`` argument is passed on to `restore_data` to normalize
        mole or mass fractions. By default, ``normalize`` is ``True``.
        """
        try:
            # pandas handles escaped entries correctly
            _import_pandas()
            df = _pandas.read_csv(filename)
            self.from_pandas(df)
            return
        except ImportError:
            pass

        # fall back to numpy; this works unless CSV file contains escaped entries
        data = np.genfromtxt(filename, delimiter=',', deletechars='',
                                dtype=None, names=True, encoding=None)
        data_dict = {label: data[label] for label in data.dtype.names}
        self.restore_data(data_dict, normalize)

    def to_pandas(self, cols=None, *args, **kwargs):
        """
        Returns the data specified by ``cols`` in a single `pandas.DataFrame`.

        Additional arguments are passed on to `collect_data`. This method works
        only with 1D `SolutionArray` objects and requires a working *pandas*
        installation. Use pip or conda to install ``pandas`` to enable this method.
        """
        if not _pandas:
            _import_pandas()

        data_dict = self.collect_data(*args, cols=cols, tabular=True, **kwargs)
        data = np.hstack([d[:, np.newaxis] for d in data_dict.values()])
        labels = list(data_dict.keys())
        return _pandas.DataFrame(data=data, columns=labels)

    def from_pandas(self, df, normalize=True):
        """
        Restores `SolutionArray` data from a `pandas.DataFrame` ``df``.

        This method is intended for loading of data that were previously
        exported by `to_pandas`. The method requires a working *pandas*
        installation. The package ``pandas`` can be installed using pip or conda.

        The ``normalize`` argument is passed on to `restore_data` to normalize
        mole or mass fractions. By default, ``normalize`` is ``True``.
        """
        data_dict = {}
        for label in list(df.columns):
            data_dict[label] = df[label].values
            if data_dict[label].dtype.type == np.object_:
                # convert object columns to string
                data_dict[label] = data_dict[label].astype('U')
        self.restore_data(data_dict, normalize)

    def save(self, fname, name=None, sub=None, description=None, *,
             overwrite=False, compression=0, basis=None):
        """
        Save current `SolutionArray` contents to a data file.

        Data can be saved either in CSV format (extension ``*.csv``), YAML container
        format (extension ``*.yaml``/``*.yml``) or HDF container format (extension
        ``*.h5``/``*.hdf5``/``*.hdf``). The output format is automatically inferred from
        the file extension.

        CSV files preserve state data and auxiliary data for a single `SolutionArray` in
        a comma-separated text format, container files may hold multiple `SolutionArray`
        entries in an internal hierarchical structure. While YAML is a human-readable
        text format, HDF is a binary format that supports compression and is recommended
        for large datasets.

        For container files (YAML and HDF), header information contains automatically
        generated time stamps, version information and an optional description.
        Container files also preserve `SolutionArray` metadata (example:
        `SolutionArray` objects generated by `Sim1D` store simulation settings).

        :param fname:
            Name of output file (CSV, YAML or HDF)
        :param name:
            Identifier of storage location within the container file; this node/group
            contains header information and a subgroup holding actual `SolutionArray`
            data (YAML/HDF only).
        :param sub:
            Name identifier for the subgroup holding the `SolutionArray` data and
            metadata objects. If `None`, the subgroup name defaults to ``data``
            (YAML/HDF only).
        :param description:
            Custom comment describing the dataset to be stored (YAML/HDF only).
        :param overwrite:
            Force overwrite if file/name exists; optional (default=`False`)
        :param compression:
            Compression level (0-9); optional (default=0; HDF only)
        :param basis:
            Output mass (``Y``/``mass``) or mole (``Y``/``mass``) fractions;
            if not specified (`None`), the native basis of the underlying `ThermoPhase`
            manager is used.

        .. versionadded:: 3.0
        """
        self._cxx_save(fname, name, sub, description, overwrite, compression, basis)

    def restore(self, fname, name=None, sub=None):
        """
        Restore `SolutionArray` data and header information from a container file.

        This method retrieves data from a YAML or HDF files that were previously saved
        using the `save` method.

        :param fname:
            Name of container file (YAML or HDF)
        :param name:
            Identifier of location within the container file; this node/group contains
            header information and a subgroup holding actual `SolutionArray` data
        :param sub:
            Name identifier for the subgroup holding the `SolutionArray` data and
            metadata objects. If `None`, the subgroup name defaults to ``data``
        :return:
            Dictionary holding `SolutionArray` meta data.

        .. versionadded:: 3.0
        """
        meta = self._cxx_restore(fname, name, sub)

        # ensure self._indices and self._output_dummy are set
        self.shape = self._api_shape()
        return meta

    def __reduce__(self):
        raise NotImplementedError('SolutionArray object is not picklable')

    def __copy__(self):
        raise NotImplementedError('SolutionArray object is not copyable')


def _state2_prop(name, doc_source):
    # Factory for creating properties which consist of a tuple of two variables,
    # such as 'TP' or 'SV'
    def getter(self):
        a = np.empty(self.shape)
        b = np.empty(self.shape)
        for loc, index in enumerate(self._indices):
            self._set_loc(loc)
            a[index], b[index] = getattr(self._phase, name)
        return a, b

    def setter(self, AB):
        if len(AB) != 2:
            raise ValueError("Expected 2 elements, got {}".format(len(AB)))
        A, B, _ = np.broadcast_arrays(AB[0], AB[1], self._output_dummy)
        for loc, index in enumerate(self._indices):
            self._set_loc(loc)
            setattr(self._phase, name, (A[index], B[index]))
            self._update_state(loc)

    return getter, setter


def _state3_prop(name, doc_source, scalar=False):
    # Factory for creating properties which consist of a tuple of three
    # variables, such as 'TPY' or 'UVX'
    def getter(self):
        a = np.empty(self.shape)
        b = np.empty(self.shape)
        if scalar:
            c = np.empty(self.shape)
        else:
            c = np.empty(self.shape + (self._phase.n_selected_species,))
        for loc, index in enumerate(self._indices):
            self._set_loc(loc)
            a[index], b[index], c[index] = getattr(self._phase, name)
        return a, b, c

    def setter(self, ABC):
        if len(ABC) != 3:
            raise ValueError("Expected 3 elements, got {}".format(len(ABC)))
        A, B, _ = np.broadcast_arrays(ABC[0], ABC[1], self._output_dummy)
        XY = ABC[2] # composition
        if len(np.shape(XY)) < 2:
            # composition is a single array (or string or dict)
            for loc, index in enumerate(self._indices):
                self._set_loc(loc)
                setattr(self._phase, name, (A[index], B[index], XY))
                self._update_state(loc)
        else:
            # composition is an array with trailing dimension n_species
            C = np.empty(self.shape + (self._phase.n_selected_species,))
            C[:] = XY
            for loc, index in enumerate(self._indices):
                self._set_loc(loc)
                setattr(self._phase, name, (A[index], B[index], C[index]))
                self._update_state(loc)

    return getter, setter

def _make_functions():
    # this is wrapped in a function to avoid polluting the module namespace

    names = []
    for ph, ext in [(ThermoPhase, 'XY'), (PureFluid, 'Q')]:

        # all state setters/getters are combination of letters
        setters = 'TDPUVHS' + ext
        scalar = ext == 'Q'

        # obtain setters/getters from thermo objects
        all_states = [k for k in ph.__dict__
                      if not set(k) - set(setters) and len(k)>1]

        # state setters (copy from ThermoPhase objects)
        for name in all_states:
            if name in names:
                continue
            names.append(name)
            if len(name) == 2:
                getter, setter = _state2_prop(name, ph)
            elif len(name) == 3:
                getter, setter = _state3_prop(name, ph, scalar)
            doc = getattr(ph, name).__doc__
            prop = property(getter, setter, doc=doc)
            setattr(SolutionArray, name, prop)

    # Functions which define empty output arrays of an appropriate size for
    # different properties
    def empty_scalar(self):
        return np.empty(self.shape)

    def empty_strings(self):
        # The maximum length of strings assigned by built-in methods is
        # currently limited to 50 characters; an attempt to assign longer
        # character arrays will result in truncated strings.
        return np.empty(self.shape, dtype='U50')

    def empty_species(self):
        return np.empty(self.shape + (self._phase.n_selected_species,))

    def empty_total_species(self):
        n_tot = self._phase.n_total_species
        # account for deselected species
        n_tot -= self._phase.n_species - self._phase.n_selected_species
        return np.empty(self.shape + (n_tot,))

    def empty_species2(self):
        return np.empty(self.shape + (self._phase.n_species, self._phase.n_species))

    def empty_reactions(self):
        return np.empty(self.shape + (self._phase.n_reactions,))

    # Factory for creating read-only properties
    def make_prop(name, get_container, doc_source, block_interface=False):
        def getter(self):
            if block_interface and isinstance(self._phase, Interface):
                # used to block Interface methods that require synchronized updates of
                # linked phases
                raise NotImplementedError(
                    "Method not implemented for SolutionArray containing Interface.")
            v = get_container(self)
            for loc, index in enumerate(self._indices):
                self._set_loc(loc)
                v[index] = getattr(self._phase, name)
            return v
        return property(getter, doc=getattr(doc_source, name).__doc__)

    for name in SolutionArray._scalar:
        setattr(SolutionArray, name, make_prop(name, empty_scalar, Solution))

    for name in SolutionArray._strings:
        setattr(SolutionArray, name, make_prop(name, empty_strings, Solution))

    for name in SolutionArray._n_species:
        setattr(SolutionArray, name, make_prop(name, empty_species, Solution))

    for name in SolutionArray._interface_n_species:
        setattr(SolutionArray, name, make_prop(name, empty_species, Interface))

    for name in SolutionArray._purefluid_scalar:
        setattr(SolutionArray, name, make_prop(name, empty_scalar, PureFluid))

    for name in SolutionArray._n_total_species:
        setattr(SolutionArray, name,
                make_prop(name, empty_total_species, Solution, True))

    for name in SolutionArray._n_species2:
        setattr(SolutionArray, name, make_prop(name, empty_species2, Solution, True))

    for name in SolutionArray._n_reactions:
        setattr(SolutionArray, name, make_prop(name, empty_reactions, Solution, True))

    # Factory for creating wrappers for functions which return a value
    def caller(name, get_container):
        def wrapper(self, *args, **kwargs):
            v = get_container(self)
            for loc, index in enumerate(self._indices):
                self._set_loc(loc)
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

    def passthrough_method(orig):
        def f(self, *args, **kwargs):
            return orig(self._phase, *args, **kwargs)
        f.__doc__ = orig.__doc__
        return f

    for name in SolutionArray._passthrough:
        orig = getattr(Solution, name)
        if hasattr(orig, "__call__"):
            setattr(SolutionArray, name, passthrough_method(orig))
        else:
            setattr(SolutionArray, name, passthrough_prop(name, Solution))

    for name in SolutionArray._interface_passthrough:
        orig = getattr(Interface, name)
        if hasattr(orig, "__call__"):
            setattr(SolutionArray, name, passthrough_method(orig))
        else:
            setattr(SolutionArray, name, passthrough_prop(name, Interface))

_make_functions()
