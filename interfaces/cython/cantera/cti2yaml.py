#!/usr/bin/env python

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""Convert legacy CTI input files to YAML.

There are two main entry points to this script, `main` and `convert`. The former is
used from the command line interface and parses the arguments passed. The latter
accepts either the name of the CTI input file or a string containing the CTI
content.
"""

import sys
import re
import pathlib
from collections import OrderedDict
import numpy as np
from email.utils import formatdate
import argparse

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml


def _printerr(*args):
    # All debug and error output should go to stderr
    print(*args, file=sys.stderr)


class InputError(Exception):
    """
    Exception raised if an error is encountered while parsing the input file.
    """
    def __init__(self, msg, *args):
        if args:
            msg = msg.format(*args)
        super().__init__(msg)

BlockMap = yaml.comments.CommentedMap

def FlowMap(*args, **kwargs):
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m

def FlowList(*args, **kwargs):
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst

# Improved float formatting requires Numpy >= 1.14
if hasattr(np, 'format_float_positional'):
    def float2string(data):
        if data == 0:
            return '0.0'
        elif 0.01 <= abs(data) < 10000:
            return np.format_float_positional(data, trim='0')
        else:
            return np.format_float_scientific(data, trim='0')
else:
    def float2string(data):
        return repr(data)

def represent_float(self, data):
    # type: (Any) -> Any
    if data != data:
        value = '.nan'
    elif data == self.inf_value:
        value = '.inf'
    elif data == -self.inf_value:
        value = '-.inf'
    else:
        value = float2string(data)

    return self.represent_scalar(u'tag:yaml.org,2002:float', value)

yaml.RoundTripRepresenter.add_representer(float, represent_float)


def applyUnits(value):
    if isinstance(value, (float, int)):
        return value
    else:
        units = value[1]
        units = re.sub(r'([A-Za-z])-([A-Za-z])', r'\1*\2', units)
        units = re.sub(r'([A-Za-z])([-\d])', r'\1^\2', units)
        return '{} {}'.format(float2string(value[0]), units)


# map of old CTI/XML names to the new YAML names
_newNames = {
    'GasKinetics': 'gas',
    'Interface': 'surface',
    'Edge': 'edge',
    'Mix': 'mixture-averaged',
    'Multi': 'multicomponent',
    'Ion': 'ionized-gas',
    'molar_volume': 'species-molar-volume',
    'solvent_volume': 'solvent-molar-volume',
    'unity': 'unity'
}

# constants that can be used in .cti files
OneAtm = 1.01325e5
OneBar = 1.0e5
# Conversion from eV to J/kmol (electron charge * Avogadro constant)
eV = 9.64853364595687e7
# Electron Mass in kg
ElectronMass = 9.10938291e-31

# default units
_ulen = 'm'
_umol = 'kmol'
_umass = 'kg'
_utime = 's'
_ue = 'J/kmol'
_uenergy = 'J'
_upres = 'Pa'

# default std state pressure
_pref = OneAtm

_name = 'noname'

# these lists store top-level entries
_elements = []
_species = []
_speciesnames = []
_phases = []
_reactions = {'reactions': []}

# default for Motz & Wise correction
_motz_wise = None

def enable_motz_wise():
    """
    Enable the Motz & Wise correction by default for all sticking reactions.
    """
    global _motz_wise
    _motz_wise = True

def disable_motz_wise():
    """
    Disable the Motz & Wise correction by default for all sticking reactions.
    """
    global _motz_wise
    _motz_wise = False

def validate(species = 'yes', reactions = 'yes'):
    pass

def dataset(nm):
    "Set the dataset name. Invoke this to change the name of the YAML file."
    global _name
    _name = nm

def standard_pressure(p0):
    """Set the default standard-state pressure."""
    global _pref
    _pref = p0

def units(length = '', quantity = '', mass = '', time = '',
          act_energy = '', energy = '', pressure = ''):
    """
    Set the default units.

    :param length:
        The default units for length. Default: ``'m'``
    :param mass:
        The default units for mass. Default: ``'kg'``
    :param quantity:
        The default units to specify number of molecules. Default: ``'kmol'``
    :param time:
        The default units for time. Default: ``'s'``
    :param energy:
        The default units for energies. Default: ``'J'``
    :param act_energy:
        The default units for activation energies. Default: ``'K'``
    :param pressure:
        The default units for pressure. Default: ``'Pa'``
    """
    global _ulen, _umol, _ue, _utime, _umass, _uenergy, _upres
    if length: _ulen = length
    if quantity: _umol = quantity
    if act_energy: _ue = act_energy
    if time: _utime = time
    if mass: _umass = mass
    if energy: _uenergy = energy
    if pressure: _upres = pressure


def get_composition(atoms):
    if isinstance(atoms, dict): return atoms
    a = atoms.replace(',',' ')
    toks = a.split()
    d = OrderedDict()
    for t in toks:
        b = t.split(':')
        try:
            d[b[0]] = int(b[1])
        except ValueError:
            d[b[0]] = float(b[1])
    return d


class element:
    """ An atomic element or isotope. """
    def __init__(self, symbol='', atomic_mass=0.01, atomic_number=None):
        """
        :param symbol:
            The symbol for the element or isotope.
        :param atomic_mass:
            The atomic mass in amu.
        :param atomic_number:
            The atomic number of the element or isotope. Optional.
        """
        self.symbol = symbol
        self.atomic_weight = atomic_mass
        self.atomic_number = atomic_number
        _elements.append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('symbol', node.symbol),
                        ('atomic-weight', node.atomic_weight)])
        if node.atomic_number is not None:
            out['atomic-number'] = node.atomic_number
        return representer.represent_dict(out)


class species:
    """A constituent of a phase or interface."""

    def __init__(self, name, atoms='', note = '', thermo=None, transport=None,
                 charge=None, size=1.0, standardState=None):
        """
        :param name:
            The species name (or formula). The name may be arbitrarily long,
            although usually a relatively short, abbreviated name is most
            convenient. Required parameter.
        :param atoms:
            The atomic composition, specified by a string containing
            space-delimited ``<element>:<atoms>`` pairs. The number of atoms may be
            either an integer or a floating-point number.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param thermo:
            The parameterization to use to compute the reference-state
            thermodynamic properties. This must be one of the entry types
            described in `Thermodynamic Property Models
            <https://cantera.org/science/science-species.html#sec-thermo-models>`__.
            To specify multiple parameterizations, each for a different temperature range,
            group them in parentheses.
        :param transport:
            An entry specifying parameters to compute this species'
            contribution to the transport properties. This must be one of the
            entry types described in `Species Transport Coefficients
            <https://cantera.org/science/science-species.html#species-transport-coefficients>`__,
            and must be consistent with the transport model of the phase into which
            the species is imported. To specify parameters for multiple
            transport models, group the entries in parentheses.
        :param charge:
            The charge, in multiples of :math:`|e|`. If not specified, the
            charge will be calculated from the number of "atoms" of element
            ``E``, which represents an electron.
        :param size:
            The species "size". Currently used only for surface species,
            where it represents the number of sites occupied.
        :param standardState:
            The species standard state model. Currently used only for
            `IdealSolidSolution` and derived classes where it is used to
            calculate the phase density.
        """
        self.name = name
        self.atoms = get_composition(atoms)
        if charge is not None and 'E' not in self.atoms:
            self.atoms['E'] = -charge
        self.size = size
        self.comment = note

        if isinstance(thermo, (list, tuple)):
            if isinstance(thermo[0], (NASA, NASA9, Shomate)):
                self.thermo = MultiPolyThermo(thermo)
        elif isinstance(thermo, (NASA, NASA9, Shomate)):
            self.thermo = MultiPolyThermo([thermo])
        elif thermo is not None:
            self.thermo = thermo
        else:
            self.thermo = const_cp()

        self.transport = transport
        self.standard_state = standardState

        self.rk_pure = {}
        self.rk_binary = {}
        self.density = None

        _species.append(self)
        _speciesnames.append(name)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('name', node.name),
                        ('composition', FlowMap(node.atoms.items()))])
        if node.size != 1:
            out['sites'] = node.size
        out['thermo'] = node.thermo
        if node.density:
            out['equation-of-state'] = {
                'model': 'constant-volume',
                'density': applyUnits(node.density)
            }

        if node.rk_pure:
            a = node.rk_pure['a']
            if isinstance(a, (tuple, list)):
                a = FlowList([applyUnits(ai) for ai in a])
            else:
                a = applyUnits(a)
            out['equation-of-state'] = {
                'model': 'Redlich-Kwong',
                'a': a,
                'b': applyUnits(node.rk_pure['b'])
            }

        if node.rk_binary:
            rkbin = BlockMap()
            for species, a in node.rk_binary.items():
                if isinstance(a, (tuple, list)):
                    rkbin[species] = FlowList([applyUnits(ai) for ai in a])
                else:
                    rkbin[species] = applyUnits(a)
            out['equation-of-state']['binary-a'] = rkbin

        if node.standard_state:
            out['equation-of-state'] = {
                'model': 'constant-volume',
                'molar-volume': applyUnits(node.standard_state.molar_volume)
            }

        if node.transport:
            out['transport'] = node.transport
        if node.comment:
            comment = node.comment.split("\n")
            comment = " ".join(c.strip() for c in comment)
            out['note'] = comment
        return representer.represent_dict(out)


class thermo:
    """Base class for species thermodynamic properties."""
    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['model'] = self.model
        pref = self.pref or _pref
        if pref != OneAtm:
            out['reference-pressure'] = pref


class NASA(thermo):
    """The 7-coefficient NASA polynomial parameterization."""
    def __init__(self, Trange=(0.0, 0.0), coeffs=(), p0=None):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required.
        :param coeffs:
            List of seven coefficients :math:`(a_0, \ldots , a_6)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value is used, which is set by the ``standard_pressure``
            directive.
        """
        self.model = 'NASA7'
        self.T_range = Trange
        self.pref = p0
        if len(coeffs) != 7:
            raise InputError('NASA coefficient list must have length = 7')
        self.coeffs = coeffs


class NASA9(thermo):
    """NASA9 polynomial parameterization for a single temperature region."""

    def __init__(self, Trange=(0.0, 0.0), coeffs=(), p0=None):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required.
        :param coeffs:
            List of nine coefficients :math:`(a_0, \ldots , a_8)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value is used, which is set by the ``standard_pressure``
            directive.
        """
        self.model = 'NASA9'
        self.T_range = Trange
        self.pref = p0
        if len(coeffs) != 9:
            raise InputError('NASA9 coefficient list must have length = 9')
        self.coeffs = coeffs


class MultiPolyThermo(thermo):
    def __init__(self, regions):
        regions = sorted(regions, key=lambda r: r.T_range[0])
        self.pref = regions[0].pref
        self.Tranges = [regions[0].T_range[0]]
        self.model = regions[0].model
        self.data = []
        for r in regions:
            self.Tranges.append(r.T_range[1])
            self.data.append(r.coeffs)

    def get_yaml(self, out):
        super().get_yaml(out)
        out['temperature-ranges'] = FlowList(self.Tranges)
        out['data'] = [FlowList(coeffs) for coeffs in self.data]


class Shomate(thermo):
    """Shomate polynomial parameterization."""
    def __init__(self, Trange=(0.0, 0.0), coeffs=(), p0=None):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required input.
        :param coeffs:
            Sequence of seven coefficients :math:`(A, \ldots , G)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value set by the ``standard_pressure`` directive is used.
        """
        self.model = 'Shomate'
        self.T_range = Trange
        self.pref = p0
        if len(coeffs) != 7:
            raise InputError('Shomate coefficient list must have length = 7')
        self.coeffs = coeffs


class const_cp(thermo):
    """Constant specific heat."""

    def __init__(self, t0=None, cp0=None, h0=None, s0=None, tmax=None,
                 tmin=None):
        """
        :param t0:
            Temperature parameter T0. Default: 298.15 K.
        :param cp0:
            Reference-state molar heat capacity (constant). Default: 0.0.
        :param h0:
            Reference-state molar enthalpy at temperature T0. Default: 0.0.
        :param s0:
            Reference-state molar entropy at temperature T0. Default: 0.0.
        """
        self.model = 'constant-cp'
        self.pref = None
        self.t0 = t0
        self.h0 = h0
        self.s0 = s0
        self.cp0 = cp0

    def get_yaml(self, out):
        super().get_yaml(out)
        if self.t0 is not None:
            out['T0'] = applyUnits(self.t0)
        if self.h0 is not None:
            out['h0'] = applyUnits(self.h0)
        if self.s0 is not None:
            out['s0'] = applyUnits(self.s0)
        if self.cp0 is not None:
            out['cp0'] = applyUnits(self.cp0)


class gas_transport:
    """
    Species-specific Transport coefficients for gas-phase transport models.
    """
    def __init__(self, geom, diam, well_depth, dipole=0.0, polar=0.0,
                 rot_relax=0.0, acentric_factor=None, disp_coeff=0.0,
                 quad_polar=0.0):
        """
        :param geom:
            A string specifying the molecular geometry. One of ``atom``,
            ``linear``, or ``nonlinear``. Required.
        :param diam:
            The Lennard-Jones collision diameter in Angstroms. Required.
        :param well_depth:
            The Lennard-Jones well depth in Kelvin. Required.
        :param dipole:
            The permanent dipole moment in Debye. Default: 0.0
        :param polar:
            The polarizability in A^3. Default: 0.0
        :param rot_relax:
            The rotational relaxation collision number at 298 K. Dimensionless.
            Default: 0.0
        :param acentric_factor:
            Pitzer's acentric factor.  Dimensionless.
            Default: 0.0
        :param disp_coeff:
            The dispersion coefficient in A^5
            Default: 0.0
        :param quad_polar:
            The quadrupole polarizability
            Default: 0.0
        """
        self.geometry = geom
        self.diameter = diam
        self.well_depth = well_depth
        self.dipole = dipole
        self.polarizability = polar
        self.rot_relax = rot_relax
        self.acentric_factor = acentric_factor
        self.disp_coeff = disp_coeff
        self.quad_polar = quad_polar

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('model', 'gas'),
                        ('geometry', node.geometry),
                        ('diameter', node.diameter),
                        ('well-depth', node.well_depth)])
        if node.dipole:
            out['dipole'] = node.dipole
        if node.polarizability:
            out['polarizability'] = node.polarizability
        if node.rot_relax:
            out['rotational-relaxation'] = node.rot_relax
        if node.acentric_factor:
            out['acentric-factor'] = node.acentric_factor
        if node.disp_coeff:
            out['dispersion-coefficient'] = node.disp_coeff
        if node.quad_polar:
            out['quadrupole-polarizability'] = node.quad_polar
        return representer.represent_dict(out)


class Arrhenius:
    def __init__(self, A=0.0, b=0.0, E=0.0, coverage=()):
        """
        :param A:
            The pre-exponential coefficient. Required input. If entered without
            units, the units will be computed considering all factors that
            affect the units. The resulting units string is written to the CTML
            file individually for each reaction pre-exponential coefficient.
        :param b:
            The temperature exponent. Dimensionless. Default: 0.0.
        :param E:
            Activation energy. Default: 0.0.
        :param coverage:
            For a single coverage dependency, a list with four elements: the
            species name followed by the three coverage parameters. For multiple
            coverage dependencies, a list of lists containing the individual
            sets of coverage parameters. Only used for surface and edge
            reactions.
        """

        self.A = A
        self.b = b
        self.E = E

        if coverage:
            if isinstance(coverage[0], str):
                self.coverage = [coverage]
            else:
                self.coverage = coverage
            for cov in self.coverage:
                if len(cov) != 4:
                    raise InputError("Incorrect number of coverage parameters")
        else:
            self.coverage = None

    @classmethod
    def to_yaml(cls, representer, node):
        out = FlowMap([('A', applyUnits(node.A)),
                       ('b', applyUnits(node.b)),
                       ('Ea', applyUnits(node.E))])
        return representer.represent_dict(out)


class stick(Arrhenius):
    """
    A rate expression for a surface reaction given as a sticking probability,
    parameterized using a modified Arrhenius expression.
    """
    def __init__(self, *args, **kwargs):
        """
        :param motz_wise:
            ``True`` if the Motz & Wise correction should be used, ``False`` if
            not. If unspecified, use the mechanism default (set using the
            functions `enable_motz_wise` or `disable_motz_wise`).
        """
        self.motz_wise = kwargs.pop('motz_wise', None)
        Arrhenius.__init__(self, *args, **kwargs)


class reaction:
    """
    A homogeneous chemical reaction with pressure-independent rate coefficient
    and mass-action kinetics.
    """
    def __init__(self, equation, kf, id='', order='', options=()):
        r"""
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function :math:`A T^b exp(-E/\hat{R}T)`.
        :param id:
            An optional identification string.
        :param order:
            Override the default reaction orders implied by the reactant
            stoichiometric coefficients. Given as a string of key:value pairs,
            e.g., ``"CH4:0.25 O2:1.5"``.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
            May be one or more (as a list) of the following: ``'duplicate'``,
            ``'negative_A'``,`` 'negative_orders'``, ``'nonreactant_orders'``.
        """
        self.equation = equation
        self.order = get_composition(order)
        self.number = len(_reactions['reactions']) + 1
        self.id = id
        self.options = [options] if isinstance(options, str) else options
        self.kf = Arrhenius(*kf) if isinstance(kf, (list, tuple)) else kf
        self.type = 'elementary'
        _reactions['reactions'].append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['equation'] = self.equation
        out.yaml_add_eol_comment('Reaction {}'.format(self.number), 'equation')
        if self.type not in ('elementary', 'edge', 'surface'):
            out['type'] = self.type

        if self.id:
            out['id'] = self.id

        if self.type in ('elementary', 'three-body', 'edge', 'surface'):
            out['rate-constant'] = self.kf

        if 'duplicate' in self.options:
            out['duplicate'] = True
        if 'negative_A' in self.options:
            out['negative-A'] = True

        if self.order:
            out['orders'] = FlowMap(self.order.items())
        if 'negative_orders' in self.options:
            out['negative-orders'] = True
        if 'nonreactant_orders' in self.options:
            out['nonreactant-orders'] = True


class three_body_reaction(reaction):
    """
    A three-body reaction.
    """
    def __init__(self, equation, kf, efficiencies='', id='', options=()):
        """
        :param equation:
            A string specifying the chemical equation. The reaction can be
            written in either the association or dissociation directions, and
            may be reversible or irreversible.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies.
            The efficiencies for unspecified species are set to 1.0.
        :param id:
            An optional identification string.
        :param options: Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        """
        super().__init__(equation, kf, id, '', options)
        self.type = 'three-body'
        self.efficiencies = get_composition(efficiencies)

    def get_yaml(self, out):
        super().get_yaml(out)
        if self.efficiencies:
            out['efficiencies'] = FlowMap(self.efficiencies)


class falloff_base(reaction):
    """ Base class for falloff_reaction and chemically_activated_reaction """
    def __init__(self, equation, klow, khigh, efficiencies, falloff, id, options):
        super().__init__(equation, None, id, '', options)
        self.k_low = Arrhenius(*klow) if isinstance(klow, (list, tuple)) else klow
        self.k_high = Arrhenius(*khigh) if isinstance(khigh, (list, tuple)) else khigh
        self.falloff = falloff
        self.efficiencies = get_composition(efficiencies)

    def get_yaml(self, out):
        super().get_yaml(out)

        out['low-P-rate-constant'] = self.k_low
        out['high-P-rate-constant'] = self.k_high

        if self.falloff:
            self.falloff.get_yaml(out)

        if self.efficiencies:
            out['efficiencies'] = FlowMap(self.efficiencies)


class falloff_reaction(falloff_base):
    """ A gas-phase falloff reaction. """
    def __init__(self, equation, kf0, kf, efficiencies='', falloff=None, id='',
                 options=()):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param kf0:
            The rate coefficient for the forward direction in the low-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies. The
            efficiency for unspecified species is set to 1.0.
        :param falloff:
            An embedded entry specifying a falloff function. If omitted, a
            unity falloff function (Lindemann form) will be used.
        :param id:
            An optional identification string.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        """
        super().__init__(equation, kf0, kf, efficiencies, falloff, id, options)
        self.type = 'falloff'


class chemically_activated_reaction(falloff_base):
    """ A gas-phase, chemically activated reaction. """

    def __init__(self, equation, kLow, kHigh,
                 efficiencies='', falloff=None, id='', options=()):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kLow:
            The rate coefficient for the forward direction in the low-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param kHigh:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies. The
            efficiency for unspecified species is set to 1.0.
        :param falloff:
            An embedded entry specifying a falloff function. If omitted, a
            unity falloff function (Lindemann form) will be used.
        :param id:
            An optional identification string.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        """
        super().__init__(equation, kLow, kHigh, efficiencies, falloff, id,
                         options)
        self.type = 'chemically-activated'


class pdep_arrhenius(reaction):
    """
    Pressure-dependent rate calculated by interpolating between Arrhenius
    expressions at different pressures.
    """
    def __init__(self, equation, *args, **kwargs):
        """
        :param equation:
            A string specifying the chemical equation.
        :param args:
            Each additional argument is a sequence of four elements specifying
            the pressure and the Arrhenius parameters at that pressure.
        :param kwargs:
            ``id``, ``order``, and ``options`` may be specified as keyword
            arguments, with the same meanings as for class `reaction`.
        """
        super().__init__(equation, None, **kwargs)
        self.arrhenius = args
        self.type = 'pressure-dependent-Arrhenius'

    def get_yaml(self, out):
        super().get_yaml(out)
        rates = []
        for p, A, b, Ea in self.arrhenius:
            rates.append(FlowMap([('P', applyUnits(p)),
                                 ('A', applyUnits(A)),
                                 ('b', applyUnits(b)),
                                 ('Ea', applyUnits(Ea))]))
        out['rate-constants'] = rates


class chebyshev_reaction(reaction):
    """
    Pressure-dependent rate calculated in terms of a bivariate Chebyshev
    polynomial.
    """
    def __init__(self, equation, Tmin=300.0, Tmax=2500.0, Pmin=(0.001, 'atm'),
                 Pmax=(100.0, 'atm'), coeffs=(), **kwargs):
        """
        :param equation:
            A string specifying the chemical equation.
        :param Tmin:
            The minimum temperature at which the rate expression is defined
        :param Tmax:
            the maximum temperature at which the rate expression is defined
        :param Pmin:
            The minimum pressure at which the rate expression is defined
        :param Pmax:
            The maximum pressure at which the rate expression is defined
        :param coeffs:
            A 2D array of the coefficients defining the rate expression. For a
            polynomial with M points in temperature and N points in pressure,
            this should be a list of M lists each with N elements.
        :param kwargs:
            ``id``, ``order``, and ``options`` may be specified as keyword
            arguments, with the same meanings as for class `reaction`.
        """
        super().__init__(equation, None, **kwargs)
        self.type = 'Chebyshev'
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.coeffs = coeffs

    def get_yaml(self, out):
        super().get_yaml(out)
        out['temperature-range'] = FlowList([applyUnits(self.Tmin),
                                             applyUnits(self.Tmax)])
        out['pressure-range'] = FlowList([applyUnits(self.Pmin),
                                          applyUnits(self.Pmax)])
        out['data'] = [FlowList(line) for line in self.coeffs]


class surface_reaction(reaction):
    """
    A heterogeneous chemical reaction with pressure-independent rate
    coefficient and mass-action kinetics.
    """
    def __init__(self, equation, kf, id='', order='', beta=None, options=(),
                 rate_coeff_type=''):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function.
        :param id:
            An optional identification string.
        :param beta:
            Charge transfer coefficient: A number between 0 and 1 which, for a
            charge transfer reaction, determines how much of the electric
            potential difference between two phases is applied to the
            activation energy of the fwd reaction.  The remainder is applied to
            the reverse reaction.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`__.
        """
        super().__init__(equation, kf, id, order, options)
        self.type = 'surface'
        self.sticking = isinstance(kf, stick)
        self.beta = beta
        self.rate_coeff_type = rate_coeff_type

    def get_yaml(self, out):
        super().get_yaml(out)
        if self.sticking:
            del out['rate-constant']
            out.insert(1, 'sticking-coefficient', self.kf)
            if self.kf.motz_wise is not None:
                out['Motz-Wise'] = self.kf.motz_wise
        if self.rate_coeff_type == 'exchangecurrentdensity':
            out['exchange-current-density-formulation'] = True

        if self.kf.coverage is not None:
            cov = {c[0]: FlowMap([('a', c[1]), ('m', c[2]), ('E', c[3])])
                   for c in self.kf.coverage}
            out['coverage-dependencies'] = cov
        if self.beta is not None:
            out['beta'] = self.beta


class edge_reaction(surface_reaction):
    def __init__(self, equation, kf, id='', order='', beta=None, options=(),
                 rate_coeff_type=''):
        super().__init__(equation, kf, id, order, beta, options, rate_coeff_type)
        self.type = 'edge'


class state:
    """
    An embedded entry that specifies the thermodynamic state of a phase
    or interface.
    """
    def __init__(self, temperature=None, pressure=None, mole_fractions=None,
                 mass_fractions=None, density=None, coverages=None,
                 solute_molalities=None):
        """
        :param temperature:
            The temperature.
        :param pressure:
            The pressure.
        :param mole_fractions:
            A string specifying the species mole fractions. Unspecified species
            are set to zero.
        :param mass_fractions:
            A string specifying the species mass fractions. Unspecified species
            are set to zero.
        :param density:
            The density. Cannot be specified if the phase is incompressible.
        :param coverages:
            A string specifying the species coverages. Unspecified species are
            set to zero. Can only be specified for interfaces.
        :param solute_molalities:
            A string specifying the solute molalities. Unspecified molalities
            are set to zero. Only applies to molality-based thermodynamic
            models.
        """
        self.t = temperature
        self.rho = density
        self.p = pressure
        self.X = mole_fractions
        self.Y = mass_fractions
        self.coverages = coverages
        self.molalities = solute_molalities

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        if node.t is not None:
            out['T'] = applyUnits(node.t)
        if node.p is not None:
            out['P'] = applyUnits(node.p)
        if node.rho is not None:
            out['density'] = applyUnits(node.rho)
        if node.X is not None:
            out['X'] = FlowMap(get_composition(node.X).items())
        if node.Y is not None:
            out['Y'] = FlowMap(get_composition(node.Y).items())
        if node.coverages is not None:
            out['coverages'] = FlowMap(get_composition(node.coverages).items())
        if node.molalities is not None:
            out['molalities'] = FlowMap(get_composition(node.molalities).items())
        return representer.represent_dict(out)


class phase:
    """Base class for phases of matter."""

    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', initial_state=None, options=()):
        """
        :param name:
            A string to identify the phase. Must be unique among the phase
            names within the file.
        :param elements:
            The elements. A string of element symbols.
        :param species:
            The species. A string or sequence of strings in the format
            described in `Defining the Species
            <https://cantera.org/tutorials/cti/phases.html#defining-the-species>`__.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param reactions:
            The homogeneous reactions. If omitted, no reactions will be
            included. A string or sequence of strings in the format described
            in `Declaring the Reactions
            <https://cantera.org/tutorials/cti/phases.html#declaring-the-reactions>`__.
            This field is not allowed for ``stoichiometric_solid`` and
            ``stoichiometric_liquid`` entries.
        :param initial_state:
            Initial thermodynamic state, specified with an embedded state entry.
        :param options:
            Special processing options. Optional.
        """

        self.name = name
        self.elements = elements
        self.species = []
        self.reactions = []
        self.kinetics = None
        self.transport = None
        self.comment = note
        self.options = [options] if isinstance(options, str) else options

        #--------------------------------
        #        process species
        #--------------------------------

        # if a single string is entered, make it a list
        if isinstance(species, str):
            species = [species]

        # for each species string, check whether or not the species
        # are imported or defined locally. If imported, the string
        # contains a colon (:)
        for sp in species:
            foundColon = False
            allLocal = True
            for token in sp.split():
                if ':' in sp:
                    foundColon = True
                if token not in _speciesnames:
                    allLocal = False

            if foundColon and not allLocal:
                icolon = sp.find(':')
                datasrc = sp[:icolon].strip()
                spnames = sp[icolon+1:].strip()
                if spnames != 'all':
                    spnames = FlowList(spnames.split())
                self.species.append((datasrc + '.yaml/species', spnames))

            else:
                spnames = sp
                self.species.append(('species', FlowList(spnames.split())))

        if isinstance(reactions, str):
            reactions = [reactions]

        # for each reaction string, check whether or not the reactions
        # are imported or defined locally. If imported, the string
        # contains a colon (:)
        for r in reactions:
            icolon = r.find(':')
            if icolon > 0:
                datasrc = r[:icolon].strip() + '.yaml/reactions'
                rnum = r[icolon+1:].strip()
            else:
                datasrc = 'reactions'
                rnum = r.strip()
            if rnum == 'all' and 'skip_undeclared_species' in self.options:
                rnum = 'declared-species'
            if rnum != 'none':
                self.reactions.append([datasrc, rnum])
            if rnum in ('all', 'declared-species', 'none'):
                continue
            if '*' in rnum:
                if datasrc != 'reactions':
                    _printerr("WARNING: Reaction id-pattern matching from remote"
                        " files not supported ({}: {})".format(datasrc, rnum))
            else:
                _printerr("WARNING: Reaction specification"
                          " '{}' not supported".format(rnum))

        self.initial_state = initial_state

        # add this phase to the global phase list
        _phases.append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['name'] = self.name
        out['thermo'] = self.thermo_model
        if self.elements:
            out['elements'] = FlowList(self.elements.split())

        if len(self.species) == 1 and self.species[0][0] == 'species':
            # all local species
            out['species'] = self.species[0][1]
        else:
            out['species'] = [BlockMap([(sp[0], sp[1])]) for sp in self.species]

        if 'skip_undeclared_elements' in self.options:
            out['skip-undeclared-elements'] = True

        if 'skip_undeclared_third_bodies' in self.options:
            out['skip-undeclared-third-bodies'] = True

        # Convert reaction pattern matching to use of multiple reaction sections
        for i in range(len(self.reactions)):
            spec = self.reactions[i][1]
            name = self.name + '-reactions'
            if '*' in spec and name not in _reactions:
                pattern = re.compile(spec.replace('*', '.*'))
                misses = []
                hits = []
                for reaction in _reactions['reactions']:
                    if pattern.match(reaction.id):
                        hits.append(reaction)
                    else:
                        misses.append(reaction)
                _reactions[name] = hits
                _reactions['reactions'] = misses
                self.reactions[i] = [name, 'all']

        if self.kinetics and self.reactions:
            out['kinetics'] = _newNames[self.kinetics]
            if len(self.reactions) == 1 and self.reactions[0][0] == 'reactions':
                out['reactions'] = self.reactions[0][1]
            elif all(r[1] == 'all' for r in self.reactions):
                out['reactions'] = FlowList(r[0] for r in self.reactions)
            else:
                out['reactions'] = [BlockMap([(r[0], r[1])]) for r in self.reactions]

        if self.transport:
            out['transport'] = _newNames[self.transport]

        if self.comment:
            out['note'] = self.comment

        if self.initial_state:
            out['state'] = self.initial_state


class ideal_gas(phase):
    """An ideal gas mixture."""
    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', kinetics='GasKinetics', transport=None,
                 initial_state=None, options=()):
        """
        The parameters correspond to those of :class:`phase`, with the
        following modifications:

        :param kinetics:
            The kinetics model. Usually this field is omitted, in which case
            kinetics model GasKinetics, appropriate for reactions in ideal gas
            mixtures, is used.
        :param transport:
            The transport property model. One of the strings ``'none'``,
            ``'multi'``, or ``'mix'``. Default: ``'none'``.
        """

        phase.__init__(self, name, elements, species, note, reactions,
                       initial_state, options)
        self.kinetics = kinetics
        self.transport = transport
        self.thermo_model = 'ideal-gas'


class stoichiometric_solid(phase):
    """
    A solid compound or pure element. Stoichiometric solid phases contain
    exactly one species, which always has unit activity. The solid is assumed
    to have constant density. Therefore the rates of reactions involving these
    phases do not contain any concentration terms for the (one) species in the
    phase, since the concentration is always the same."""
    def __init__(self, name='', elements='', species='', note='', density=None,
                 transport='None', initial_state=None, options=()):
        """
        See :class:`phase` for descriptions of the parameters.
        """

        phase.__init__(self, name, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'fixed-stoichiometry'
        self.density = density
        if self.density is None:
            raise InputError('density must be specified.')
        self.transport = None if transport == 'None' else transport

    def get_yaml(self, out):
        super().get_yaml(out)
        for section, names in self.species:
            if section != 'species':
                out['density'] = applyUnits(self.density)
            else:
                species = [S for S in _species if S.name == names[0]][0]
                species.density = self.density


class stoichiometric_liquid(stoichiometric_solid):
    """
    An incompressible stoichiometric liquid. Currently, there is no
    distinction between stoichiometric liquids and solids.
    """


class metal(phase):
    """A metal."""
    def __init__(self, name='', elements='', species='', note='', density=-1.0,
                 transport='None', initial_state=None, options=()):

        phase.__init__(self, name, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'electron-cloud'
        self.density = density

    def get_yaml(self, out):
        super().get_yaml(out)
        out['density'] = applyUnits(self.density)


class incompressible_solid(phase):
    """An incompressible solid."""
    def __init__(self, name='', elements='', species='', note='', density=None,
                 transport='None', initial_state=None, options=()):

        phase.__init__(self, name, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'constant-density'
        self.density = density
        if self.density is None:
            raise InputError('density must be specified.')

    def get_yaml(self, out):
        super().get_yaml(out)
        out['density'] = applyUnits(self.density)


class liquid_vapor(phase):
    """
    A fluid with a complete liquid/vapor equation of state. This entry type
    selects one of a set of predefined fluids with built-in liquid/vapor
    equations of state. The substance_flag parameter selects the fluid. See
    liquidvapor.cti and liquidvapor.py for the usage of this entry type.
    """
    pure_fluids = {
        0: 'water',
        1: 'nitrogen',
        2: 'methane',
        3: 'hydrogen',
        4: 'oxygen',
        5: 'HFC-134a',
        7: 'carbon-dioxide',
        8: 'heptane'
    }

    def __init__(self, name='', elements='', species='', note='',
                 substance_flag=0, initial_state=None, options=()):

        phase.__init__(self, name, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'pure-fluid'
        self.substance_flag = substance_flag

    def get_yaml(self, out):
        super().get_yaml(out)
        if self.substance_flag in self.pure_fluids:
            out['pure-fluid-name'] = self.pure_fluids[self.substance_flag]
        else:
            raise InputError('liquid_vapor: unrecognized value "{}" for '
                '"substance_flag"', self.substance_flag)


class pureFluidParameters:
    def __init__(self, species=None, a_coeff=(), b_coeff=0):
        self.species = species
        self.a_coeff = a_coeff
        self.b_coeff = b_coeff


class crossFluidParameters:
    def __init__(self, species=None, a_coeff=(), b_coeff=()):
        self.species1, self.species2 = species.split(' ')
        self.a_coeff = a_coeff
        self.b_coeff = b_coeff


class RedlichKwongMFTP(phase):
    """
    A multi-component fluid model for non-ideal gas fluids.
    """

    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', kinetics='GasKinetics', initial_state=None,
                 activity_coefficients=None, transport='None', options=()):

        phase.__init__(self,name, elements, species, note, reactions,
                       initial_state,options)
        self.thermo_model = 'Redlich-Kwong'
        self.kinetics = kinetics
        self.transport = None if transport == 'None' else transport
        self.activity_coefficients = activity_coefficients

    def get_yaml(self, out):
        super().get_yaml(out)
        for section, names in self.species:
            if section != 'species':
                _printerr("WARNING: Converting Redlich-Kwong species from"
                    " different input files ({}) is not supported.".format(section))

        spdict = {sp.name: sp for sp in _species}
        for params in self.activity_coefficients:
            if isinstance(params, pureFluidParameters):
                sp = spdict[params.species]
                sp.rk_pure = {'a': params.a_coeff, 'b': params.b_coeff}
            elif isinstance(params, crossFluidParameters):
                sp1 = spdict[params.species1]
                sp1.rk_binary[params.species2] = params.a_coeff
                sp2 = spdict[params.species2]
                sp2.rk_binary[params.species1] = params.a_coeff


class constantIncompressible:
    """Constant molar volume."""
    def __init__(self, molarVolume=0.0):
        """
        :param molarVolume:
            Reference-state molar volume. Default: 0.0.
        """
        self.molar_volume = molarVolume


class IdealSolidSolution(phase):
    """An IdealSolidSolution phase."""
    def __init__(self, name='', elements='', species='', note='',
                 transport='None', initial_state=None,
                 standard_concentration=None, options=()):
        phase.__init__(self, name, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'ideal-condensed'
        self.standard_concentration = standard_concentration
        if self.standard_concentration is None:
            raise InputError('In phase {}: standard_concentration must be specified.', name)
        self.transport = None if transport == 'None' else transport

    def get_yaml(self, out):
        super().get_yaml(out)
        out['standard-concentration-basis'] = _newNames[self.standard_concentration]


class table:
    """User provided thermo table for BinarySolutionTabulatedThermo"""
    def __init__(self, moleFraction=([],''), enthalpy=([],''), entropy=([],'')):
        """
        :param moleFraction:
            The mole fraction of the tabulated species. Required parameter.
        :param enthalpy:
            The enthalpy of the tabulated species. Required parameter.
        :param entropy:
            The entropy of the tabulated species. Required parameter.
        """
        self.x = moleFraction
        self.h = enthalpy
        self.s = entropy


class BinarySolutionTabulatedThermo(IdealSolidSolution):
    """A BinarySolutionTabulatedThermo phase."""
    def __init__(self, name='', elements='', species='', note='',
                 transport='None', initial_state=None,
                 standard_concentration=None, tabulated_species=None,
                 tabulated_thermo=None, options=()):
        """
        The parameters correspond to those of :class:`phase`, with the
        following modifications:

        :param standard_concentration:
            Basis for the standard concentration. One of ``unity``,
            ``molar_volume``, or ``solvent_volume``.
        :param tabulated_species:
            The name of the species for to which the ``tabulated_thermo`` is
            added.
        :param tabulated_thermo:
            An instance of class `table` providing the excess enthalpy
            and entropy to be added to the ``tabulated_species``.

        """
        super().__init__(name, elements, species, note, transport,
                         initial_state, standard_concentration, options)
        self.thermo_model = 'binary-solution-tabulated'
        self.tabulated_species = tabulated_species
        self.tabulated_thermo = tabulated_thermo
        if tabulated_species is None:
            raise InputError('In phase {}: tabulated_species must be specified.', name)
        if tabulated_thermo is None:
            raise InputError('In phase {}: Thermo data must be provided for the tabulated_species.', name)

    def get_yaml(self, out):
        super().get_yaml(out)
        out['tabulated-species'] = self.tabulated_species
        energy_units, quantity_units = self.tabulated_thermo.h[1].split('/')
        tabThermo = BlockMap()
        if energy_units != _uenergy or quantity_units != _umol:
            tabThermo['units'] = FlowMap([('energy', energy_units),
                                          ('quantity', quantity_units)])
        tabThermo['mole-fractions'] = FlowList(self.tabulated_thermo.x[0])
        tabThermo['enthalpy'] = FlowList(self.tabulated_thermo.h[0])
        tabThermo['entropy'] = FlowList(self.tabulated_thermo.s[0])
        out['tabulated-thermo'] = tabThermo

class lattice(phase):
    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', transport='None', initial_state=None,
                 options=(), site_density=None):
        phase.__init__(self, name, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'lattice'
        self.site_density = site_density

        if name == '':
            raise InputError('sublattice name must be specified')
        if species == '':
            raise InputError('sublattice species must be specified')
        if site_density is None:
            raise InputError('sublattice '+name
                            +' site density must be specified')

    def get_yaml(self, out):
        super().get_yaml(out)
        out['site-density'] = applyUnits(self.site_density)

class ideal_interface(phase):
    """A chemically-reacting ideal surface solution of multiple species."""
    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', site_density=0.0, phases=(),
                 kinetics='Interface', transport='None', initial_state=None,
                 options=()):
        """
        The parameters correspond to those of :class:`phase`, with the
        following modifications:

        :param reactions:
            The heterogeneous reactions at this interface. If omitted, no
            reactions will be included. A string or sequence of strings in the
            format described in `Declaring the Reactions
            <https://cantera.org/tutorials/cti/phases.html#declaring-the-reactions>`__.
        :param site_density:
            The number of adsorption sites per unit area.
        :param phases:
            A string listing the bulk phases that participate in reactions
            at this interface.
        """
        phase.__init__(self, name, elements, species, note, reactions,
                       initial_state, options)
        self.thermo_model = 'ideal-surface'
        self.kinetics = kinetics
        self.transport = None if transport == 'None' else transport
        self.site_density = site_density

    def get_yaml(self, out):
        super().get_yaml(out)
        out['site-density'] = applyUnits(self.site_density)
        if _motz_wise is not None:
            out['Motz-Wise'] = _motz_wise


class edge(ideal_interface):
    """A 1D boundary between two surface phases."""
    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', site_density=0.0, phases=(), kinetics='Edge',
                 transport='None', initial_state=None, options=()):

        ideal_interface.__init__(self, name, elements, species, note, reactions,
            site_density, phases, kinetics, transport, initial_state, options)
        self.thermo_model = 'edge'


# Falloff parameterizations

class Troe:
    """The Troe falloff function."""
    def __init__(self, A=0.0, T3=0.0, T1=0.0, T2=None):
        """
        Parameters: *A*, *T3*, *T1*, *T2*. These must be entered as pure
        numbers with no attached dimensions.
        """
        self.A = A
        self.T3 = T3
        self.T1 = T1
        self.T2 = T2

    def get_yaml(self, out):
        troe = FlowMap([('A', self.A), ('T3', self.T3), ('T1', self.T1)])
        if self.T2 is not None:
            troe['T2'] = self.T2
        out['Troe'] = troe


class SRI:
    """ The SRI falloff function."""
    def __init__(self, A=0.0, B=0.0, C=0.0, D=None, E=None):
        """
        Parameters: *A*, *B*, *C*, *D*, *E*. These must be entered as
        pure numbers without attached dimensions.
        """
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def get_yaml(self, out):
        sri = FlowMap([('A', self.A), ('B', self.B), ('C', self.C)])
        if self.D is not None:
            sri['D'] = self.D
        if self.E is not None:
            sri['E'] = self.E
        out['SRI'] = sri


class Lindemann:
    """The Lindemann falloff function."""
    def get_yaml(self, out):
        pass


def convert(filename=None, output_name=None, text=None):
    # Reset global state, in case cti2yaml is being used as a module and convert
    # is being called multiple times.
    units('m', 'kmol', 'kg', 's', 'J/kmol', 'J', 'Pa')
    standard_pressure(OneAtm)
    global _motz_wise
    _motz_wise = None
    _elements.clear()
    _species.clear()
    _speciesnames.clear()
    _phases.clear()
    _reactions.clear()
    _reactions['reactions'] = []

    if filename is None and text is None:
        raise ValueError("One of filename or text must be specified")
    elif filename is not None and text is not None:
        raise ValueError("Only one of filename or text should be specified")

    if filename is not None:
        filename = pathlib.Path(filename).expanduser()
        base = filename.name
        root = filename.stem
        dataset(root)

    if output_name is None and _name != 'noname':
        output_name = pathlib.Path(_name + '.yaml')
    else:
        output_name = pathlib.Path(output_name)

    try:
        if filename is not None:
            with filename.open('r', encoding='latin-1') as f:
                code = compile(f.read(), str(filename), 'exec')
        else:
            code = compile(text, '<string>', 'exec')
        exec(code)
    except SyntaxError as err:
        # Show more context than the default SyntaxError message
        # to help see problems in multi-line statements
        if filename:
            text = pathlib.Path(filename).read_text()
        text = text.split('\n')
        _printerr('{} in "{}" on line {}:\n'.format(
            err.__class__.__name__, err.filename, err.lineno))
        _printerr('|  Line |')
        for i in range(max(err.lineno-6, 0),
                       min(err.lineno+3, len(text))):
            _printerr('| {: 5d} |'.format(i+1), text[i].rstrip())
            if i == err.lineno-1:
                _printerr(' '* (err.offset+9) + '^')
        _printerr()
        sys.exit(3)
    except Exception as err:
        import traceback

        if filename:
            text = pathlib.Path(filename).read_text()
        else:
            filename = '<string>'
        text = text.split('\n')
        tb = traceback.extract_tb(sys.exc_info()[2])
        lineno = tb[-1][1]
        if tb[-1][0] == filename:
            # Error in input file
            _printerr('{} on line {} of {}:'.format(
                err.__class__.__name__, lineno, filename))
            _printerr(err)
            _printerr('\n| Line |')

            for i in range(max(lineno-6, 0),
                           min(lineno+3, len(text))):
                if i == lineno-1:
                    _printerr('> {: 4d} >'.format(i+1), text[i].rstrip())
                else:
                    _printerr('| {: 4d} |'.format(i+1), text[i].rstrip())
        else:
            # Error in cti2yaml or elsewhere
            traceback.print_exc()

        sys.exit(4)

    # write the YAML file
    emitter = yaml.YAML()
    emitter.width = 70

    for name, cls in globals().items():
        if hasattr(cls, 'to_yaml'):
            emitter.register_class(cls)

    with output_name.open('w') as dest:
        # information regarding conversion
        metadata = BlockMap([
            ('generator', 'cti2yaml'),
            ('cantera-version', '2.5.0a4'),
            ('date', formatdate(localtime=True)),
        ])
        if filename is not None:
            metadata['input-files'] = FlowList([base])
        emitter.dump(metadata, dest)

        out_units = FlowMap([])
        if _umass != 'kg':
            out_units['mass'] = _umass
        if _ulen != 'm':
            out_units['length'] = _ulen
        if _utime != 's':
            out_units['time'] = _utime
        if _upres != 'Pa':
            out_units['pressure'] = _upres
        if _uenergy != 'J':
            out_units['energy'] = _uenergy
        if _umol != 'kmol':
            out_units['quantity'] = _umol
        if _ue != 'J/kmol':
            out_units['activation-energy'] = _ue

        if out_units:
            units_map = BlockMap([('units', out_units)])
            units_map.yaml_set_comment_before_after_key('units', before='\n')
            emitter.dump(units_map, dest)

        if _elements:
            elements_map = BlockMap([('elements', _elements)])
            elements_map.yaml_set_comment_before_after_key('elements', before='\n')
            emitter.dump(elements_map, dest)

        if _phases:
            phases_map = BlockMap([('phases', _phases)])
            phases_map.yaml_set_comment_before_after_key('phases', before='\n')
            emitter.dump(phases_map, dest)

        if _species:
            species_map = BlockMap([('species', _species)])
            species_map.yaml_set_comment_before_after_key('species', before='\n')
            emitter.dump(species_map, dest)

        for name, reactions in _reactions.items():
            if reactions:
                reactions_map = BlockMap([(name, reactions)])
                reactions_map.yaml_set_comment_before_after_key(name, before='\n')
                emitter.dump(reactions_map, dest)


def main():
    """Parse command line arguments and pass them to `convert`."""
    parser = argparse.ArgumentParser(
        description="Convert legacy CTI input files to YAML format",
        epilog=("The 'output' argument is optional. If it is not given, an output "
                "file with the same name as the input file is used, with the extension "
                "changed to '.yaml'.")
    )
    parser.add_argument("input", help="The input CTI filename. Must be specified.")
    parser.add_argument("output", nargs="?", help="The output YAML filename. Optional.")
    if len(sys.argv) not in [2, 3]:
        if len(sys.argv) > 3:
            print(
                "cti2yaml.py: error: unrecognized arguments:",
                ' '.join(sys.argv[3:]),
                file=sys.stderr,
            )
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_file = pathlib.Path(args.input)
    if args.output is None:
        output_file = input_file.with_suffix(".yaml")
    else:
        output_file = pathlib.Path(args.output)

    convert(input_file, output_file)


if __name__ == "__main__":
    main()
