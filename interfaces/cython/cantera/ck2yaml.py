#!/usr/bin/env python
# encoding: utf-8

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
ck2yaml.py: Convert Chemkin-format mechanisms to Cantera YAML input files

Usage:
    ck2yaml [--input=<filename>]
            [--thermo=<filename>]
            [--transport=<filename>]
            [--surface=<filename>]
            [--name=<name>]
            [--extra=<filename>]
            [--output=<filename>]
            [--permissive]
            [--quiet]
            [--no-validate]
            [-d | --debug]

Example:
    ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.yaml'.

An input file containing only species definitions (which can be referenced from
phase definitions in other input files) can be created by specifying only a
thermo file.

For the case of a surface mechanism, the gas phase input file should be
specified as 'input' and the surface phase input file should be specified as
'surface'.

The '--permissive' option allows certain recoverable parsing errors (e.g.
duplicate transport data) to be ignored. The '--name=<name>' option
is used to override default phase names (i.e. 'gas').

The '--extra=<filename>' option takes a YAML file as input. This option can be
used to add to the file description, or to define custom fields that are
included in the YAML output.
"""

from collections import defaultdict, OrderedDict
import logging
import os.path
import sys
import numpy as np
import re
import itertools
import getopt
import textwrap
from email.utils import formatdate

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

BlockMap = yaml.comments.CommentedMap

logger = logging.getLogger(__name__)
loghandler = logging.StreamHandler(sys.stdout)
logformatter = logging.Formatter('%(message)s')
loghandler.setFormatter(logformatter)
logger.handlers.clear()
logger.addHandler(loghandler)
logger.setLevel(logging.INFO)

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

QUANTITY_UNITS = {'MOL': 'mol',
                  'MOLE': 'mol',
                  'MOLES': 'mol',
                  'MOLEC': 'molec',
                  'MOLECULES': 'molec'}

ENERGY_UNITS = {'CAL/': 'cal/mol',
                'CAL/MOL': 'cal/mol',
                'CAL/MOLE': 'cal/mol',
                'EVOL': 'eV',
                'EVOLTS': 'eV',
                'JOUL': 'J/mol',
                'JOULES/MOL': 'J/mol',
                'JOULES/MOLE': 'J/mol',
                'KCAL': 'kcal/mol',
                'KCAL/MOL': 'kcal/mol',
                'KCAL/MOLE': 'kcal/mol',
                'KELV': 'K',
                'KELVIN': 'K',
                'KELVINS': 'K',
                'KJOU': 'kJ/mol',
                'KJOULES/MOL': 'kJ/mol',
                'KJOULES/MOLE': 'kJ/mol'}

def strip_nonascii(s):
    return s.encode('ascii', 'ignore').decode()


def compatible_quantities(quantity_basis, units):
    if quantity_basis == 'mol':
        return 'molec' not in units
    elif quantity_basis == 'molec':
        return 'molec' in units or 'mol' not in units
    else:
        raise ValueError('Unknown quantity basis: "{}"'.format(quantity_basis))


class InputError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin-format
    mechanism files. Pass a string describing the circumstances that caused
    the exceptional behavior.
    """
    def __init__(self, message, *args, **kwargs):
        if args or kwargs:
            super().__init__(message.format(*args, **kwargs))
        else:
            super().__init__(message)


class Species:
    def __init__(self, label, sites=None):
        self.label = label
        self.thermo = None
        self.transport = None
        self.sites = sites
        self.composition = None
        self.note = None

    def __str__(self):
        return self.label

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('name', node.label),
                        ('composition', FlowMap(node.composition.items()))])
        if node.thermo:
            out['thermo'] = node.thermo
        if node.transport:
            out['transport'] = node.transport
        if node.sites:
            out['sites'] = node.sites
        if node.note:
            out['note'] = node.note
        return representer.represent_dict(out)


class Nasa7:
    """
    Thermodynamic data parameterized as two seven-coefficient NASA
    polynomials.
    See https://cantera.org/science/science-species.html#the-nasa-7-coefficient-polynomial-parameterization
    """
    def __init__(self, *, Tmin, Tmax, Tmid, low_coeffs, high_coeffs, note=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Tmid = Tmid
        self.low_coeffs = low_coeffs
        self.high_coeffs = high_coeffs
        self.note = note

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('model', 'NASA7')])
        out['temperature-ranges'] = FlowList([node.Tmin, node.Tmid, node.Tmax])
        out['data'] = [FlowList(node.low_coeffs), FlowList(node.high_coeffs)]
        if node.note:
            note = textwrap.dedent(node.note.rstrip())
            if '\n' in note:
                note = yaml.scalarstring.PreservedScalarString(note)
            out['note'] = note
        return representer.represent_dict(out)


class Nasa9:
    """
    Thermodynamic data parameterized as any number of nine-coefficient NASA
    polynomials.
    See https://cantera.org/science/science-species.html#the-nasa-9-coefficient-polynomial-parameterization

    :param data:
        List of polynomials, where each polynomial is written as
        ```
        [(T_low, T_high), [a_0, a_1, ..., a_8]]
        ```
    """
    def __init__(self, *, data, note=''):
        self.note = note
        self.data = list(sorted(data))
        self.Tranges = [self.data[0][0][0]]
        for i in range(1, len(data)):
            if abs(self.data[i-1][0][1] - self.data[i][0][0]) > 0.01:
                raise ValueError('NASA9 polynomials contain non-adjacent temperature ranges')
            self.Tranges.append(self.data[i][0][0])
        self.Tranges.append(self.data[-1][0][1])

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('model', 'NASA9')])
        out['temperature-ranges'] = FlowList(node.Tranges)
        out['data'] = [FlowList(poly) for (trange, poly) in node.data]
        if node.note:
            out['note'] = node.note
        return representer.represent_dict(out)


class Reaction:
    """
    :param index:
        A unique nonnegative integer index
    :param reactants:
        A list of `(stoichiometry, species name)` tuples
    :param products:
        A list of `(stoichiometry, species name)` tuples
    :param kinetics:
        A `KineticsModel` instance which describes the rate constant
    :param reversible:
        Boolean indicating whether the reaction is reversible
    :param duplicate:
        Boolean indicating whether the reaction is a known (permitted) duplicate
    :param forward_orders:
        A dictionary specifying a non-default reaction order (value) for each
        specified species (key)
    :param third_body:
        A string name used for the third-body species written in
        pressure-dependent reaction types (usually "M")
    """

    def __init__(self, parser, index=-1, reactants=None, products=None,
                 kinetics=None, reversible=True, duplicate=False,
                 forward_orders=None, third_body=None):
        self.parser = parser
        self.index = index
        self.reactants = reactants  # list of (stoichiometry, species) tuples
        self.products = products  # list of (stoichiometry, species) tuples
        self.kinetics = kinetics
        self.reversible = reversible
        self.duplicate = duplicate
        self.forward_orders = forward_orders or {}
        self.third_body = ''
        self.comment = ''

    def _coeff_string(self, coeffs):
        L = []
        for stoichiometry, species in coeffs:
            if stoichiometry != 1:
                L.append('{0} {1}'.format(stoichiometry, species))
            else:
                L.append(str(species))
        expression = ' + '.join(L)
        expression += self.kinetics.reaction_string_suffix(self.third_body)
        return expression

    def __str__(self):
        """
        Return a string representation of the reaction, e.g. 'A + B <=> C + D'.
        """
        return '{}{}{}'.format(self._coeff_string(self.reactants),
                               ' <=> ' if self.reversible else ' => ',
                               self._coeff_string(self.products))

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('equation', str(node))])
        out.yaml_add_eol_comment('Reaction {}'.format(node.index), 'equation')
        if node.duplicate:
            out['duplicate'] = True
        node.kinetics.reduce(out)
        if node.forward_orders:
            out['orders'] = FlowMap(node.forward_orders)
            if any((float(x) < 0 for x in node.forward_orders.values())):
                out['negative-orders'] = True
                node.parser.warn('Negative reaction order for reaction {} ({}).'.format(
                    node.index, str(node)))
            reactant_names = {r[1].label for r in node.reactants}
            if any((species not in reactant_names for species in node.forward_orders)):
                out['nonreactant-orders'] = True
                node.parser.warn('Non-reactant order for reaction {} ({}).'.format(
                    node.index, str(node)))
        if node.comment:
            comment = textwrap.dedent(node.comment.rstrip())
            if '\n' in comment:
                comment = yaml.scalarstring.PreservedScalarString(comment)
            out['note'] = comment
        return representer.represent_dict(out)


class KineticsModel:
    """
    A base class for kinetics models
    """
    pressure_dependent = None  # overloaded in derived classes

    def __init__(self):
        self.efficiencies = {}

    def reaction_string_suffix(self, species):
        """
        Suffix for reactant and product strings, used for pressure-dependent
        reactions
        """
        return ''

    def reduce(self, output):
        """
        Assign data from this object to the YAML mapping ``output``
        """
        raise InputError('reduce is not implemented for objects of class {}',
                         self.__class__.__name__)


class Arrhenius:
    """
    Represent a modified Arrhenius rate.

    :param A:
        The pre-exponential factor, given as a tuple consisting of a floating
        point value and a units string
    :param b:
        The temperature exponent
    :param Ea:
        The activation energy, given as a tuple consisting of a floating
        point value and a units string
    """
    def __init__(self, A=0.0, b=0.0, Ea=0.0, *, parser):
        self.A = A
        self.b = b
        self.Ea = Ea
        self.parser = parser

    def as_yaml(self, extra=()):
        out = FlowMap(extra)
        if compatible_quantities(self.parser.output_quantity_units, self.A[1]):
            out['A'] = self.A[0]
        else:
            out['A'] = "{0:e} {1}".format(*self.A)

        out['b'] = self.b

        if self.Ea[1] == self.parser.output_energy_units:
            out['Ea'] = self.Ea[0]
        else:
            out['Ea'] = "{0} {1}".format(*self.Ea)

        return out


class ElementaryRate(KineticsModel):
    """
    A reaction rate described by a single Arrhenius expression.
    See https://cantera.org/science/reactions.html#reactions-with-a-pressure-independent-rate

    :param rate:
        The Arrhenius expression describing this reaction rate.
    """
    pressure_dependent = False

    def __init__(self, rate, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.rate = rate

    def reduce(self, output):
        output['rate-constant'] = self.rate.as_yaml()
        if self.rate.A[0] < 0:
            output['negative-A'] = True

class SurfaceRate(KineticsModel):
    """
    An Arrhenius-like reaction occurring on a surface
    See https://cantera.org/science/reactions.html#surface-reactions

    :param rate:
        The Arrhenius expression describing this reaction rate.
    :param coverages:
        A list of tuples where each tuple specifies the coverage dependencies
        for a species, in the form `(species_name, a_k, m_k, E_k)`
    :param is_sticking:
        True if the Arrhenius expression is a parameterization of a sticking
        coefficient, rather than the rate constant itself.
    :param motz_wise:
        True if the sticking coefficient should be translated into a rate
        coefficient using the correction factor developed by Motz & Wise for
        reactions with high (near-unity) sticking coefficients
    """
    pressure_dependent = False

    def __init__(self, *, rate, coverages, is_sticking, motz_wise, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.rate = rate
        self.coverages = coverages
        self.is_sticking = is_sticking
        self.motz_wise = motz_wise

    def reduce(self, output):
        if self.is_sticking:
            output['sticking-coefficient'] = self.rate.as_yaml()
        else:
            output['rate-constant'] = self.rate.as_yaml()

        if self.motz_wise is not None:
            output['Motz-Wise'] = self.motz_wise

        if self.coverages:
            covdeps = BlockMap()
            for species,A,m,E in self.coverages:
                # Energy units for coverage modification match energy units for
                # base reaction
                if self.rate.Ea[1] != self.rate.parser.output_energy_units:
                    E = '{} {}'.format(E, self.rate.Ea[1])
                    covdeps[species] = FlowList([A, m, E])
            output['coverage-dependencies'] = covdeps


class PDepArrhenius(KineticsModel):
    """
    A rate calculated by interpolating between Arrhenius expressions at
    various pressures.
    See https://cantera.org/science/reactions.html#pressure-dependent-arrhenius-rate-expressions-p-log

    :param pressures:
        A list of pressures at which Arrhenius expressions are given.
    :param pressure_units:
        A string indicating the units used for the pressures
    :param arrhenius:
        A list of `Arrhenius` objects at each given pressure
    """
    pressure_dependent = True

    def __init__(self, *, pressures, pressure_units, arrhenius, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.pressures = pressures
        self.pressure_units = pressure_units
        self.arrhenius = arrhenius or []

    def reduce(self, output):
        output['type'] = 'pressure-dependent-Arrhenius'
        rates = []
        for pressure, arrhenius in zip(self.pressures, self.arrhenius):
            rates.append(arrhenius.as_yaml(
                [('P', '{0} {1}'.format(pressure, self.pressure_units))]))
        output['rate-constants'] = rates


class Chebyshev(KineticsModel):
    """
    A rate calculated in terms of a bivariate Chebyshev polynomial.
    See https://cantera.org/science/reactions.html#chebyshev-reaction-rate-expressions

    :param coeffs:
        Matrix of Chebyshev coefficients, dimension N_T by N_P
    :param Tmin:
        Minimum temperature for which the parameterization is valid
    :param Tmax:
        Maximum temperature for which the parameterization is valid
    :param Pmin:
        Minimum pressure for which the parameterization is valid, given as a
        `(value, units)` tuple
    :param Pmax:
        Maximum pressure for which the parameterization is valid, given as a
        `(value, units)` tuple
    :param quantity_units:
        Quantity units for the rate constant
    """
    pressure_dependent = True

    def __init__(self, coeffs, *, Tmin, Tmax, Pmin, Pmax, quantity_units,
                 **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.coeffs = coeffs
        self.quantity_units = quantity_units

    def reaction_string_suffix(self, species):
        return ' (+{})'.format(species if species else 'M')

    def reduce(self, output):
        output['type'] = 'Chebyshev'
        output['temperature-range'] = FlowList([self.Tmin, self.Tmax])
        output['pressure-range'] = FlowList(['{0} {1}'.format(*self.Pmin),
                                             '{0} {1}'.format(*self.Pmax)])
        if self.quantity_units is not None:
            output['units'] = FlowMap([('quantity', self.quantity_units)])
        output['data'] = [FlowList(float(v) for v  in row) for row in self.coeffs]


class ThreeBody(KineticsModel):
    """
    A rate calculated for a reaction which includes a third-body collider.
    See https://cantera.org/science/reactions.html#three-body-reactions

    :param high_rate:
        The Arrhenius kinetics (high-pressure limit)
    :param efficiencies:
        A mapping of species names to collider efficiencies
    """
    pressure_dependent = True

    def __init__(self, high_rate=None, efficiencies=None, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.high_rate = high_rate
        self.efficiencies = efficiencies or {}

    def reaction_string_suffix(self, species):
        return ' + M'

    def reduce(self, output):
        output['type'] = 'three-body'
        output['rate-constant'] = self.high_rate.as_yaml()
        if self.high_rate.A[0] < 0:
            output['negative-A'] = True
        if self.efficiencies:
            output['efficiencies'] = FlowMap(self.efficiencies)


class Falloff(ThreeBody):
    """
    A rate for a pressure-dependent falloff reaction.
    See https://cantera.org/science/reactions.html#falloff-reactions

    :param low_rate:
        The Arrhenius kinetics at the low-pressure limit
    :param high_rate:
        The Arrhenius kinetics at the high-pressure limit
    :param efficiencies:
        A mapping of species names to collider efficiencies
    :param F:
        Falloff function parameterization
    """
    def __init__(self, low_rate=None, F=None, **kwargs):
        ThreeBody.__init__(self, **kwargs)
        self.low_rate = low_rate
        self.F = F

    def reaction_string_suffix(self, species):
        return ' (+{})'.format(species)

    def reduce(self, output):
        output['type'] = 'falloff'
        output['low-P-rate-constant'] = self.low_rate.as_yaml()
        output['high-P-rate-constant'] = self.high_rate.as_yaml()
        if self.high_rate.A[0] < 0 and self.low_rate.A[0] < 0:
            output['negative-A'] = True
        if self.F:
            self.F.reduce(output)
        if self.efficiencies:
            output['efficiencies'] = FlowMap(self.efficiencies)


class ChemicallyActivated(ThreeBody):
    """
    A rate for a chemically-activated reaction.
    See https://cantera.org/science/reactions.html#chemically-activated-reactions

    :param low_rate:
        The Arrhenius kinetics at the low-pressure limit
    :param high_rate:
        The Arrhenius kinetics at the high-pressure limit
    :param efficiencies:
        A mapping of species names to collider efficiencies
    :param F:
        Falloff function parameterization
    """
    def __init__(self, low_rate=None, F=None, **kwargs):
        ThreeBody.__init__(self, **kwargs)
        self.low_rate = low_rate
        self.F = F

    def reaction_string_suffix(self, species):
        return ' (+{})'.format(species)

    def reduce(self, output):
        output['type'] = 'chemically-activated'
        output['low-P-rate-constant'] = self.low_rate.as_yaml()
        output['high-P-rate-constant'] = self.high_rate.as_yaml()
        if self.high_rate.A[0] < 0 and self.low_rate.A[0] < 0:
            output['negative-A'] = True
        if self.F:
            self.F.reduce(output)
        if self.efficiencies:
            output['efficiencies'] = FlowMap(self.efficiencies)


class Troe:
    """
    The Troe falloff function, described with either 3 or 4 parameters.
    See https://cantera.org/science/reactions.html#the-troe-falloff-function
    """
    def __init__(self, A=0.0, T3=0.0, T1=0.0, T2=None):
        self.A = A
        self.T3 = T3
        self.T1 = T1
        self.T2 = T2

    def reduce(self, output):
        troe = FlowMap([('A', self.A), ('T3', self.T3), ('T1', self.T1)])
        if self.T2 is not None:
            troe['T2'] = self.T2
        output['Troe'] = troe


class Sri:
    """
    The SRI falloff function, described with either 3 or 5 parameters.
    See https://cantera.org/science/reactions.html#the-sri-falloff-function
    """
    def __init__(self, *, A, B, C, D=None, E=None):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def reduce(self, output):
        sri = FlowMap([('A', self.A), ('B', self.B), ('C', self.C)])
        if self.D:
            sri['D'] = self.D
        if self.E:
            sri['E'] = self.E

        output['SRI'] = sri


class TransportData:
    geometry_flags = ['atom', 'linear', 'nonlinear']

    def __init__(self, label, geometry, well_depth, collision_diameter,
                 dipole_moment, polarizability, z_rot, note=''):

        try:
            geometry = int(geometry)
        except ValueError:
            raise InputError(
                "Bad geometry flag '{}' for species '{}', is the flag a float "
                "or character? It should be an integer.", geometry, label)
        if geometry not in (0, 1, 2):
            raise InputError("Bad geometry flag '{}' for species '{}'",
                             geometry, label)

        self.geometry = self.geometry_flags[int(geometry)]
        self.well_depth = float(well_depth)
        self.collision_diameter = float(collision_diameter)
        self.dipole_moment = float(dipole_moment)
        self.polarizability = float(polarizability)
        self.z_rot = float(z_rot)
        self.note = note.strip()

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('model', 'gas'),
                        ('geometry', node.geometry),
                        ('well-depth', node.well_depth),
                        ('diameter', node.collision_diameter)])
        if node.dipole_moment:
            out['dipole'] = node.dipole_moment
        if node.polarizability:
            out['polarizability'] = node.polarizability
        if node.z_rot:
            out['rotational-relaxation'] = node.z_rot
        if node.note:
            out['note'] = node.note
        return representer.represent_dict(out)


def fortFloat(s):
    """
    Convert a string representation of a floating point value to a float,
    allowing for some of the peculiarities of allowable Fortran representations.
    """
    return float(s.strip().lower().replace('d', 'e').replace('e ', 'e+'))


def get_index(seq, value):
    """
    Find the first location in *seq* which contains a case-insensitive,
    whitespace-insensitive match for *value*. Returns *None* if no match is
    found.
    """
    if isinstance(seq, str):
        seq = seq.split()
    value = value.lower().strip()
    for i, item in enumerate(seq):
        if item.lower() == value:
            return i
    return None


def contains(seq, value):
    if isinstance(seq, str):
        return value.lower() in seq.lower()
    else:
        return get_index(seq, value) is not None


class Surface:
    def __init__(self, name, site_density):
        self.name = name
        self.site_density = site_density
        self.species_list = []
        self.reactions = []


class Parser:
    def __init__(self):
        self.processed_units = False
        self.energy_units = 'cal/mol'  # for the current REACTIONS section
        self.output_energy_units = 'cal/mol'  # for the output file
        self.quantity_units = 'mol'  # for the current REACTIONS section
        self.output_quantity_units = 'mol'  # for the output file
        self.motz_wise = None
        self.warning_as_error = True

        self.elements = []
        self.element_weights = {}  # for custom elements only
        self.species_list = []  # bulk species only
        self.species_dict = {}  # bulk and surface species
        self.surfaces = []
        self.reactions = []
        self.header_lines = []
        self.extra = {}  # for extra entries
        self.files = []  # input file names

    def warn(self, message):
        if self.warning_as_error:
            raise InputError(message)
        else:
            logger.warning(message)

    @staticmethod
    def parse_composition(elements, nElements, width):
        """
        Parse the elemental composition from a 7 or 9 coefficient NASA polynomial
        entry.
        """
        composition = {}
        for i in range(nElements):
            symbol = elements[width*i:width*i+2].strip()
            count = elements[width*i+2:width*i+width].strip()
            if not symbol:
                continue
            try:
                # Convert to float first for cases where ``count`` is a string
                # like "2.00".
                count = int(float(count))
                if count:
                    composition[symbol.capitalize()] = count
            except ValueError:
                pass
        return composition

    @staticmethod
    def get_rate_constant_units(length_dims, length_units, quantity_dims,
                                quantity_units, time_dims=1, time_units='s'):

        units = ''
        if length_dims:
            units += length_units
        if length_dims > 1:
            units += '^' + str(length_dims)
        if quantity_dims:
            units += '/' + quantity_units
        if quantity_dims > 1:
            units += '^' + str(quantity_dims)
        if time_dims:
            units += '/' + time_units
        if time_dims > 1:
            units += '^' + str(time_dims)
        if units.startswith('/'):
            units = '1' + units
        return units

    def add_element(self, element_string):
        if '/' in element_string:
            name, weight, _ = element_string.split('/')
            weight = fortFloat(weight)
            name = name.capitalize()
            self.elements.append(name)
            self.element_weights[name] = weight
        else:
            self.elements.append(element_string.capitalize())

    def read_NASA7_entry(self, lines, TintDefault, comments):
        """
        Read a thermodynamics entry for one species in a Chemkin-format file
        (consisting of two 7-coefficient NASA polynomials). Returns the label of
        the species, the thermodynamics model as a :class:`Nasa7` object, and
        the elemental composition of the species.

        For more details on this format, see `Debugging common errors in CK files
        <https://cantera.org/tutorials/ck2cti-tutorial.html#debugging-common-errors-in-ck-files>`__.
        """
        identifier = lines[0][0:24].split()
        species = identifier[0].strip()

        if len(identifier) > 1:
            note = ''.join(identifier[1:]).strip()
        else:
            note = ''

        comments = '\n'.join(c.rstrip() for c in comments if c.strip())
        if comments and note:
            note = '\n'.join((note, comments))
        elif comments:
            note = comments

        # Normal method for specifying the elemental composition
        composition = self.parse_composition(lines[0][24:44], 4, 5)

        # Chemkin-style extended elemental composition: additional lines
        # indicated by '&' continuation character on preceding lines. Element
        # names and abundances are separated by whitespace (not fixed width)
        if lines[0].rstrip().endswith('&'):
            complines = []
            for i in range(len(lines)-1):
                if lines[i].rstrip().endswith('&'):
                    complines.append(lines[i+1])
                else:
                    break
            lines = [lines[0]] + lines[i+1:]
            comp = ' '.join(line.rstrip('&\n') for line in complines).split()
            composition = {}
            for i in range(0, len(comp), 2):
                composition[comp[i].capitalize()] = int(comp[i+1])

        # Non-standard extended elemental composition data may be located beyond
        # column 80 on the first line of the thermo entry
        if len(lines[0]) > 80:
            elements = lines[0][80:]
            composition2 = self.parse_composition(elements, len(elements)//10, 10)
            composition.update(composition2)

        if not composition:
            raise InputError("Error parsing elemental composition for "
                             "species '{}'", species)

        # Extract the NASA polynomial coefficients
        # Remember that the high-T polynomial comes first!
        Tmin = fortFloat(lines[0][45:55])
        Tmax = fortFloat(lines[0][55:65])
        try:
            Tint = fortFloat(lines[0][65:75])
        except ValueError:
            Tint = TintDefault

        high_coeffs = [fortFloat(lines[i][j:k])
                       for i,j,k in [(1,0,15), (1,15,30), (1,30,45), (1,45,60),
                                     (1,60,75), (2,0,15), (2,15,30)]]
        low_coeffs = [fortFloat(lines[i][j:k])
                       for i,j,k in [(2,30,45), (2,45,60), (2,60,75), (3,0,15),
                                     (3,15,30), (3,30,45), (3,45,60)]]

        # Duplicate the valid set of coefficients if only one range is provided
        if all(c == 0 for c in low_coeffs) and Tmin == Tint:
            low_coeffs = high_coeffs
        elif all(c == 0 for c in high_coeffs) and Tmax == Tint:
            high_coeffs = low_coeffs

        # Construct and return the thermodynamics model
        thermo = Nasa7(Tmin=Tmin, Tmax=Tmax, Tmid=Tint,
                           low_coeffs=low_coeffs, high_coeffs=high_coeffs,
                           note=note)

        return species, thermo, composition

    def read_NASA9_entry(self, entry, comments):
        """
        Read a thermodynamics ``entry`` for one species given as one or more
        9-coefficient NASA polynomials, written in the format described in
        Appendix A of NASA Reference Publication 1311 (McBride and Gordon, 1996).
        Returns the label of the species, the thermodynamics model as a
        :class:`Nasa9` object, and the elemental composition of the species
        """
        tokens = entry[0].split()
        species = tokens[0]
        note = ' '.join(tokens[1:])
        N = int(entry[1][:2])
        note2 = entry[1][3:9].strip()
        if note and note2:
            note = '{0} [{1}]'.format(note, note2)
        elif note2:
            note = note2

        comments = '\n'.join(c.rstrip() for c in comments if c.strip())
        if comments and note:
            note = '\n'.join((note, comments))
        elif comments:
            note = comments

        composition = self.parse_composition(entry[1][10:50], 5, 8)

        polys = []
        try:
            for i in range(N):
                A, B, C = entry[2+3*i:2+3*(i+1)]
                Trange = [fortFloat(A[1:11]), fortFloat(A[11:21])]
                coeffs = [fortFloat(B[0:16]), fortFloat(B[16:32]),
                          fortFloat(B[32:48]), fortFloat(B[48:64]),
                          fortFloat(B[64:80]), fortFloat(C[0:16]),
                          fortFloat(C[16:32]), fortFloat(C[48:64]),
                          fortFloat(C[64:80])]
                polys.append((Trange, coeffs))
        except (IndexError, ValueError) as err:
            raise InputError('Error while reading thermo entry for species {}:\n{}',
                             species, err)

        thermo = Nasa9(data=polys, note=note)

        return species, thermo, composition

    def setup_kinetics(self):
        # We look for species including the next permissible character. '\n' is
        # appended to the reaction string to identify the last species in the
        # reaction string. Checking this character is necessary to correctly
        # identify species with names ending in '+' or '='.
        self.species_tokens = set()
        for next_char in ('<', '=', '(', '+', '\n'):
            self.species_tokens.update(k + next_char for k in self.species_dict)
        self.other_tokens = {'M': 'third-body', 'm': 'third-body',
                             '(+M)': 'falloff3b', '(+m)': 'falloff3b',
                             '<=>': 'equal', '=>': 'equal', '=': 'equal',
                             'HV': 'photon', 'hv': 'photon'}
        self.other_tokens.update(('(+{})'.format(k), 'falloff3b: {}'.format(k))
                                 for k in self.species_dict)
        self.Slen = max(map(len, self.other_tokens))

    def read_kinetics_entry(self, entry, surface):
        """
        Read a kinetics ``entry`` for a single reaction as loaded from a
        Chemkin-format file. Returns a :class:`Reaction` object with the
        reaction and its associated kinetics.
        """

        # Handle non-default units which apply to this entry
        energy_units = self.energy_units
        quantity_units = self.quantity_units
        if 'units' in entry.lower():
            for units in sorted(QUANTITY_UNITS, key=lambda k: -len(k)):
                pattern = re.compile(r'units *\/ *{} *\/'.format(re.escape(units)),
                                     flags=re.IGNORECASE)
                m = pattern.search(entry)
                if m:
                    entry = pattern.sub('', entry)
                    quantity_units = QUANTITY_UNITS[units]
                    break

            for units in sorted(ENERGY_UNITS, key=lambda k: -len(k)):
                pattern = re.compile(r'units *\/ *{} *\/'.format(re.escape(units)),
                                     re.IGNORECASE)
                m = pattern.search(entry)
                if m:
                    entry = pattern.sub('', entry)
                    energy_units = ENERGY_UNITS[units]
                    break

        lines = entry.strip().splitlines()

        # The first line contains the reaction equation and a set of modified Arrhenius parameters
        tokens = lines[0].split()
        A = float(tokens[-3])
        b = float(tokens[-2])
        Ea = float(tokens[-1])
        reaction = ''.join(tokens[:-3]) + '\n'
        original_reaction = reaction  # for use in error messages

        # Identify tokens in the reaction expression in order of
        # decreasing length
        locs = {}
        for i in range(self.Slen, 0, -1):
            for j in range(len(reaction)-i+1):
                test = reaction[j:j+i]
                if test in self.species_tokens:
                    reaction = reaction[:j] + ' '*(i-1) + reaction[j+i-1:]
                    locs[j] = test[:-1], 'species'
                elif test in self.other_tokens:
                    reaction = reaction[:j] + '\n'*i + reaction[j+i:]
                    locs[j] = test, self.other_tokens[test]

        # Anything that's left should be a stoichiometric coefficient or a '+'
        # between species
        for token in reaction.split():
            j = reaction.find(token)
            i = len(token)
            reaction = reaction[:j] + ' '*i + reaction[j+i:]
            if token == '+':
                continue

            try:
                locs[j] = int(token), 'coeff'
            except ValueError:
                try:
                    locs[j] = float(token), 'coeff'
                except ValueError:
                    raise InputError('Unexpected token "{}" in reaction expression "{}".',
                                     token, original_reaction)

        reactants = []
        products = []
        stoichiometry = 1
        lhs = True
        for token, kind in [v for k,v in sorted(locs.items())]:
            if kind == 'equal':
                reversible = token in ('<=>', '=')
                lhs = False
            elif kind == 'coeff':
                stoichiometry = token
            elif lhs:
                reactants.append((stoichiometry, token, kind))
                stoichiometry = 1
            else:
                products.append((stoichiometry, token, kind))
                stoichiometry = 1

        if lhs:
            raise InputError("Failed to find reactant/product delimiter in reaction string.")

        # Create a new Reaction object for this reaction
        reaction = Reaction(reactants=[], products=[], reversible=reversible,
                            parser=self)

        def parse_expression(expression, dest):
            third_body_name = None
            third_body = False  # simple third body reaction (non-falloff)
            photon = False
            for stoichiometry, species, kind in expression:
                if kind == 'third-body':
                    third_body = True
                    third_body_name = 'M'
                elif kind == 'falloff3b':
                    third_body_name = 'M'
                elif kind.startswith('falloff3b:'):
                    third_body_name = kind.split()[1]
                elif kind == 'photon':
                    photon = True
                else:
                    dest.append((stoichiometry, self.species_dict[species]))

            return third_body_name, third_body, photon

        third_body_name_r, third_body, photon_r = parse_expression(reactants, reaction.reactants)
        third_body_name_p, third_body, photon_p = parse_expression(products, reaction.products)

        if third_body_name_r != third_body_name_p:
            raise InputError('Third bodies do not match: "{}" and "{}" in'
                ' reaction entry:\n\n{}', third_body_name_r, third_body_name_p, entry)

        if photon_r:
            raise InputError('Reactant photon not supported. '
                             'Found in reaction:\n{}', entry.strip())
        if photon_p and reversible:
            self.warn('Found reversible reaction containing a product photon:'
                '\n{0}\nIf the "--permissive" option was specified, this will '
                'be converted to an irreversible reaction with the photon '
                'removed.'.format(entry.strip()))
            reaction.reversible = False

        reaction.third_body = third_body_name_r

        # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
        # This assumes elementary kinetics for all reactions
        rStoich = sum(r[0] for r in reaction.reactants) + (1 if third_body else 0)
        if rStoich < 1:
            raise InputError('No reactant species for reaction {}.', reaction)

        length_dim = 3 * (rStoich - 1)
        quantity_dim = rStoich - 1
        kunits = self.get_rate_constant_units(length_dim, 'cm',
                                              quantity_dim, quantity_units)
        klow_units = self.get_rate_constant_units(length_dim + 3, 'cm',
                                                  quantity_dim + 1, quantity_units)

        # The rest of the first line contains Arrhenius parameters
        arrhenius = Arrhenius(
            A=(A, kunits),
            b=b,
            Ea=(Ea, energy_units),
            parser=self
        )

        low_rate = None
        high_rate = None
        falloff = None
        pdep_arrhenius = []
        efficiencies = {}
        coverages = []
        cheb_coeffs = []
        revReaction = None
        is_sticking = None
        motz_wise = None
        Tmin = Tmax = Pmin = Pmax = None  # Chebyshev parameters
        degreeT = degreeP = None

        # Note that the subsequent lines could be in any order
        for line in lines[1:]:
            if not line.strip():
                continue
            tokens = line.split('/')
            parsed = False

            if 'stick' in line.lower():
                parsed = True
                is_sticking = True

            if 'mwon' in line.lower():
                parsed = True
                motz_wise = True

            if 'mwoff' in line.lower():
                parsed = True
                motz_wise = False

            if 'dup' in line.lower():
                # Duplicate reaction
                parsed = True
                reaction.duplicate = True

            if 'low' in line.lower():
                # Low-pressure-limit Arrhenius parameters for "falloff" reaction
                parsed = True
                tokens = tokens[1].split()
                low_rate = Arrhenius(
                    A=(float(tokens[0].strip()), klow_units),
                    b=float(tokens[1].strip()),
                    Ea=(float(tokens[2].strip()), energy_units),
                    parser=self
                )

            elif 'high' in line.lower():
                # High-pressure-limit Arrhenius parameters for "chemically
                # activated" reaction
                parsed = True
                tokens = tokens[1].split()
                high_rate = Arrhenius(
                    A=(float(tokens[0].strip()), kunits),
                    b=float(tokens[1].strip()),
                    Ea=(float(tokens[2].strip()), energy_units),
                    parser=self
                )
                # Need to fix units on the base reaction:
                arrhenius.A = (arrhenius.A[0], klow_units)

            elif 'rev' in line.lower():
                parsed = True
                reaction.reversible = False
                tokens = tokens[1].split()
                # If the A factor in the rev line is zero, don't create the reverse reaction
                if float(tokens[0].strip()) != 0.0:
                    # Create a reaction proceeding in the opposite direction
                    revReaction = Reaction(reactants=reaction.products,
                                           products=reaction.reactants,
                                           third_body=reaction.third_body,
                                           reversible=False,
                                           parser=self)

                    rev_rate = Arrhenius(
                        A=(float(tokens[0].strip()), klow_units),
                        b=float(tokens[1].strip()),
                        Ea=(float(tokens[2].strip()), energy_units),
                        parser=self
                    )
                    if third_body:
                        revReaction.kinetics = ThreeBody(rev_rate)
                    else:
                        revReaction.kinetics = ElementaryRate(rev_rate)

            elif 'ford' in line.lower():
                parsed = True
                tokens = tokens[1].split()
                reaction.forward_orders[tokens[0].strip()] = float(tokens[1])

            elif 'troe' in line.lower():
                # Troe falloff parameters
                parsed = True
                tokens = tokens[1].split()
                falloff = Troe(A=float(tokens[0].strip()),
                               T3=float(tokens[1].strip()),
                               T1=float(tokens[2].strip()),
                               T2=float(tokens[3].strip()) if len(tokens) > 3 else None)
            elif 'sri' in line.lower():
                # SRI falloff parameters
                parsed = True
                tokens = tokens[1].split()
                A = float(tokens[0].strip())
                B = float(tokens[1].strip())
                C = float(tokens[2].strip())
                try:
                    D = float(tokens[3].strip())
                    E = float(tokens[4].strip())
                except (IndexError, ValueError):
                    D = None
                    E = None

                if D is None or E is None:
                    falloff = Sri(A=A, B=B, C=C)
                else:
                    falloff = Sri(A=A, B=B, C=C, D=D, E=E)

            elif 'cov' in line.lower():
                parsed = True
                C = tokens[1].split()
                coverages.append(
                    [C[0], fortFloat(C[1]), fortFloat(C[2]), fortFloat(C[3])])

            elif 'cheb' in line.lower():
                # Chebyshev parameters
                parsed = True
                tokens = [t.strip() for t in tokens]
                if contains(tokens, 'TCHEB'):
                    index = get_index(tokens, 'TCHEB')
                    tokens2 = tokens[index+1].split()
                    Tmin = float(tokens2[0].strip())
                    Tmax = float(tokens2[1].strip())
                if contains(tokens, 'PCHEB'):
                    index = get_index(tokens, 'PCHEB')
                    tokens2 = tokens[index+1].split()
                    Pmin = (float(tokens2[0].strip()), 'atm')
                    Pmax = (float(tokens2[1].strip()), 'atm')
                if contains(tokens, 'TCHEB') or contains(tokens, 'PCHEB'):
                    pass
                elif degreeT is None or degreeP is None:
                    tokens2 = tokens[1].split()
                    degreeT = int(float(tokens2[0].strip()))
                    degreeP = int(float(tokens2[1].strip()))
                    cheb_coeffs.extend([float(t.strip()) for t in tokens2[2:]])
                else:
                    tokens2 = tokens[1].split()
                    cheb_coeffs.extend([float(t.strip()) for t in tokens2])

            elif 'plog' in line.lower():
                # Pressure-dependent Arrhenius parameters
                parsed = True
                tokens = tokens[1].split()
                pdep_arrhenius.append([float(tokens[0].strip()), Arrhenius(
                    A=(float(tokens[1].strip()), kunits),
                    b=float(tokens[2].strip()),
                    Ea=(float(tokens[3].strip()), energy_units),
                    parser=self
                )])
            elif len(tokens) >= 2:
                # Assume a list of collider efficiencies
                parsed = True
                for collider, efficiency in zip(tokens[0::2], tokens[1::2]):
                    efficiencies[collider.strip()] = float(efficiency.strip())

            if not parsed:
                raise InputError('Unparsable line:\n"""\n{}\n"""', line)

        # Decide which kinetics to keep and store them on the reaction object.
        # At most one of the special cases should be true
        tests = [cheb_coeffs, pdep_arrhenius, low_rate, high_rate, third_body,
                 surface]
        if sum(bool(t) for t in tests) > 1:
            raise InputError('Reaction {} contains parameters for more than '
                'one reaction type.', original_reaction)

        if cheb_coeffs:
            if Tmin is None or Tmax is None:
                raise InputError('Missing TCHEB line for reaction {}', reaction)
            if Pmin is None or Pmax is None:
                raise InputError('Missing PCHEB line for reaction {}', reaction)
            if len(cheb_coeffs) != degreeT * degreeP:
                raise InputError('Incorrect number of Chebyshev coefficients. '
                    'Expected {}*{} = {} but got {}', degreeT, degreeP,
                    degreeT * degreeP, len(cheb_coeffs))
            if quantity_units == self.quantity_units:
                quantity_units = None
            reaction.kinetics = Chebyshev(
                Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax,
                quantity_units=quantity_units,
                coeffs=np.array(cheb_coeffs, np.float64).reshape((degreeT, degreeP)))
        elif pdep_arrhenius:
            reaction.kinetics = PDepArrhenius(
                pressures=[P for P, arrh in pdep_arrhenius],
                pressure_units="atm",
                arrhenius=[arrh for P, arrh in pdep_arrhenius]
            )
        elif low_rate is not None:
            reaction.kinetics = Falloff(high_rate=arrhenius,
                                        low_rate=low_rate,
                                        F=falloff,
                                        efficiencies=efficiencies)
        elif high_rate is not None:
            reaction.kinetics = ChemicallyActivated(high_rate=high_rate,
                                                    low_rate=arrhenius,
                                                    F=falloff,
                                                    efficiencies=efficiencies)
        elif third_body:
            reaction.kinetics = ThreeBody(high_rate=arrhenius,
                                          efficiencies=efficiencies)
        elif reaction.third_body:
            raise InputError('Reaction equation implies pressure '
                'dependence but no alternate rate parameters (i.e. HIGH or '
                'LOW) were given for reaction {}', reaction)
        elif surface:
            reaction.kinetics = SurfaceRate(rate=arrhenius,
                                            coverages=coverages,
                                            is_sticking=is_sticking,
                                            motz_wise=motz_wise)
        else:
            reaction.kinetics = ElementaryRate(arrhenius)

        if revReaction:
            revReaction.duplicate = reaction.duplicate
            revReaction.kinetics.efficiencies = reaction.kinetics.efficiencies

        return reaction, revReaction

    def load_extra_file(self, path):
        """
        Load YAML-formatted entries from ``path`` on disk.
        """
        with open(path, 'rt', encoding="utf-8") as stream:
            yml = yaml.round_trip_load(stream)

        # do not overwrite reserved field names
        reserved = {'generator', 'input-files', 'cantera-version', 'date',
                    'units', 'phases', 'species', 'reactions'}
        reserved &= set(yml.keys())
        if reserved:
            raise InputError("The YAML file '{}' provided as '--extra' input "
                "must not redefine reserved field name: "
                "'{}'".format(path, reserved))

        # replace header lines
        if 'description' in yml:
            if isinstance(yml['description'], str):
                if self.header_lines:
                    self.header_lines += ['']
                self.header_lines += yml.pop('description').split('\n')
            else:
                raise InputError("The alternate description provided in "
                    "'{}' needs to be a string".format(path))

        # remainder
        self.extra = yml

    def load_chemkin_file(self, path, skip_undeclared_species=True, surface=False):
        """
        Load a Chemkin-format input file from ``path`` on disk.
        """
        transportLines = []
        self.line_number = 0

        with open(path, 'r', errors='ignore') as ck_file:

            def readline():
                self.line_number += 1
                line = strip_nonascii(ck_file.readline())
                if '!' in line:
                    return line.split('!', 1)
                elif line:
                    return line, ''
                else:
                    return None, None

            # @TODO: This loop is a bit of a mess, and could probably be cleaned
            # up by refactoring it into a set of methods for processing each
            # input file section.
            line, comment = readline()
            advance = True
            inHeader = True
            header = []
            indent = 80
            while line is not None:
                tokens = line.split() or ['']
                if inHeader and not line.strip():
                    header.append(comment.rstrip())
                    if comment.strip() != '': # skip indent calculation if empty
                        indent = min(indent, re.search('[^ ]', comment).start())

                if tokens[0].upper().startswith('ELEM'):
                    inHeader = False
                    tokens = tokens[1:]
                    while line is not None and get_index(line, 'END') is None:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('SPEC', 'SPECIES'):
                            self.warn('"ELEMENTS" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            tokens.pop()
                            break

                        line, comment = readline()
                        # Normalize custom atomic weights
                        line = re.sub(r'\s*/\s*([0-9\.EeDd+-]+)\s*/', r'/\1/ ', line)
                        tokens.extend(line.split())

                    for token in tokens:
                        if token.upper() == 'END':
                            break
                        self.add_element(token)

                elif tokens[0].upper().startswith('SPEC'):
                    # List of species identifiers
                    species = tokens[1:]
                    inHeader = False
                    comments = {}
                    while line is not None and get_index(line, 'END') is None:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('REAC', 'REACTIONS', 'TRAN',
                                                  'TRANSPORT', 'THER', 'THERMO'):
                            self.warn('"SPECIES" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            species.pop()
                            # Fix the case where there THERMO ALL or REAC UNITS
                            # ends the species section
                            if (species[-1].upper().startswith('THER') or
                                    species[-1].upper().startswith('REAC')):
                                species.pop()
                            break

                        line, comment = readline()
                        comment = comment.strip()
                        line_species = line.split()
                        if len(line_species) == 1 and comment:
                            comments[line_species[0]] = comment
                        species.extend(line_species)

                    for token in species:
                        if token.upper() == 'END':
                            break
                        if token in self.species_dict:
                            species = self.species_dict[token]
                            self.warn('Found additional declaration of species {}'.format(species))
                        else:
                            species = Species(label=token)
                            if token in comments:
                                species.note = comments[token]
                            self.species_dict[token] = species
                            self.species_list.append(species)

                elif tokens[0].upper().startswith('SITE'):
                    # List of species identifiers for surface species
                    if '/' in tokens[0]:
                        surf_name = tokens[0].split('/')[1]
                    else:
                        surf_name = 'surface{}'.format(len(self.surfaces)+1)
                    tokens = tokens[1:]
                    site_density = None
                    for token in tokens[:]:
                        if token.upper().startswith('SDEN/'):
                            site_density = fortFloat(token.split('/')[1])
                            tokens.remove(token)

                    if site_density is None:
                        raise InputError('SITE section defined with no site density')
                    self.surfaces.append(Surface(name=surf_name,
                                                 site_density=site_density))
                    surf = self.surfaces[-1]

                    inHeader = False
                    while line is not None and get_index(line, 'END') is None:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('REAC', 'REACTIONS', 'THER',
                                                  'THERMO'):
                            self.warn('"SITE" section implicitly ended by start of '
                                      'next section on line {}.'.format(self.line_number))
                            advance = False
                            tokens.pop()
                            # Fix the case where there THERMO ALL or REAC UNITS
                            # ends the species section
                            if (tokens[-1].upper().startswith('THER') or
                                    tokens[-1].upper().startswith('REAC')):
                                tokens.pop()
                            break

                        line, comment = readline()
                        tokens.extend(line.split())

                    for token in tokens:
                        if token.upper() == 'END':
                            break
                        if token.count('/') == 2:
                            # species occupies a specific number of sites
                            token, sites, _ = token.split('/')
                            sites = float(sites)
                        else:
                            sites = None
                        if token in self.species_dict:
                            species = self.species_dict[token]
                            self.warn('Found additional declaration of species {0}'.format(species))
                        else:
                            species = Species(label=token, sites=sites)
                            self.species_dict[token] = species
                            surf.species_list.append(species)

                elif tokens[0].upper().startswith('THER') and contains(line, 'NASA9'):
                    inHeader = False
                    entryLength = None
                    entry = []
                    # Gather comments on lines preceding and within this entry
                    comments = []
                    while line is not None and get_index(line, 'END') != 0:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('REAC', 'REACTIONS', 'TRAN', 'TRANSPORT'):
                            self.warn('"THERMO" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            tokens.pop()
                            break

                        line, comment = readline()
                        comments.append(comment)
                        if not line:
                            continue

                        if entryLength is None:
                            entryLength = 0
                            # special case if (redundant) temperature ranges are
                            # given as the first line
                            try:
                                s = line.split()
                                float(s[0]), float(s[1]), float(s[2])
                                continue
                            except (IndexError, ValueError):
                                pass

                        entry.append(line)
                        if len(entry) == 2:
                            entryLength = 2 + 3 * int(line.split()[0])

                        if len(entry) == entryLength:
                            label, thermo, comp = self.read_NASA9_entry(entry, comments)
                            comments = []
                            entry = []
                            if label not in self.species_dict:
                                if skip_undeclared_species:
                                    logger.info('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
                                    continue
                                else:
                                    # Add a new species entry
                                    species = Species(label=label)
                                    self.species_dict[label] = species
                                    self.species_list.append(species)
                            else:
                                species = self.species_dict[label]

                            # use the first set of thermo data found
                            if species.thermo is not None:
                                self.warn('Found additional thermo entry for species {0}. '
                                          'If --permissive was given, the first entry is used.'.format(label))
                            else:
                                species.thermo = thermo
                                species.composition = comp

                elif tokens[0].upper().startswith('THER'):
                    # List of thermodynamics (hopefully one per species!)
                    inHeader = False
                    line, comment = readline()
                    if line is not None and get_index(line, 'END') is None:
                        TintDefault = float(line.split()[1])
                    thermo = []
                    current = []
                    # Gather comments on lines preceding and within this entry
                    comments = [comment]
                    while line is not None and get_index(line, 'END') != 0:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('REAC', 'REACTIONS', 'TRAN', 'TRANSPORT'):
                            self.warn('"THERMO" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            tokens.pop()
                            break

                        if comment:
                            current.append('!'.join((line, comment)))
                        else:
                            current.append(line)
                        if len(line) >= 80 and line[79] in ['1', '2', '3', '4']:
                            thermo.append(line)
                            if line[79] == '4':
                                try:
                                    label, thermo, comp = self.read_NASA7_entry(thermo, TintDefault, comments)
                                except Exception as e:
                                    error_line_number = self.line_number - len(current) + 1
                                    error_entry = ''.join(current).rstrip()
                                    logger.info(
                                        'Error while reading thermo entry starting on line {0}:\n'
                                        '"""\n{1}\n"""'.format(error_line_number, error_entry)
                                    )
                                    raise

                                if label not in self.species_dict:
                                    if skip_undeclared_species:
                                        logger.info(
                                            'Skipping unexpected species "{0}" while'
                                            ' reading thermodynamics entry.'.format(label))
                                        thermo = []
                                        line, comment = readline()
                                        current = []
                                        comments = [comment]
                                        continue
                                    else:
                                        # Add a new species entry
                                        species = Species(label=label)
                                        self.species_dict[label] = species
                                        self.species_list.append(species)
                                else:
                                    species = self.species_dict[label]

                                # use the first set of thermo data found
                                if species.thermo is not None:
                                    self.warn('Found additional thermo entry for species {0}. '
                                              'If --permissive was given, the first entry is used.'.format(label))
                                else:
                                    species.thermo = thermo
                                    species.composition = comp

                                thermo = []
                                current = []
                                comments = []
                        elif thermo and thermo[-1].rstrip().endswith('&'):
                            # Include Chemkin-style extended elemental composition
                            thermo.append(line)
                        line, comment = readline()
                        comments.append(comment)

                elif tokens[0].upper().startswith('REAC'):
                    # Reactions section
                    inHeader = False
                    for token in tokens[1:]:
                        token = token.upper()
                        if token in ENERGY_UNITS:
                            self.energy_units = ENERGY_UNITS[token]
                            if not self.processed_units:
                                self.output_energy_units = ENERGY_UNITS[token]
                        elif token in QUANTITY_UNITS:
                            self.quantity_units = QUANTITY_UNITS[token]
                            if not self.processed_units:
                                self.output_quantity_units = QUANTITY_UNITS[token]
                        elif token == 'MWON':
                            self.motz_wise = True
                        elif token == 'MWOFF':
                            self.motz_wise = False
                        else:
                            raise InputError("Unrecognized token on REACTIONS line, {0!r}", token)

                    self.processed_units = True

                    kineticsList = []
                    commentsList = []
                    startLines = []
                    kinetics = ''
                    comments = ''

                    line, comment = readline()
                    if surface:
                        reactions = self.surfaces[-1].reactions
                    else:
                        reactions = self.reactions
                    while line is not None and get_index(line, 'END') is None:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('TRAN', 'TRANSPORT'):
                            self.warn('"REACTIONS" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            break

                        lineStartsWithComment = not line and comment
                        line = line.rstrip()
                        comment = comment.rstrip()

                        if '=' in line and not lineStartsWithComment:
                            # Finish previous record
                            if comment:
                                # End of line comment belongs with this reaction
                                comments += comment + '\n'
                                comment = ''
                            kineticsList.append(kinetics)
                            commentsList.append(comments)
                            startLines.append(self.line_number)
                            kinetics = ''
                            comments = ''

                        if line.strip():
                            kinetics += line + '\n'
                        if comment:
                            comments += comment + '\n'

                        line, comment = readline()

                    # Don't forget the last reaction!
                    if kinetics.strip() != '':
                        kineticsList.append(kinetics)
                        commentsList.append(comments)

                    # We don't actually know whether comments belong to the
                    # previous or next reaction, but to keep them positioned
                    # correctly, we associate them with the next reaction. A
                    # comment after the last reaction is associated with that
                    # reaction
                    if kineticsList and kineticsList[0] == '':
                        kineticsList.pop(0)
                        final_comment = commentsList.pop()
                        if final_comment and commentsList[-1]:
                            commentsList[-1] = commentsList[-1].rstrip() + '\n' + final_comment
                        elif final_comment:
                            commentsList[-1] = final_comment

                    self.setup_kinetics()
                    for kinetics, comment, line_number in zip(kineticsList, commentsList, startLines):
                        try:
                            reaction, revReaction = self.read_kinetics_entry(kinetics, surface)
                        except Exception as e:
                            self.line_number = line_number
                            logger.info('Error reading reaction starting on '
                                'line {0}:\n"""\n{1}\n"""'.format(
                                    line_number, kinetics.rstrip()))
                            raise
                        reaction.line_number = line_number
                        reaction.comment = comment
                        reactions.append(reaction)
                        if revReaction is not None:
                            revReaction.line_number = line_number
                            reactions.append(revReaction)

                elif tokens[0].upper().startswith('TRAN'):
                    inHeader = False
                    line, comment = readline()
                    transport_start_line = self.line_number
                    while line is not None and get_index(line, 'END') is None:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('REAC', 'REACTIONS'):
                            self.warn('"TRANSPORT" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            tokens.pop()
                            break

                        if comment:
                            transportLines.append('!'.join((line, comment)))
                        else:
                            transportLines.append(line)
                        line, comment = readline()

                elif line.strip():
                    raise InputError('Section starts with unrecognized keyword'
                        '\n"""\n{}\n"""', line.rstrip())

                if advance:
                    line, comment = readline()
                else:
                    advance = True

        for h in header:
            self.header_lines.append(h[indent:])

        for index, reaction in enumerate(self.reactions):
            reaction.index = index + 1

        if transportLines:
            self.parse_transport_data(transportLines, path, transport_start_line)

    def parse_transport_data(self, lines, filename, line_offset):
        """
        Parse the Chemkin-format transport data in ``lines`` (a list of strings)
        and add that transport data to the previously-loaded species.
        """

        for i,line in enumerate(lines):
            original_line = line
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            if get_index(line, 'END') == 0:
                break

            if '!' in line:
                line, comment = line.split('!', 1)
            else:
                comment = ''

            data = line.split()

            speciesName = data[0]
            if speciesName in self.species_dict:
                if len(data) != 7:
                    raise InputError('Unable to parse line {} of {}:\n"""\n{}"""\n'
                        '6 transport parameters expected, but found {}.',
                            line_offset + i, filename, original_line, len(data)-1)

                if self.species_dict[speciesName].transport is None:
                    self.species_dict[speciesName].transport = TransportData(*data, note=comment)
                else:
                    self.warn('Ignoring duplicate transport data'
                         ' for species "{}" on line {} of "{}".'.format(
                            speciesName, line_offset + i, filename))


    def write_yaml(self, name='gas', out_name='mech.yaml'):
        emitter = yaml.YAML()
        emitter.width = 70

        emitter.register_class(Species)
        emitter.register_class(Nasa7)
        emitter.register_class(Nasa9)
        emitter.register_class(TransportData)
        emitter.register_class(Reaction)

        with open(out_name, 'w') as dest:
            have_transport = True
            for s in self.species_list:
                if not s.transport:
                    have_transport = False

            surface_names = []
            n_reacting_phases = 0
            if self.reactions:
                n_reacting_phases += 1
            for surf in self.surfaces:
                surface_names.append(surf.name)
                if surf.reactions:
                    n_reacting_phases += 1

            # Write header lines
            desc = '\n'.join(line.rstrip() for line in self.header_lines)
            desc = desc.strip('\n')
            desc = textwrap.dedent(desc)
            if desc.strip():
                emitter.dump({'description': yaml.scalarstring.PreservedScalarString(desc)}, dest)

            # Additional information regarding conversion
            files = [os.path.basename(f) for f in self.files]
            metadata = BlockMap([
                ('generator', 'ck2yaml'),
                ('input-files', FlowList(files)),
                ('cantera-version', '2.5.0a4'),
                ('date', formatdate(localtime=True)),
            ])
            if desc.strip():
                metadata.yaml_set_comment_before_after_key('generator', before='\n')
            emitter.dump(metadata, dest)

            # Write extra entries
            if self.extra:
                extra = BlockMap(self.extra)
                key = list(self.extra.keys())[0]
                extra.yaml_set_comment_before_after_key(key, before='\n')
                emitter.dump(extra, dest)

            units = FlowMap([('length', 'cm'), ('time', 's')])
            units['quantity'] = self.output_quantity_units
            units['activation-energy'] = self.output_energy_units
            units_map = BlockMap([('units', units)])
            units_map.yaml_set_comment_before_after_key('units', before='\n')
            emitter.dump(units_map, dest)

            phases = []
            reactions = []
            if name is not None:
                phase = BlockMap()
                phase['name'] = name
                phase['thermo'] = 'ideal-gas'
                phase['elements'] = FlowList(self.elements)
                phase['species'] = FlowList(S.label for S in self.species_list)
                if self.reactions:
                    phase['kinetics'] = 'gas'
                    if n_reacting_phases == 1:
                        reactions.append(('reactions', self.reactions))
                    else:
                        rname = '{}-reactions'.format(name)
                        phase['reactions'] = [rname]
                        reactions.append((rname, self.reactions))
                if have_transport:
                    phase['transport'] = 'mixture-averaged'
                phase['state'] = FlowMap([('T', 300.0), ('P', '1 atm')])
                phases.append(phase)

            for surf in self.surfaces:
                # Write definitions for surface phases
                phase = BlockMap()
                phase['name'] = surf.name
                phase['thermo'] = 'ideal-surface'
                phase['elements'] = FlowList(self.elements)
                phase['species'] = FlowList(S.label for S in surf.species_list)
                phase['site-density'] = surf.site_density
                if self.motz_wise is not None:
                    phase['Motz-Wise'] = self.motz_wise
                if surf.reactions:
                    phase['kinetics'] = 'surface'
                    if n_reacting_phases == 1:
                        reactions.append(('reactions', surf.reactions))
                    else:
                        rname = '{}-reactions'.format(surf.name)
                        phase['reactions'] = [rname]
                        reactions.append((rname, surf.reactions))
                phase['state'] = FlowMap([('T', 300.0), ('P', '1 atm')])
                phases.append(phase)

            if phases:
                phases_map = BlockMap([('phases', phases)])
                phases_map.yaml_set_comment_before_after_key('phases', before='\n')
                emitter.dump(phases_map, dest)

            # Write data on custom elements
            if self.element_weights:
                elements = []
                for name, weight in sorted(self.element_weights.items()):
                    E = BlockMap([('symbol', name), ('atomic-weight', weight)])
                    elements.append(E)
                elementsMap = BlockMap([('elements', elements)])
                elementsMap.yaml_set_comment_before_after_key('elements', before='\n')
                emitter.dump(elementsMap, dest)

            # Write the individual species data
            all_species = list(self.species_list)
            for species in all_species:
                if species.composition is None:
                    raise InputError('No thermo data found for '
                                     'species {!r}'.format(species.label))

            for surf in self.surfaces:
                all_species.extend(surf.species_list)
            speciesMap = BlockMap([('species', all_species)])
            speciesMap.yaml_set_comment_before_after_key('species', before='\n')
            emitter.dump(speciesMap, dest)

            # Write the reactions section(s)
            for label, R in reactions:
                reactionsMap = BlockMap([(label, R)])
                reactionsMap.yaml_set_comment_before_after_key(label, before='\n')
                emitter.dump(reactionsMap, dest)

        # Names of surface phases need to be returned so they can be imported as
        # part of mechanism validation
        return surface_names

    @staticmethod
    def convert_mech(input_file, thermo_file=None, transport_file=None,
                     surface_file=None, phase_name='gas', extra_file=None,
                     out_name=None, quiet=False, permissive=None):

        parser = Parser()
        if quiet:
            logger.setLevel(level=logging.ERROR)
        else:
            logger.setLevel(level=logging.INFO)

        if permissive is not None:
            parser.warning_as_error = not permissive

        if input_file:
            parser.files.append(input_file)
            input_file = os.path.expanduser(input_file)
            if not os.path.exists(input_file):
                raise IOError('Missing input file: {0!r}'.format(input_file))
            try:
                # Read input mechanism files
                parser.load_chemkin_file(input_file)
            except Exception as err:
                logger.warning("\nERROR: Unable to parse '{0}' near line {1}:\n{2}\n".format(
                                input_file, parser.line_number, err))
                raise
        else:
            phase_name = None

        if thermo_file:
            parser.files.append(thermo_file)
            thermo_file = os.path.expanduser(thermo_file)
            if not os.path.exists(thermo_file):
                raise IOError('Missing thermo file: {0!r}'.format(thermo_file))
            try:
                parser.load_chemkin_file(thermo_file,
                                       skip_undeclared_species=bool(input_file))
            except Exception:
                logger.warning("\nERROR: Unable to parse '{0}' near line {1}:\n".format(
                               thermo_file, parser.line_number))
                raise

        if transport_file:
            parser.files.append(transport_file)
            transport_file = os.path.expanduser(transport_file)
            if not os.path.exists(transport_file):
                raise IOError('Missing transport file: {0!r}'.format(transport_file))
            with open(transport_file, 'r', errors='ignore') as f:
                lines = [strip_nonascii(line) for line in f]
            parser.parse_transport_data(lines, transport_file, 1)

            # Transport validation: make sure all species have transport data
            for s in parser.species_list:
                if s.transport is None:
                    raise InputError("No transport data for species '{}'.", s)

        if surface_file:
            parser.files.append(surface_file)
            surface_file = os.path.expanduser(surface_file)
            if not os.path.exists(surface_file):
                raise IOError('Missing input file: {0!r}'.format(surface_file))
            try:
                # Read input mechanism files
                parser.load_chemkin_file(surface_file, surface=True)
            except Exception as err:
                logger.warning("\nERROR: Unable to parse '{0}' near line {1}:\n{2}\n".format(
                               surface_file, parser.line_number, err))
                raise

        if extra_file:
            parser.files.append(extra_file)
            extra_file = os.path.expanduser(extra_file)
            if not os.path.exists(extra_file):
                raise IOError('Missing input file: {0!r}'.format(extra_file))
            try:
                # Read input mechanism files
                parser.load_extra_file(extra_file)
            except Exception as err:
                logger.warning("\nERROR: Unable to parse '{0}':\n{1}\n".format(
                               extra_file, err))
                raise

        if out_name:
            out_name = os.path.expanduser(out_name)
        else:
            out_name = os.path.splitext(input_file)[0] + '.yaml'

        # Write output file
        surface_names = parser.write_yaml(name=phase_name, out_name=out_name)
        if not quiet:
            nReactions = len(parser.reactions) + sum(len(surf.reactions) for surf in parser.surfaces)
            logger.info('Wrote YAML mechanism file to {0!r}.'.format(out_name))
            logger.info('Mechanism contains {0} species and {1} reactions.'.format(
                        len(parser.species_list), nReactions))
        return parser, surface_names

    def show_duplicate_reactions(self, error_message):
        # Find the reaction numbers of the duplicate reactions by looking at
        # the YAML file lines shown in the error message generated by
        # Kinetics::checkDuplicates.
        reactions = []
        for line in error_message.split('\n'):
            match = re.match('>.*# Reaction ([0-9]+)', line)
            if match:
                reactions.append(int(match.group(1))-1)

        if len(reactions) != 2:
            # Something went wrong while parsing the error message, so just
            # display it as-is instead of trying to be clever.
            logger.warning(error_message)
            return

        # Give an error message that references the line numbers in the
        # original input file.
        equation = str(self.reactions[reactions[0]])
        lines = [self.reactions[i].line_number for i in reactions]
        logger.warning('Undeclared duplicate reaction {}\nfound on lines {} and {} of '
              'the kinetics input file.'.format(equation, lines[0], lines[1]))


def convert_mech(input_file, thermo_file=None, transport_file=None,
                 surface_file=None, phase_name='gas', extra_file=None,
                 out_name=None, quiet=False, permissive=None):
    _, surface_names = Parser.convert_mech(
        input_file, thermo_file, transport_file, surface_file, phase_name,
        extra_file, out_name, quiet, permissive)
    return surface_names


def main(argv):

    longOptions = ['input=', 'thermo=', 'transport=', 'surface=', 'name=',
                   'extra=', 'output=', 'permissive', 'help', 'debug', 'quiet',
                   'no-validate', 'id=']

    try:
        optlist, args = getopt.getopt(argv, 'dh', longOptions)
        options = dict()
        for o,a in optlist:
            options[o] = a

        if args:
            raise getopt.GetoptError('Unexpected command line option: ' +
                                     repr(' '.join(args)))

    except getopt.GetoptError as e:
        logger.error('ck2yaml.py: Error parsing arguments:')
        logger.error(e)
        logger.error('Run "ck2yaml.py --help" to see usage help.')
        sys.exit(1)

    if not options or '-h' in options or '--help' in options:
        logger.info(__doc__)
        sys.exit(0)

    input_file = options.get('--input')
    thermo_file = options.get('--thermo')
    permissive = '--permissive' in options
    quiet = '--quiet' in options
    transport_file = options.get('--transport')
    surface_file = options.get('--surface')

    if '--id' in options:
        phase_name = options.get('--id', 'gas')
        logger.warning("\nFutureWarning: "
                       "Option '--id=...' will be replaced by '--name=...'")
    else:
        phase_name = options.get('--name', 'gas')

    if not input_file and not thermo_file:
        logger.error('At least one of the arguments "--input=..." or "--thermo=..."'
                     ' must be provided.\nRun "ck2yaml.py --help" to see usage help.')
        sys.exit(1)

    extra_file = options.get('--extra')

    if '--output' in options:
        out_name = options['--output']
        if not out_name.endswith('.yaml') and not out_name.endswith('.yml'):
            out_name += '.yaml'
    elif input_file:
        out_name = os.path.splitext(input_file)[0] + '.yaml'
    else:
        out_name = os.path.splitext(thermo_file)[0] + '.yaml'

    parser, surfaces = Parser.convert_mech(input_file, thermo_file,
            transport_file, surface_file, phase_name, extra_file, out_name,
            quiet, permissive)

    # Do full validation by importing the resulting mechanism
    if not input_file:
        # Can't validate input files that don't define a phase
        return

    if '--no-validate' in options:
        return

    try:
        import cantera as ct
    except ImportError:
        logger.warning('WARNING: Unable to import Cantera Python module. '
                        'Output mechanism has not been validated')
        sys.exit(0)

    try:
        logger.info('Validating mechanism...')
        gas = ct.Solution(out_name)
        for surf_name in surfaces:
            phase = ct.Interface(out_name, surf_name, [gas])
        logger.info('PASSED')
    except RuntimeError as e:
        logger.info('FAILED')
        msg = str(e)
        if 'Undeclared duplicate reactions' in msg:
            parser.show_duplicate_reactions(msg)
        else:
            logger.warning(e)
        sys.exit(1)


def script_entry_point():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
