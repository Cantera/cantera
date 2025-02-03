#!/usr/bin/env python3
# encoding: utf-8

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""Convert Chemkin-format mechanism files to YAML.

There are two main entry points to this script, `main` and `convert`. The former is
used from the command line interface and parses the arguments passed. The latter uses
arguments that correspond to options of the command line interface.
"""

import logging
import os.path
import sys
from typing import Callable
import numpy as np
import re
import warnings
import argparse
import textwrap
from email.utils import formatdate
from pathlib import Path
from ruamel import yaml

# yaml.version_info is a tuple with the three parts of the version
yaml_version = yaml.version_info
# We choose ruamel.yaml 0.17.16 as the minimum version since it is the highest version
# available in the Ubuntu 22.04 repositories.
yaml_min_version = (0, 17, 16)
if yaml_version < yaml_min_version:
    raise RuntimeError(
        "The minimum supported version of ruamel.yaml is 0.17.16. If you "
        "installed ruamel.yaml from your operating system's package manager, "
        "please install an updated version using pip or conda."
    )

BlockMap = yaml.comments.CommentedMap

class ErrorFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord):
        message = super().format(record)
        if record.levelno >= logging.ERROR:
            return '*' * 79 + '\n' + message
        elif '\n' in message and not message.endswith('\n'):
            return message + '\n'
        else:
            return message

logger = logging.getLogger(__name__)

def FlowMap(*args, **kwargs):
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m

def FlowList(*args, **kwargs):
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst

def represent_float(self, data):
    if data != data:
        value = '.nan'
    elif data == self.inf_value:
        value = '.inf'
    elif data == -self.inf_value:
        value = '-.inf'
    elif data == 0:
        value = '0.0'
    elif 0.01 <= abs(data) < 10000:
        value = np.format_float_positional(data, trim='0')
    else:
        value = np.format_float_scientific(data, trim='0')

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
    def __init__(self, message):
        message += ("\nPlease check https://cantera.org/stable/userguide/"
                   "ck2yaml-tutorial.html#debugging-common-errors-in-ck-files"
                   "\nfor the correct Chemkin syntax.")
        super().__init__(message)


class ParserLogger(logging.Handler):
    def __init__(self, parser):
        self.parser = parser
        self.errors = []
        super().__init__()

    def emit(self, record: logging.LogRecord):
        self.parser.max_loglevel = max(self.parser.max_loglevel, record.levelno)
        self.errors.append(self.format(record))


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
    polynomials. See :ref:`sec-thermo-nasa7`.
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
        if node.Tmid is not None:
            out['temperature-ranges'] = FlowList([node.Tmin, node.Tmid, node.Tmax])
            out['data'] = [FlowList(node.low_coeffs), FlowList(node.high_coeffs)]
        else:
            out['temperature-ranges'] = FlowList([node.Tmin, node.Tmax])
            out['data'] = [FlowList(node.low_coeffs)]
        if node.note:
            note = textwrap.dedent(node.note.rstrip())
            if '\n' in note:
                note = yaml.scalarstring.PreservedScalarString(note)
            out['note'] = note
        return representer.represent_dict(out)


class Nasa9:
    """
    Thermodynamic data parameterized as any number of nine-coefficient NASA
    polynomials. See :ref:`sec-thermo-nasa9`.

    :param data:
        List of polynomials, where each polynomial is written as
        ```
        [(T_low, T_high), [a_0, a_1, ..., a_8]]
        ```
    """
    def __init__(self, *, parser, data, note=''):
        self.note = note
        self.data = list(sorted(data))
        self.Tranges = [self.data[0][0][0]]
        for i in range(1, len(data)):
            if abs(self.data[i-1][0][1] - self.data[i][0][0]) > 0.01:
                logger.error(parser.entry('thermo entry') +
                    'NASA9 polynomials contain non-adjacent temperature ranges')
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

    def __init__(self, index=-1, reactants=None, products=None,
                 kinetics=None, reversible=True, duplicate=False,
                 forward_orders=None, third_body=None):
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
        Return a string representation of the reaction, such as 'A + B <=> C + D'.
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
                logger.info("Negative reaction order for reaction "
                            f"{node.index} ({node!s}).")
            reactant_names = {r[1].label for r in node.reactants}
            if any((species not in reactant_names for species in node.forward_orders)):
                out['nonreactant-orders'] = True
                logger.info(f"Non-reactant order for reaction {node.index} ({node!s}).")
        if node.comment:
            comment = textwrap.dedent(node.comment)
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
    See :ref:`sec-arrhenius-rate`.

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
    See :ref:`sec-surface-rate`.

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
    various pressures. See :ref:`sec-plog-rate`.

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
    See :ref:`sec-chebyshev-rate`.

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
    See :ref:`sec-three-body-reaction`.

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
    See :ref:`sec-falloff-rate`.

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
    See :ref:`sec-chemically-activated-rate`.

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
    See :ref:`sec-troe-falloff`.
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
    See :ref:`sec-sri-falloff`.
    """
    def __init__(self, *, A, B, C, D=None, E=None):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def reduce(self, output):
        sri = FlowMap([('A', self.A), ('B', self.B), ('C', self.C)])
        if self.D is not None:
            sri['D'] = self.D
        if self.E is not None:
            sri['E'] = self.E

        output['SRI'] = sri


class TransportData:
    geometry_flags = ['atom', 'linear', 'nonlinear']

    def __init__(self, parser, label, geometry, well_depth, collision_diameter,
                 dipole_moment, polarizability, z_rot, note=''):

        try:
            geometry = int(geometry)
        except ValueError:
            try:
                geometry = float(geometry)
            except ValueError:
                logger.error(parser.entry("transport data") +
                    f"Invalid geometry flag '{geometry}' for species '{label}'. "
                    "Flag must be an integer.")
                return
            if geometry == int(geometry):
                # This is a minor issue, so even at the default level it's just a
                # warning, and in permissive mode we won't even mention it.
                if not parser.permissive:
                    logger.info(f"Incorrect geometry flag syntax '{geometry}' for "
                        f"species '{label}'. The flag was automatically converted to "
                        "an integer.")
                geometry = int(geometry)
            else:
                logger.error(parser.entry("transport data") +
                    f"Invalid float geometry flag '{geometry}' for species '{label}'. "
                    "Flag must be an integer.")
                return
        if geometry not in (0, 1, 2):
            logger.error(parser.entry("transport data") +
                f"Invalid geometry flag value '{geometry}' for species '{label}'. "
                "Flag value must be 0, 1, or 2.")
            return

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
        self.max_loglevel = logging.NOTSET
        self.handler = ParserLogger(self)
        logger.addHandler(self.handler)
        self.processed_units = False
        self.energy_units = 'cal/mol'  # for the current REACTIONS section
        self.output_energy_units = 'cal/mol'  # for the output file
        self.quantity_units = 'mol'  # for the current REACTIONS section
        self.output_quantity_units = 'mol'  # for the output file
        self.motz_wise = None
        self.single_intermediate_temperature = False
        self.skip_undeclared_species = True
        self.exit_on_error = False
        self.permissive = False
        self.verbose = False

        self.elements = []
        self.element_weights = {}  # for custom elements only
        self.species_list = []  # bulk species only
        self.species_dict = {}  # bulk and surface species
        self.surfaces = []
        self.reactions = []
        self.header_lines = []
        self.extra = {}  # for extra entries
        self.files = []  # input file names
        self.raw_lines = []  # lines in current input file; used for error messages
        self.current_range = [0, 0]  # range of current entry; used in error messages

    def __del__(self):
        logger.removeHandler(self.handler)

    def entry(self, where="entry", kind="Error"):
        lines = self.raw_lines[self.current_range[0]:self.current_range[1]+1]
        error_entry = '\n'.join(lines).rstrip()
        return (f'{kind} while reading {where} in {self.files[-1]} starting on '
                f'line {self.current_range[0]+1}:\n"""\n{error_entry}\n"""\n')

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

    def read_NASA7_entry(self, lines, TintDefault, comments):
        """
        Read a thermodynamics entry for one species in a Chemkin-format file
        (consisting of two 7-coefficient NASA polynomials). Returns the label of
        the species, the thermodynamics model as a :class:`Nasa7` object, and
        the elemental composition of the species.

        For more details on this format, see :ref:`sec-debugging-chemkin`.
        """
        identifier = lines[0][0:24].split(maxsplit=1)
        species = identifier[0].strip()

        if len(identifier) > 1:
            note = identifier[1]
        else:
            note = ''

        comments = '\n'.join(c for c in comments if c.strip())
        if comments and note:
            note = '\n'.join((note, comments))
        elif comments:
            note = comments

        # Normal method for specifying the elemental composition
        composition = self.parse_composition(lines[0][24:44], 4, 5)

        # Chemkin-style extended elemental composition: additional lines
        # indicated by '&' continuation character on preceding lines. Element
        # names and abundances are separated by whitespace (not fixed width)
        if lines[0].endswith('&'):
            complines = []
            for i in range(len(lines)-1):
                if lines[i].endswith('&'):
                    complines.append(lines[i+1])
                else:
                    break
            lines = [lines[0]] + lines[i+1:]
            comp = ' '.join(line.rstrip('&') for line in complines).split()
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
            logger.error(self.entry('thermo entry') + "Error parsing elemental "
                         f"composition for species '{species}'.")
            return species, None, {}

        for symbol in composition.keys():
            # Some CHEMKIN input files may have quantities of elements with
            # more than 3 digits. This violates the column-based input format
            # standard, so the entry cannot be read and we need to raise a
            # more useful error message.
            if any(map(str.isdigit, symbol)) and symbol not in self.elements:
                logger.error(self.entry('thermo entry') +
                             "Error parsing elemental composition for species thermo "
                             "entry. Element amounts\ncan have no more than 3 digits.")
                return species, None, {}

        # Extract the NASA polynomial coefficients
        # Remember that the high-T polynomial comes first!
        Tmin = fortFloat(lines[0][45:55])
        Tmax = fortFloat(lines[0][55:65])
        if self.single_intermediate_temperature:
            # Intermediate temperature is shared across all species, except if the
            # species only has one temperature range
            Tint = TintDefault if Tmin < TintDefault < Tmax else None
        else:
            # Non-default intermediate temperature can be provided
            try:
                Tint = fortFloat(lines[0][65:75])
            except ValueError:
                Tint = TintDefault if Tmin < TintDefault < Tmax else None

        high_coeffs = [fortFloat(lines[i][j:k])
                       for i,j,k in [(1,0,15), (1,15,30), (1,30,45), (1,45,60),
                                     (1,60,75), (2,0,15), (2,15,30)]]
        low_coeffs = [fortFloat(lines[i][j:k])
                      for i,j,k in [(2,30,45), (2,45,60), (2,60,75), (3,0,15),
                                    (3,15,30), (3,30,45), (3,45,60)]]

        # Cases where only one temperature range is needed
        if Tint == Tmin or Tint == Tmax or high_coeffs == low_coeffs:
            Tint = None

        # Duplicate the valid set of coefficients if only one range is provided
        if Tint is None:
            if all(c == 0 for c in low_coeffs):
                # Use the first set of coefficients if the second is all zeros
                coeffs = high_coeffs
            elif all(c == 0 for c in high_coeffs):
                # Use the second set of coefficients if the first is all zeros
                coeffs = low_coeffs
            elif high_coeffs == low_coeffs:
                # If the coefficients are duplicated, that's fine too
                coeffs = low_coeffs
            else:
                logger.error(self.entry("thermo entry") +
                    "Only one temperature range defined but two distinct "
                    "sets of coefficients given\nin species thermo entry.")
                return species, None, {}

            thermo = Nasa7(Tmin=Tmin, Tmax=Tmax, Tmid=None, low_coeffs=coeffs,
                           high_coeffs=None, note=note)
        else:
            thermo = Nasa7(Tmin=Tmin, Tmax=Tmax, Tmid=Tint, low_coeffs=low_coeffs,
                           high_coeffs=high_coeffs, note=note)

        return species, thermo, composition

    def parse_nasa7_section(self, lines):
        """
        Parse a THERMO section containing species defined using the 7-coefficient
        NASA polynomial format.

        :param lines:
            A list of ``(line number, section name, line content, comment)`` tuples
        """
        TintDefault = float(lines[1, 2].split()[1])
        lines = lines[2:]

        marker = '1'
        comments = []
        redundant_count = 0
        for i, (input_lineno, _, line, comment) in enumerate(lines):
            comments.append(comment)
            self.current_range[1] = input_lineno
            if len(line) >= 80 and line[79] == marker == '1':
                self.current_range[0] = input_lineno
                entry = [line]
                entry_lines = [i]
                if line.endswith('&'):
                    marker = None
                else:
                    marker = '2'
            elif marker is None:
                entry.append(line)
                entry_lines.append(i)
                if not line.endswith('&'):
                    marker = '2'
            elif len(line) >= 80 and line[79] == marker and marker in ('2', '3'):
                    entry.append(line)
                    entry_lines.append(i)
                    marker = chr(ord(marker) + 1)
            elif len(line) >= 80 and line[79] == marker == '4':
                entry.append(line)
                entry_lines.append(i)
                lines[entry_lines, 1] = 'USED'
                try:
                    name, thermo, comp = self.read_NASA7_entry(
                        entry, TintDefault, comments)
                except ValueError as err:
                    logger.error(self.entry("thermo entry") + str(err))
                    comments = []
                    marker = '1'
                    continue

                if name not in self.species_dict:
                    if self.skip_undeclared_species:
                        logger.debug(f"Skipping unexpected species '{name}' while "
                            "reading thermodynamics entry.")
                        marker = '1'
                        comments = []
                        continue
                    else:
                        # Add a new species entry
                        species = Species(label=name)
                        self.species_dict[name] = species
                        self.species_list.append(species)
                else:
                    species = self.species_dict[name]

                # use the first set of thermo data found
                if species.thermo is not None:
                    redundant_count += 1
                    if redundant_count <= 5 or self.verbose:
                        if self.permissive:
                            logger.warning("Ignoring redundant thermo data for species "
                                f"'{name}' starting on line {input_lineno+1} of "
                                f"{self.files[-1]}.")
                        else:
                            logger.error(self.entry("thermo entry") +
                                f"Found additional thermo entry for species '{name}'. "
                                "Run ck2yaml again with the\n'--permissive' option to "
                                "ignore this redundant entry.")
                else:
                    species.thermo = thermo
                    species.composition = comp
                comments = []
                marker = '1'
            elif line.strip():
                marker = '1'
                comments = []

        if redundant_count > 5 and not self.verbose:
            kind = "warnings" if self.permissive else "errors"
            logger.warning(f"Suppressed {redundant_count - 5} additional {kind} "
                "about redundant thermo data.\nRun ck2yaml again with the "
                f"'--verbose' option to see all {kind}.")

        # Check for blocks of consecutive lines that couldn't be parsed as a valid
        # thermo entry
        unexpected = lines[lines[:,1] != "USED"]
        if not any(unexpected[:,2]):
            return

        blocks = []
        start = 0
        i = -1
        for i in range(len(unexpected) - 1):
            if unexpected[i+1][0] - unexpected[i][0] > 1:
                blocks.append((start, i+1))
                start = i+1
        blocks.append((start, i+2))
        for block in blocks:
            entry = unexpected[block[0]:block[1]]
            if any(entry[:,2]):
                self.current_range = [entry[0,0], entry[-1,0]]
                if self.permissive:
                    logger.warning(self.entry("thermo data", "Unparsable lines") +
                                   "Lines could not be parsed as a NASA7 entry.")
                else:
                    logger.error(self.entry("thermo data", "Unparsable lines") +
                        "Lines could not be parsed as a NASA7 entry. Run ck2yaml again "
                        "with the\n'--permissive' option to continue without parsing "
                        "these lines.")

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

        comments = '\n'.join(c for c in comments if c.strip())
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
            logger.error(self.entry("thermo entry") + str(err))

        thermo = Nasa9(data=polys, note=note, parser=self)

        return species, thermo, composition

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
                    probable_species = re.sub(r'\+?(.+?)(\+\d*)?', r'\1', token)
                    logger.error(self.entry("reaction") +
                        f"Unexpected token '{token}' in reaction expression "
                        f"'{original_reaction.strip()}'.\nMay be due to undeclared "
                        f"species '{probable_species}'.")
                    return Reaction(), None

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
            logger.error(self.entry("reaction") +
                "Failed to find reactant/product delimiter in reaction expression "
                f"'{original_reaction.strip()}'.")

        # Create a new Reaction object for this reaction
        reaction = Reaction(reactants=[], products=[], reversible=reversible)

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
            logger.error(self.entry("reaction") +
                f"Third bodies do not match: '{third_body_name_r}' and "
                f"'{third_body_name_p}'.")

        if photon_r:
            logger.error(self.entry("reaction") + "Reactant photon not supported.")
        if photon_p and reversible:
            if self.permissive:
                logger.warning(self.entry("reaction", "Issue") +
                    "Found a reversible reaction containing a product photon. "
                    "Converting to an\nirreversible reaction with the photon removed.")
            else:
                logger.error(self.entry("reaction") +
                    "Found a reversible reaction containing a product photon. "
                    "To automatically\nconvert this reaction to an irreversible "
                    "reaction with the photon removed,\nrun ck2yaml again with the "
                    "'--permissive' option.")
            reaction.reversible = False

        reaction.third_body = third_body_name_r

        # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
        # This assumes elementary kinetics for all reactions
        rStoich = sum(r[0] for r in reaction.reactants) + (1 if third_body else 0)
        if rStoich < 1:
            logger.error(self.entry("reaction") +
                         "No reactant species found in reaction equation.")

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
                                           reversible=False)

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
                third_body = False # strip optional third-body collider
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
                logger.error(self.entry("reaction") + f"Unparsable line: '{line}'.")

        # Decide which kinetics to keep and store them on the reaction object.
        # At most one of the special cases should be true
        tests = [cheb_coeffs, pdep_arrhenius, low_rate, high_rate, third_body,
                 surface]
        if sum(bool(t) for t in tests) > 1:
            logger.error(self.entry("reaction") +
                "Reaction contains parameters for more than one reaction type.")
            return reaction, revReaction

        if cheb_coeffs:
            if Tmin is None or Tmax is None:
                logger.error(self.entry("reaction") +
                             "Missing TCHEB entry for Chebyshev reaction")
                return reaction, revReaction
            if Pmin is None or Pmax is None:
                logger.error(self.entry("reaction") +
                             "Missing PCHEB entry for Chebyshev reaction")
                return reaction, revReaction
            if len(cheb_coeffs) != degreeT * degreeP:
                logger.error(self.entry("reaction") +
                    "Incorrect number of Chebyshev coefficients. Expected "
                    f"{degreeT}*{degreeP} = {degreeP*degreeT} but got "
                    f"{len(cheb_coeffs)}.")
                return reaction, revReaction
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
            logger.error(self.entry("reaction") +
                "Reaction equation implies pressure dependence but no alternate rate\n"
                "parameters (such as HIGH or LOW) were given.")
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
        try:
            yaml_ = yaml.YAML(typ="rt")
            with open(path, 'rt', encoding="utf-8") as stream:
                yml = yaml_.load(stream)
        except yaml.constructor.ConstructorError:
            with open(path, "rt", encoding="utf-8") as stream:
                # Ensure that the loader remains backward-compatible with legacy
                # ruamel.yaml versions (prior to 0.17.0).
                yml = yaml.round_trip_load(stream)

        # do not overwrite reserved field names
        reserved = {'generator', 'input-files', 'cantera-version', 'date',
                    'units', 'phases', 'species', 'reactions'}
        reserved &= set(yml.keys())
        if reserved:
            logger.error(f"The YAML file '{Path(path).name}' provided as '--extra' "
                         "input\nmust not redefine reserved field name(s): "
                         + ", ".join(f'{item!r}' for item in reserved))
            return

        # replace header lines
        if 'description' in yml:
            if isinstance(yml['description'], str):
                if self.header_lines:
                    self.header_lines += ['']
                self.header_lines += yml.pop('description').split('\n')
            else:
                logger.error("The alternate description provided in "
                             f"'{Path(path).name}'\nneeds to be a string.")

        # remainder
        self.extra = yml

    def load_data_file(self, path: str, load_method: Callable, kind: str, **kwargs):
        if not path:
            return

        self.files.append(Path(path).name)
        path = os.path.expanduser(path)
        if not os.path.exists(path):
            if self.exit_on_error:
                logger.error(f"Missing input file: {path!r}")
                sys.exit(1)
            else:
                raise IOError(f"Missing input file: {path!r}")
        try:
            load_method(path, **kwargs)
        except Exception as err:
            logger.warning(f"\nERROR: Unable to parse '{path}' near line "
                           f"{self.current_range[0]+1}:\n{err}\n")
            raise

    def load_chemkin_file(self, path, surface=False):
        """
        Load a Chemkin-format input file from ``path`` on disk.
        """
        lines = []
        self.raw_lines = Path(path).read_text(errors='ignore').splitlines()
        for i, line in enumerate(self.raw_lines):
            if '!' in line:
                line, comment = line.split('!', 1)
            else:
                comment = ''

            lines.append((i, None, line, comment))

        lines = np.array(lines, dtype=object)

        # Identify file sections
        current_section = None
        sections = {
            "ELEMENTS": re.compile(r"\s*ELEM(?:ENTS)?\b(.*)", re.I),
            "SPECIES": re.compile(r"\s*SPEC(?:IES)?\b(.*)", re.I),
            "SITE": re.compile(r"\s*SITE\b(.*)", re.I),
            "THERMO NASA9": re.compile(r"\s*THER(?:M|MO)\s+NASA9(.*)", re.I),
            "THERMO": re.compile(r"\s*THER(?:M|MO)\b(.*)", re.I),
            "REACTIONS": re.compile(r"\s*REAC(?:TION|TIONS)?\b(.*)", re.I),
            "TRANSPORT": re.compile(r"\s*TRAN(?:SPORT)?\b(.*)", re.I),
        }
        found_sections = set()
        in_header = True
        for i, row in enumerate(lines):
            line = row[2]
            for section, pattern in sections.items():
                if m := pattern.match(line):
                    if current_section and not self.permissive:
                        logger.warning(f"{current_section} section implicitly ended "
                            f"by start of {section} section on line {row[0]+1} of "
                            f"{self.files[-1]}.")
                    current_section = section
                    in_header = False
                    found_sections.add(section)
                    lines[i:, 1] = section
                    line = lines[i, 2] = m.group(1)
                    del sections[section]
                    break
            if current_section in ("SPECIES", "ELEMENTS", "SITE"):
                if m := re.match(r"(.*)\s?END\b", line, re.I):
                    current_section = None
                    lines[i+1:, 1] = None
                    lines[i, 2] = m.group(1)
            elif m := re.match(r"\s*END\b", line, re.I):
                current_section = None
                lines[i+1:, 1] = None
                lines[i, 2] = ""
            elif in_header:
                lines[i:, 1] = "HEADER"
                if row[3].strip():
                    found_sections.add("HEADER")
            elif current_section is None and line.strip():
                self.current_range = [row[0]] * 2
                logger.error(self.entry() +
                    f"Section starts with unrecognized keyword '{line}'")
                break

        if "ELEMENTS" in found_sections:
            self.parse_elements_section(lines[lines[:,1] == "ELEMENTS"])
        if "SPECIES" in found_sections:
            self.parse_species_section(lines[lines[:,1] == "SPECIES"])
        if "SITE" in found_sections:
            if not surface:
                logger.error(self.entry("section") +
                    "Surface phase input file must be specified using the "
                    "'--surface' option\nwith the bulk (gas) phase input "
                    "specified using the '--input' option.")
                return
            self.parse_site_section(lines[lines[:,1] == "SITE"])
        if "THERMO NASA9" in found_sections:
            self.parse_nasa9_section(lines[lines[:,1] == "THERMO NASA9"])
        if "THERMO" in found_sections:
            self.parse_nasa7_section(lines[lines[:,1] == "THERMO"])
        if "REACTIONS" in found_sections:
            self.parse_reactions_section(lines[lines[:,1] == "REACTIONS"], surface)
        if "TRANSPORT" in found_sections:
            self.parse_transport_section(lines[lines[:,1] == "TRANSPORT"])

        if "HEADER" in found_sections:
            indent = 80
            header = lines[lines[:,1] == "HEADER", 3]
            for row in header:
                if row.strip():
                    indent = min(indent, re.search('[^ ]', row).start())
            for row in header:
                self.header_lines.append(row[indent:])

    def parse_elements_section(self, lines):
        """
        Parse the ELEMENTS section of a Chemkin-format input file

        :param lines:
            A list of ``(line number, section name, line content, comment)`` tuples
        """
        tokens = []
        for lineno, _, line, comment in lines:
            # Normalize custom atomic weights
            line = re.sub(r'\s*/\s*([0-9\.EeDd+-]+)\s*/', r'/\1/ ', line)
            tokens.extend(line.split())

        for element_string in tokens:
            if '/' in element_string:
                name, weight, _ = element_string.split('/')
                weight = fortFloat(weight)
                name = name.capitalize()
                self.elements.append(name)
                self.element_weights[name] = weight
            else:
                self.elements.append(element_string.capitalize())

    def parse_species_section(self, lines):
        """
        Parse the SPECIES section of a Chemkin-format input file

        :param lines:
            A list of ``(line number, section name, line content, comment)`` tuples
        """
        comments = {}
        species = []
        for line, comment in lines[:, 2:]:
            line_species = line.split()
            if len(line_species) == 1 and comment:
                comments[line_species[0]] = comment
            species.extend(line_species)

        redundant_count = 0
        for token in species:
            if token in self.species_dict:
                redundant_count += 1
                species = self.species_dict[token]
                if redundant_count > 5:
                    continue
                if self.permissive:
                    logger.warning("Ignoring redundant declaration for "
                                   f"species '{species}'")
                else:
                    logger.error(f"Found multiple declarations for species "
                        f"'{species}'. Run ck2yaml again with the\n'--permissive' "
                        "option to ignore the extra declarations.")
            else:
                species = Species(label=token)
                if token in comments:
                    species.note = comments[token]
                self.species_dict[token] = species
                self.species_list.append(species)

        if redundant_count > 5 and not self.verbose:
            kind = "warnings" if self.permissive else "errors"
            logger.warning(f"Suppressed {redundant_count - 5} additional {kind} about "
                "redundant species declarations.\nRun ck2yaml again with the "
                f"'--verbose' option to see all {kind}.")

    def parse_site_section(self, lines):
        """
        Parse the SITE (surface species) section of a Chemkin-format input file

        :param lines:
            A list of ``(line number, section name, line content, comment)`` tuples
        """
        # First line contains optional parameters, followed by species declarations
        tokens = lines[0,2].split() or ['']
        self.current_range = [lines[0,0]] * 2
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
            logger.error(self.entry("SITE section") +
                         "SITE section defined with no site density")
        self.surfaces.append(Surface(name=surf_name, site_density=site_density))
        surf = self.surfaces[-1]

        # Get the rest of the species declarations
        for line in lines[1:]:
            tokens.extend(line[2].split())

        # List of species identifiers for surface species
        for token in tokens:
            if token.count('/') == 2:
                # species occupies a specific number of sites
                token, sites, _ = token.split('/')
                sites = float(sites)
            else:
                sites = None
            if token in self.species_dict:
                species = self.species_dict[token]
                if self.permissive:
                    logger.warning("Ignoring redundant declaration for "
                                   f"species '{species}'")
                else:
                    logger.error(f"Found multiple declarations for species "
                        f"'{species}'. Run ck2yaml again with the\n'--permissive' "
                        "option to ignore the extra declarations.")
            else:
                species = Species(label=token, sites=sites)
                self.species_dict[token] = species
                surf.species_list.append(species)

    def parse_nasa9_section(self, lines):
        """
        Parse a THERMO section containing species defined using the 9-coefficient
        NASA polynomial format.

        :param lines:
            A list of ``(line number, section name, line content, comment)`` tuples
        """
        # special case if (redundant) temperature ranges are given as the first line
        try:
            s = lines[1, 2].split()
            float(s[0]), float(s[1]), float(s[2])
            lines[1, 2] = ''
        except (IndexError, ValueError):
            pass

        entryLength = 0
        redundant_count = 0
        entry = []
        # Gather comments on lines preceding and within this entry
        comments = []

        for lineno, _, line, comment in lines:
            comments.append(comment)
            if not line:
                continue
            self.current_range[-1] = lineno
            entry.append(line)
            if len(entry) == 1:
                self.current_range[0] = lineno
            elif len(entry) == 2:
                entryLength = 2 + 3 * int(line.split()[0])

            if len(entry) == entryLength:
                label, thermo, comp = self.read_NASA9_entry(entry, comments)
                comments = []
                entry = []
                if label not in self.species_dict:
                    if self.skip_undeclared_species:
                        logger.debug(f'Skipping unexpected species "{label}" while '
                                     'reading thermodynamics entry.')
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
                    redundant_count += 1
                    if redundant_count > 5:
                        continue
                    if self.permissive:
                        logger.warning("Ignoring redundant thermo data for species "
                            f"'{label}' starting on line {self.current_range[0]+1} of "
                            f"{self.files[-1]}.")
                    else:
                        logger.error(self.entry("thermo entry") +
                            f"Found additional thermo entry for species '{label}'. Run "
                            "ck2yaml again with the\n'--permissive' option to ignore "
                            "this redundant entry.")
                else:
                    species.thermo = thermo
                    species.composition = comp

        if redundant_count > 5 and not self.verbose:
            kind = "warnings" if self.permissive else "errors"
            logger.warning(f"Suppressed {redundant_count - 5} additional {kind} about "
                "redundant thermo data.\nRun ck2yaml again with the '--verbose' option "
                f"to see all {kind}.")

    def parse_reactions_section(self, lines, surface):
        """
        Parse the REACTIONS section of a Chemkin-formatted input file.

        :param lines:
            A list of ``(line number, section name, line content, comment)`` tuples
        :param surface:
            Boolean indicating whether this is a surface reaction mechanism
        """
        self.current_range = [lines[0,0]] * 2
        # Reactions section
        for token in lines[0,2].upper().split():
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
                logger.error(self.entry("section") +
                    f"Unrecognized token {token!r} on REACTIONS line")

        lines = lines[1:]
        self.processed_units = True
        kineticsList = []
        commentsList = []
        startLines = []
        kinetics = ''
        comments = ''
        reactions = self.surfaces[-1].reactions if surface else self.reactions

        for lineno, _, line, comment in lines:
            if '=' in line:
                # Finish previous record
                if comment:
                    # End of line comment belongs with this reaction
                    comments += comment + '\n'
                    comment = ''
                kineticsList.append(kinetics)
                commentsList.append(comments)
                startLines.append(lineno)
                kinetics = ''
                comments = ''

            if line.strip():
                kinetics += line + '\n'
            if comment:
                comments += comment + '\n'

        # Don't forget the last reaction!
        if kinetics.strip() != '':
            kineticsList.append(kinetics)
            commentsList.append(comments)

        # We don't actually know whether comments belong to the previous or next
        # reaction, but to keep them positioned correctly, we associate them with the
        # next reaction. A comment after the last reaction is associated with that
        # reaction
        if kineticsList and kineticsList[0] == '':
            kineticsList.pop(0)
            final_comment = commentsList.pop()
            if final_comment and commentsList[-1]:
                commentsList[-1] = commentsList[-1] + '\n' + final_comment
            elif final_comment:
                commentsList[-1] = final_comment

        # We look for species including the next permissible character. '\n' is appended
        # to the reaction string to identify the last species in the reaction string.
        # Checking this character is necessary to correctly identify species with names
        # ending in '+' or '='.
        self.species_tokens = set()
        for next_char in ('<', '=', '(', '+', '\n'):
            self.species_tokens.update(k + next_char for k in self.species_dict)
        self.other_tokens = {'M': 'third-body', 'm': 'third-body', '(+M)': 'falloff3b',
                             '(+m)': 'falloff3b', '<=>': 'equal', '=>': 'equal',
                             '=': 'equal', 'HV': 'photon', 'hv': 'photon'}
        self.other_tokens.update(('(+{})'.format(k), 'falloff3b: {}'.format(k))
                                 for k in self.species_dict)
        self.Slen = max(map(len, self.other_tokens))

        for kinetics, comment, line_number in zip(kineticsList, commentsList, startLines):
            self.current_range = [line_number, line_number + kinetics.count('\n') - 1]
            try:
                reaction, revReaction = self.read_kinetics_entry(kinetics, surface)
            except Exception as e:
                logger.error(self.entry("reaction") + str(e))
                reaction, revReaction = Reaction(), None
            reaction.line_number = line_number
            reaction.comment = comment
            reactions.append(reaction)
            if revReaction is not None:
                revReaction.line_number = line_number
                reactions.append(revReaction)

        for index, reaction in enumerate(reactions):
            reaction.index = index + 1

    def load_transport_file(self, path):
        self.raw_lines = Path(path).read_text(errors='ignore').splitlines()
        lines = []
        for i, line in enumerate(self.raw_lines):
            if '!' in line:
                line, comment = line.split('!', 1)
            else:
                comment = ''

            lines.append((i, None, line.rstrip(), comment.rstrip()))
        lines = np.array(lines, dtype=object)
        self.parse_transport_section(lines)

    def parse_transport_section(self, lines):
        """
        Parse Chemkin-format transport data in ``lines`` and add that transport data
        to the previously-loaded species.

        :param lines:
            A list of ``(line number, section name, line content, comment)`` tuples
        """

        redundant_count = 0
        for lineno, _, line, comment in lines:
            self.current_range = [lineno] * 2
            data = line.strip().split()
            if not data:
                continue

            speciesName = data[0]
            if speciesName in self.species_dict:
                if len(data) != 7:
                    logger.error(self.entry("transport data") + "6 transport "
                        f"parameters were expected, but found {len(data)-1}.")
                    self.species_dict[speciesName].transport = False
                    continue

                if self.species_dict[speciesName].transport is None:
                    self.species_dict[speciesName].transport = TransportData(self, *data, note=comment)
                else:
                    redundant_count += 1
                    if redundant_count > 5 and not self.verbose:
                        continue
                    if self.permissive:
                        logger.warning('Ignoring duplicate transport data for species '
                            f'"{speciesName}" on line {lineno+1} of "{self.files[-1]}".')
                    else:
                        logger.error(self.entry("transport data") +
                            f"Duplicate transport data for species '{speciesName}'. "
                            "Run ck2yaml again with the\n'--permissive' option to "
                            "ignore this redundant entry.")

        if redundant_count > 5 and not self.verbose:
            kind = "warnings" if self.permissive else "errors"
            logger.warning(f"Suppressed {redundant_count - 5} additional {kind} "
                "about duplicate transport data.\nRun ck2yaml again with the "
                f"'--verbose' option to see all {kind}.")

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
            desc = '\n'.join(self.header_lines)
            desc = desc.strip('\n')
            desc = textwrap.dedent(desc)
            if desc.strip():
                emitter.dump({'description': yaml.scalarstring.PreservedScalarString(desc)}, dest)

            # Additional information regarding conversion
            metadata = BlockMap([
                ("generator", "ck2yaml"),
                ("input-files", FlowList(self.files)),
                ("cantera-version", "3.2.0a1"),
                ("date", formatdate(localtime=True)),
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
                phase['adjacent-phases'] = FlowList([name])
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
            speciesMap = BlockMap([('species', self.all_species)])
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
                     out_name=None, single_intermediate_temperature=False, quiet=False,
                     permissive=None, verbose=False, exit_on_error=False):

        # Direct log messages to the console
        loghandler = logging.StreamHandler(sys.stdout)
        loghandler.setFormatter(ErrorFormatter())
        logger.handlers.clear()
        logger.addHandler(loghandler)
        logger.setLevel(logging.INFO)
        logger.propagate = False

        parser = Parser()
        parser.verbose = verbose
        parser.exit_on_error = exit_on_error
        parser.single_intermediate_temperature = single_intermediate_temperature
        if quiet:
            logger.setLevel(level=logging.ERROR if permissive else logging.WARNING)
        elif verbose:
            logger.setLevel(level=logging.DEBUG)
        else:
            logger.setLevel(level=logging.INFO)

        if permissive is not None:
            parser.permissive = permissive

        parser.load_data_file(input_file, parser.load_chemkin_file, "input")
        if surface_file and not input_file:
            msg = "Cannot specify a surface mechanism without a gas phase"
            if exit_on_error:
                logger.warning(f"\nERROR: {msg}\n")
                sys.exit(1)
            else:
                raise InputError(msg)

        parser.load_data_file(surface_file, parser.load_chemkin_file,
                              "surface input", surface=True)

        parser.skip_undeclared_species = bool(input_file)
        parser.load_data_file(thermo_file, parser.load_chemkin_file, "thermo")

        parser.all_species = list(parser.species_list)
        for surf in parser.surfaces:
            parser.all_species.extend(surf.species_list)

        for species in parser.all_species:
            if species.composition is None:
                logger.error(f'No thermo data found for species {species.label!r}')

        parser.load_data_file(transport_file, parser.load_transport_file, "transport")
        if transport_file:
            # Transport validation: make sure all species have transport data
            for s in parser.species_list:
                if s.transport is None:
                    logger.error(f"No transport data for species '{s}'.")

        parser.load_data_file(extra_file, parser.load_extra_file, "input")

        if parser.max_loglevel >= logging.ERROR:
            if exit_on_error:
                logger.error("Unable to convert mechanism due to errors. Please check\n"
                    "https://cantera.org/stable/userguide/"
                    "ck2yaml-tutorial.html#debugging-common-errors-in-ck-files"
                    "\nfor the correct Chemkin syntax.")
                sys.exit(1)
            else:
                raise InputError('\n'.join(parser.handler.errors))

        # Write output file
        if out_name:
            out_name = os.path.expanduser(out_name)
        else:
            out_name = os.path.splitext(input_file)[0] + '.yaml'

        if not input_file:
            phase_name = None

        surface_names = parser.write_yaml(name=phase_name, out_name=out_name)
        if not quiet:
            nSpecies = len(parser.species_list)
            nReactions = len(parser.reactions)
            for surf in parser.surfaces:
                nSpecies += len(surf.species_list)
                nReactions += len(surf.reactions)
            logger.info(f'Wrote YAML mechanism file to {out_name!r}.')
            logger.info(f'Mechanism contains {nSpecies} species and {nReactions} '
                        'reactions.')
        return parser, surface_names

    def show_duplicate_reactions(self, error_message):
        # Find the reaction numbers of the duplicate reactions by looking at
        # the YAML file lines shown in the error message generated by
        # Kinetics::checkDuplicates.
        duplicates = []
        for line in error_message.split('\n'):
            match = re.match('>.*# Reaction ([0-9]+)', line)
            if match:
                duplicates.append(int(match.group(1))-1)

        if not duplicates or len(duplicates) % 2:
            # Something went wrong while parsing the error message, so just
            # display it as-is instead of trying to be clever.
            logger.warning(error_message)
            return

        # Give an error message that references the line numbers in the
        # original input file.
        msg = ["Undeclared duplicate reactions found in the kinetics input file:"]
        for i in range(0, len(duplicates), 2):
            r1 = self.reactions[duplicates[i]]
            r2 = self.reactions[duplicates[i+1]]

            msg.append(f"Line {r1.line_number + 1:d}: {r1}")
            msg.append(f"Line {r2.line_number + 1:d}: {r2}\n")

        logger.warning("\n".join(msg))
        logger.debug(error_message)


def convert(input_file, thermo_file=None, transport_file=None,
            surface_file=None, phase_name='gas', extra_file=None,
            out_name=None, single_intermediate_temperature=False, quiet=False,
            permissive=None, verbose=False):
    _, surface_names = Parser.convert_mech(
        input_file, thermo_file, transport_file, surface_file, phase_name, extra_file,
        out_name, single_intermediate_temperature, quiet, permissive, verbose)
    return surface_names


def main(argv=None):
    """Parse command line arguments and pass them to `Parser.convert_mech`."""
    parser = create_argparser()
    if argv is None and len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(argv)

    input_file = args.input
    thermo_file = args.thermo
    out_name = args.output

    if not input_file and not thermo_file:
        logger.error('At least one of the arguments "--input" or "--thermo" '
                     'must be provided.\nRun with option "--help" to see usage help.')
        sys.exit(1)

    if out_name.endswith('.yaml') or out_name.endswith('.yml'):
        pass
    elif out_name:
        out_name += '.yaml'
    elif input_file:
        out_name = os.path.splitext(input_file)[0] + '.yaml'
    else:
        out_name = os.path.splitext(thermo_file)[0] + '.yaml'

    parser, surfaces = Parser.convert_mech(input_file, thermo_file,
            args.transport, args.surface, args.name, args.extra, out_name,
            args.single_intermediate_temperature, args.quiet, args.permissive,
            args.verbose, True)

    if not input_file:
        # Can't validate input files that don't define a phase
        return

    if args.no_validate:
        return

    # do full validation by importing the resulting mechanism
    try:
        from cantera import Solution, Interface
    except ImportError:
        logger.warning('WARNING: Unable to import Cantera Python module. '
                       'Output mechanism has not been validated')
        sys.exit(0)

    try:
        logger.info('Validating mechanism...')
        gas = Solution(out_name)
        for surf_name in surfaces:
            phase = Interface(out_name, surf_name, [gas])
        logger.info('PASSED')
    except RuntimeError as e:
        logger.info('FAILED')
        msg = str(e)
        if 'Undeclared duplicate reactions' in msg:
            parser.show_duplicate_reactions(msg)
        else:
            logger.warning(e)
        sys.exit(1)


def create_argparser():
    parser = argparse.ArgumentParser(
        description=(
            "Convert Chemkin-format mechanisms to Cantera YAML input files"),
        epilog=textwrap.dedent(
            """
            Example::

                ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat

            If the **ck2yaml** script is not on your path but the Cantera Python module is,
            **ck2yaml** can also be invoked by running::

                python -m cantera.ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat

            In both cases, the equal signs in the options are optional.
            """),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-p", "--permissive", action="store_true", default=False,
        help=("This option allows certain recoverable parsing errors (for example, "
              "duplicate thermo or transport data) to be ignored."))
    parser.add_argument(
        "-q", "--quiet", action="store_true", default=False,
        help="Suppresses warning messages, such as those about duplicate thermo data.")
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Show additional logging output, such as messages about unused species")
    parser.add_argument(
        "--input", default="",
        help=("Chemkin-format chemistry input file, containing a list of all the "
              "element names that are used, a list of all the species names, and a "
              "list of all the reactions to be considered between the species. This "
              "file can also optionally contain species THERMO and TRANSPORT data."))
    parser.add_argument(
        "--thermo", default="",
        help=("If the INPUT file does not contain THERMO data, a separate file "
              "containing thermodynamic information must be specified. If no INPUT "
              "file is provided, THERMO data are converted to a YAML file containing "
              "only species definitions (which can be referenced from phase "
              "definitions in other input files)."))
    parser.add_argument(
        "--transport", default="",
        help=("If the INPUT file does not contain TRANSPORT data, a separate file "
              "containing transport information may be specified. Transport data are "
              "required for Cantera capabilities that use transport properties such "
              "as one-dimensional flame calculations; they are usually not required "
              "for zero-dimensional reactor simulations."))
    parser.add_argument(
        "--surface", default="",
        help=("For surface mechanisms, the SURFACE file defines surface species "
              "and reactions occurring on the surface. Gas phase species and reactions "
              "should be defined in the INPUT file."))
    parser.add_argument(
        "--extra", default="",
        help=("This option specifies a YAML file which can be used to add to the "
              "**description** field or to define custom fields that are included in "
              "the YAML output."))
    parser.add_argument(
        "--name", default="gas",
        help=("This specifies the name of the phase in the resulting YAML file. The "
              "default is **gas**."))
    parser.add_argument(
        "--single-intermediate-temperature", action="store_true", default=False,
        help=("This option should be used with thermo data where only a single break "
              "temperature is used and the last value in the first line of each "
              "species thermo entry is the molecular weight instead."))
    parser.add_argument(
        "--no-validate", action="store_true", default=False,
        help=("Disables the validation step, where the YAML mechanism is imported in "
              "Cantera to check for errors such as unlabeled duplicate reactions and "
              "discontinuous thermodynamic data."))
    parser.add_argument(
        "--output", default="",
        help=("Specifies the OUTPUT file name. By default, the output file name is the "
              "input file name with the extension changed to **.yaml**."))

    return parser

if __name__ == '__main__':
    main()
