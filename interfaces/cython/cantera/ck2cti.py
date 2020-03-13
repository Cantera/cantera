#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains functions for converting Chemkin-format input files to
Cantera input files (CTI).

.. deprecated:: 2.5

    The CTI input file format is deprecated and will be removed in Cantera 3.0.
    Use `ck2yaml.py` to convert Chemkin-format input files to the YAML format.
"""

from __future__ import print_function

from collections import defaultdict
import logging
import os.path
import sys
import numpy as np
import re
import itertools
import getopt

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

Avogadro = 6.02214129e23  # in molec/mol; value consistent with ct_defs.h.

_open = open
if sys.version_info[0] == 2:
    string_types = (str, unicode)
    def strip_nonascii(s):
        return s.decode('ascii', 'ignore')

    def open(filename, mode, *args):
        return _open(filename, mode, *args)

else:
    string_types = (str,)
    def strip_nonascii(s):
        return s.encode('ascii', 'ignore').decode()

    def open(filename, mode, *args):
        mode = mode.replace('U', '')
        return _open(filename, mode, *args, errors='ignore')


def compatible_quantities(quantity_basis, units):
    if quantity_basis == 'mol':
        return 'molec' not in units
    elif quantity_basis == 'molec':
        return 'molec' in units or 'mol' not in units
    else:
        raise ValueError('Unknown quantity basis: "{0}"'.format(quantity_basis))


class InputParseError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin-format
    mechanism files. Pass a string describing the circumstances that caused
    the exceptional behavior.
    """
    pass


class Species(object):
    def __init__(self, label, sites=None):
        self.label = label
        self.thermo = None
        self.transport = None
        self.note = None
        self.sites = sites
        self.composition = None

    def __str__(self):
        return self.label

    def to_cti(self, indent=0):
        lines = []
        atoms = ' '.join('{0}:{1}'.format(*a)
                         for a in sorted(self.composition.items()))

        prefix = ' '*(indent+8)

        lines.append('species(name={0!r},'.format(self.label))
        lines.append(prefix + 'atoms={0!r},'.format(atoms))
        if self.thermo:
            lines.append(prefix +
                         'thermo={0},'.format(self.thermo.to_cti(15+indent)))
        if self.transport:
            lines.append(prefix +
                         'transport={0},'.format(self.transport.to_cti(14+indent)))
        if self.note:
            lines.append(prefix + 'note={0!r},'.format(self.note))
        if self.sites is not None:
            lines.append(prefix + 'size={0},'.format(self.sites))

        lines[-1] = lines[-1][:-1] + ')'
        lines.append('')

        return '\n'.join(lines)


class ThermoModel(object):
    """
    A base class for thermodynamics models, containing several attributes
    common to all models:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `Tmin`          ``float``           The minimum temperature at which the model is valid, or ``None`` if unknown or undefined
    `Tmax`          ``float``           The maximum temperature at which the model is valid, or ``None`` if unknown or undefined
    `comment`       ``str``             Information about the model (e.g. its source)
    =============== =================== ========================================

    """

    def __init__(self, Tmin=None, Tmax=None, comment=''):
        if Tmin is not None:
            self.Tmin = Tmin
        else:
            self.Tmin = None
        if Tmax is not None:
            self.Tmax = Tmax
        else:
            self.Tmax = None
        self.comment = comment


class NASA(ThermoModel):
    """
    A single NASA polynomial for thermodynamic data. The `coeffs` attribute
    stores the seven or nine polynomial coefficients
    :math:`\\mathbf{a} = \\left[a_{-2}\\ a_{-1}\\ a_0\\ a_1\\ a_2\\ a_3\\ a_4\\ a_5\\ a_6 \\right]`
    from which the relevant thermodynamic parameters are evaluated via the
    expressions

    .. math:: \\frac{C_\\mathrm{p}(T)}{R} = a_{-2} T^{-2} + a_{-1} T^{-1} + a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

    .. math:: \\frac{H(T)}{RT} = - a_{-2} T^{-2} + a_{-1} T^{-1} \\ln T + a_0 + \\frac{1}{2} a_1 T + \\frac{1}{3} a_2 T^2 + \\frac{1}{4} a_3 T^3 + \\frac{1}{5} a_4 T^4 + \\frac{a_5}{T}

    .. math:: \\frac{S(T)}{R} = -\\frac{1}{2} a_{-2} T^{-2} - a_{-1} T^{-1} + a_0 \\ln T + a_1 T + \\frac{1}{2} a_2 T^2 + \\frac{1}{3} a_3 T^3 + \\frac{1}{4} a_4 T^4 + a_6

    For the 7 coefficient form, the first two coefficients are taken to be zero.
    """

    def __init__(self, coeffs, **kwargs):
        ThermoModel.__init__(self, **kwargs)
        if len(coeffs) not in (7,9):
            raise InputParseError('Invalid number of NASA polynomial coefficients; '
                                  'should be 7 or 9.')
        self.coeffs = coeffs

    def to_cti(self, indent=0):
        prefix = ' '*indent
        if len(self.coeffs) == 7:
            vals = ['{0: 15.8E}'.format(i) for i in self.coeffs]
            lines = ['NASA([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                     prefix+'     [{0}, {1}, {2},'.format(*vals[0:3]),
                     prefix+'      {0}, {1}, {2},'.format(*vals[3:6]),
                     prefix+'      {0}]),'.format(vals[6])]
        else:
            vals = ['{0: 15.9E}'.format(i) for i in self.coeffs]
            lines = ['NASA9([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                     prefix+'      [{0}, {1}, {2},'.format(*vals[0:3]),
                     prefix+'       {0}, {1}, {2},'.format(*vals[3:6]),
                     prefix+'       {0}, {1}, {2}]),'.format(*vals[6:9])]

        return '\n'.join(lines)


class MultiNASA(ThermoModel):
    """
    A set of thermodynamic parameters given by NASA polynomials. This class
    stores a list of :class:`NASA` objects in the `polynomials`
    attribute. When evaluating a thermodynamic quantity, a polynomial that
    contains the desired temperature within its valid range will be used.
    """

    def __init__(self, polynomials=None, **kwargs):
        ThermoModel.__init__(self, **kwargs)
        self.polynomials = polynomials or []

    def to_cti(self, indent=0):
        prefix = ' '*indent
        lines = []
        for i,p in enumerate(self.polynomials):
            if i == 0:
                if len(self.polynomials) == 1:
                    lines.append('({0})'.format(p.to_cti(indent+1)))
                else:
                    lines.append('({0}'.format(p.to_cti(indent+1)))
            elif i != len(self.polynomials)-1:
                lines.append(prefix + ' {0}'.format(p.to_cti(indent+1)))
            else:
                lines.append(prefix + ' {0})'.format(p.to_cti(indent+1)[:-1]))

        return '\n'.join(lines)


class Reaction(object):
    """
    A chemical reaction. The attributes are:

    =================== =========================== ============================
    Attribute           Type                        Description
    =================== =========================== ============================
    `index`             :class:`int`                A unique nonnegative integer index
    `reactants`         :class:`list`               The reactant species (as :class:`Species` objects)
    `products`          :class:`list`               The product species (as :class:`Species` objects)
    `kinetics`          :class:`KineticsModel`      The kinetics model to use for the reaction
    `reversible`        ``bool``                    ``True`` if the reaction is reversible, ``False`` if not
    `duplicate`         ``bool``                    ``True`` if the reaction is known to be a duplicate, ``False`` if not
    `fwdOrders`         ``dict``                    Reaction order (value) for each specified species (key)
    =================== =========================== ============================

    """

    def __init__(self, parser, index=-1, reactants=None, products=None, kinetics=None,
                 reversible=True, duplicate=False, fwdOrders=None,
                 thirdBody=None, ID=''):
        self.parser = parser
        self.index = index
        self.reactants = reactants  # list of (stoichiometry, species) tuples
        self.products = products  # list of (stoichiometry, specis) tuples
        self.kinetics = kinetics
        self.reversible = reversible
        self.duplicate = duplicate
        self.fwdOrders = fwdOrders if fwdOrders is not None else {}
        self.thirdBody = thirdBody
        self.ID = ID
        self.comment = ''

    def _coeff_string(self, coeffs):
        L = []
        for stoichiometry, species in coeffs:
            if stoichiometry != 1:
                L.append('{0} {1}'.format(stoichiometry, species))
            else:
                L.append(str(species))
        expression = ' + '.join(L)
        if self.thirdBody:
            expression += ' (+ {0})'.format(self.thirdBody)

        return expression

    @property
    def reactantString(self):
        return self._coeff_string(self.reactants)

    @property
    def productString(self):
        return self._coeff_string(self.products)

    def __str__(self):
        """
        Return a string representation of the reaction, in the form 'A + B <=> C + D'.
        """
        arrow = ' <=> ' if self.reversible else ' -> '
        return arrow.join([self.reactantString, self.productString])

    def to_cti(self, indent=0):
        arrow = ' <=> ' if self.reversible else ' => '

        kinstr = self.kinetics.to_cti(self.reactantString, arrow,
                                      self.productString, indent)

        k_indent = ' ' * (kinstr.find('(') + 1)

        if self.fwdOrders:
            order = ' '.join('{0}:{1}'.format(k,v)
                             for (k,v) in sorted(self.fwdOrders.items()))
            kinstr = kinstr[:-1] + ",\n{0}order='{1}')".format(k_indent, order)

        if self.ID:
            kinstr = kinstr[:-1] + ",\n{0}id={1!r})".format(k_indent, self.ID)

        options = self.kinetics.options()
        if self.duplicate:
            options.append('duplicate')

        if any((float(x) < 0 for x in self.fwdOrders.values())):
            options.append('negative_orders')
            self.parser.warn('Negative reaction order for reaction {} ({}{}{}).'.format(
                self.index, self.reactantString, arrow, self.productString))

        if len(options) == 1:
            optStr = repr(options[0])
        else:
            optStr = repr(options)

        if options:
            kinstr = kinstr[:-1] + ",\n{0}options={1})".format(k_indent, optStr)

        return kinstr


class KineticsModel(object):
    """
    A base class for kinetics models, containing several attributes common to
    all models:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `Tmin`          :class:`Quantity`   The minimum absolute temperature in K at which the model is valid
    `Tmax`          :class:`Quantity`   The maximum absolute temperature in K at which the model is valid
    `Pmin`          :class:`Quantity`   The minimum absolute pressure in Pa at which the model is valid
    `Pmax`          :class:`Quantity`   The maximum absolute pressure in Pa at which the model is valid
    `comment`       :class:`str`        A string containing information about the model (e.g. its source)
    =============== =================== ========================================

    """

    def __init__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment='',
                 parser=None):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.comment = comment
        self.parser = parser
        self.efficiencies = {}

    def isPressureDependent(self):
        """
        Return ``True`` if the kinetics are pressure-dependent or ``False`` if
        they are pressure-independent. This method must be overloaded in the
        derived class.
        """
        raise InputParseError('Unexpected call to KineticsModel.isPressureDependent();'
                              ' you should be using a class derived from KineticsModel.')

    def to_cti(self, reactantstr, arrow, productstr):
        raise InputParseError('to_cti is not implemented for objects of class {0}'.format(self.__class__.__name__))

    def options(self):
        return []

    def efficiencyString(self):
        return ' '.join('{0}:{1}'.format(mol, eff)
                        for mol, eff in sorted(self.efficiencies.items()))


class KineticsData(KineticsModel):
    """
    A kinetics model based around a set of discrete (high-pressure limit)
    rate coefficients at various temperatures. The attributes are:

    =========== =================== ============================================
    Attribute   Type                Description
    =========== =================== ============================================
    `Tdata`     :class:`Quantity`   The temperatures at which the heat capacity data is provided
    `kdata`     :class:`Quantity`   The rate coefficients in SI units at each temperature in `Tdata`
    =========== =================== ============================================

    """

    def __init__(self, Tdata=None, kdata=None, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.Tdata = Tdata
        self.kdata = kdata

    def isPressureDependent(self):
        """
        Returns ``False`` since KineticsData kinetics are not
        pressure-dependent.
        """
        return False


class Arrhenius(KineticsModel):
    """
    Represent a set of modified Arrhenius kinetics. The kinetic expression has
    the form

    .. math:: k(T) = A \\left( \\frac{T}{T_0} \\right)^b \\exp \\left( - \\frac{E_\\mathrm{a}}{RT} \\right)

    where :math:`A`, :math:`b`, :math:`E_\\mathrm{a}`, and :math:`T_0` are the
    parameters to be set, :math:`T` is absolute temperature, and :math:`R` is
    the gas law constant. The attributes are:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `A`             :class:`Quantity`   The preexponential factor in s^-1, m^3/mol*s, etc.
    `T0`            :class:`Quantity`   The reference temperature in K
    `b`             :class:`Quantity`   The temperature exponent
    `Ea`            :class:`Quantity`   The activation energy in J/mol
    =============== =================== ========================================

    """

    def __init__(self, A=0.0, b=0.0, Ea=0.0, T0=1.0, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.A = A
        self.T0 = T0
        self.b = b
        self.Ea = Ea

    def isPressureDependent(self):
        """
        Returns ``False`` since Arrhenius kinetics are not pressure-dependent.
        """
        return False

    def rateStr(self):
        if compatible_quantities(self.parser.output_quantity_units, self.A[1]):
            A = '{0:e}'.format(self.A[0])
        else:
            A = "({0:e}, '{1}')".format(*self.A)

        if self.Ea[1] == self.parser.output_energy_units:
            Ea = str(self.Ea[0])
        else:
            Ea = "({0}, '{1}')".format(*self.Ea)

        return '[{0}, {1}, {2}]'.format(A, self.b, Ea)

    def options(self):
        if self.A[0] < 0:
            return ['negative_A']
        else:
            return []

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstring = reactantstr + arrow + productstr
        return 'reaction({0!r}, {1})'.format(rxnstring, self.rateStr())


class SurfaceArrhenius(Arrhenius):
    """
    An Arrhenius-like reaction occurring on a surface
    """
    def __init__(self, *args, **kwargs):
        Arrhenius.__init__(self, *args, **kwargs)
        self.coverages = []
        self.is_sticking = False
        self.motz_wise = None

    def rateStr(self):
        if self.is_sticking:
            if self.motz_wise is None:
                return ' stick({0})'.format(Arrhenius.rateStr(self)[1:-1])
            else:
                return ' stick({0}, motz_wise={1})'.format(
                    Arrhenius.rateStr(self)[1:-1], self.motz_wise)
        elif not self.coverages:
            return ' ' + Arrhenius.rateStr(self)

        s = Arrhenius.rateStr(self)
        s = '\n{0}Arrhenius({1},\n{2}coverage=['.format(' '*17, s[1:-1], ' '*27)
        for species,A,m,E in self.coverages:
            # Energy units for coverage modification match energy units for
            # base reaction
            if self.Ea[1] != self.parser.output_energy_units:
                E = (E, self.Ea[1])
            s += '[{0!r}, {1}, {2}, {3}],\n{4}'.format(str(species), A, m, E, ' '*37)

        return s.rstrip()[:-1] + '])'

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstring = reactantstr + arrow + productstr
        return 'surface_reaction({0!r},{1})'.format(rxnstring, self.rateStr())


class PDepArrhenius(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = A(P) T^{b(P)} \\exp \\left[ \\frac{-E_\\mathrm{a}(P)}{RT} \\right]

    where the modified Arrhenius parameters are stored at a variety of pressures
    and interpolated between on a logarithmic scale. The attributes are:

    =============== ================== ============================================
    Attribute       Type               Description
    =============== ================== ============================================
    `pressures`     :class:`list`      The list of pressures in Pa
    `arrhenius`     :class:`list`      The list of :class:`Arrhenius` objects at each pressure
    `highPlimit`    :class:`Arrhenius` The high (infinite) pressure limiting :class:`Arrhenius` expression
    =============== ================== ============================================

    Note that `highPlimit` is not used in evaluating k(T,P).
    """

    def __init__(self, pressures=None, arrhenius=None, highPlimit=None, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.pressures = pressures
        self.arrhenius = arrhenius or []
        self.highPlimit = highPlimit or None

    def isPressureDependent(self):
        """
        Returns ``True`` since PDepArrhenius kinetics are pressure-dependent.
        """
        return True

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstring = reactantstr + arrow + productstr
        lines = ['pdep_arrhenius({0!r},'.format(rxnstring)]
        prefix = ' '*(indent+15)
        template = '[({0}, {1!r}), {2}],'
        for pressure, arrhenius in zip(self.pressures[0], self.arrhenius):
            lines.append(prefix + template.format(pressure,
                                                  self.pressures[1],
                                                  arrhenius.rateStr()[1:-1]))
        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)


class Chebyshev(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: \\log k(T,P) = \\sum_{t=1}^{N_T} \\sum_{p=1}^{N_P} \\alpha_{tp} \\phi_t(\\tilde{T}) \\phi_p(\\tilde{P})

    where :math:`\\alpha_{tp}` is a constant, :math:`\\phi_n(x)` is the
    Chebyshev polynomial of degree :math:`n` evaluated at :math:`x`, and

    .. math:: \\tilde{T} \\equiv \\frac{2T^{-1} - T_\\mathrm{min}^{-1} - T_\\mathrm{max}^{-1}}{T_\\mathrm{max}^{-1} - T_\\mathrm{min}^{-1}}

    .. math:: \\tilde{P} \\equiv \\frac{2 \\log P - \\log P_\\mathrm{min} - \\log P_\\mathrm{max}}{\\log P_\\mathrm{max} - \\log P_\\mathrm{min}}

    are reduced temperature and reduced pressures designed to map the ranges
    :math:`(T_\\mathrm{min}, T_\\mathrm{max})` and
    :math:`(P_\\mathrm{min}, P_\\mathrm{max})` to :math:`(-1, 1)`.
    The attributes are:

    =============== =============== ============================================
    Attribute       Type            Description
    =============== =============== ============================================
    `coeffs`        :class:`list`   Matrix of Chebyshev coefficients
    `kunits`        ``str``         The units of the generated k(T, P) values
    `degreeT`       :class:`int`    The number of terms in the inverse temperature direction
    `degreeP`       :class:`int`    The number of terms in the log pressure direction
    =============== =============== ============================================

    """

    def __init__(self, coeffs=None, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        if coeffs is not None:
            self.coeffs = np.array(coeffs, np.float64)
            self.degreeT = self.coeffs.shape[0]
            self.degreeP = self.coeffs.shape[1]
        else:
            self.coeffs = None
            self.degreeT = 0
            self.degreeP = 0

    def isPressureDependent(self):
        """
        Returns ``True`` since Chebyshev polynomial kinetics are
        pressure-dependent.
        """
        return True

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + arrow + productstr
        prefix = ' '*(indent+19)
        lines = ['chebyshev_reaction({0!r},'.format(rxnstr),
                 prefix + 'Tmin={0.Tmin}, Tmax={0.Tmax},'.format(self),
                 prefix + 'Pmin={0.Pmin}, Pmax={0.Pmax},'.format(self)]
        for i in range(self.degreeT):
            coeffline = ', '.join('{0: 12.5e}'.format(self.coeffs[i,j]) for j in range(self.degreeP))
            if i == 0:
                lines.append(prefix + 'coeffs=[[{0}],'.format(coeffline))
            else:
                lines.append(prefix + '        [{0}],'.format(coeffline))

        lines[-1] = lines[-1][:-1] + '])'
        return '\n'.join(lines)


class ThirdBody(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k(T) [\\ce{M}]

    where :math:`k(T)` is an Arrhenius expression and
    :math:`[\\ce{M}] \\approx P/RT` is the concentration of the third body
    (i.e. the bath gas). A collision efficiency can be used to further correct
    the value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusHigh=None, efficiencies=None, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.arrheniusHigh = arrheniusHigh
        self.efficiencies = {}
        if efficiencies is not None:
            for mol, eff in efficiencies.items():
                self.efficiencies[mol] = eff

    def isPressureDependent(self):
        """
        Returns ``True`` since third-body kinetics are pressure-dependent.
        """
        return True

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' + M' + arrow + productstr + ' + M'
        prefix = ' '*(indent + 20)
        lines = ['three_body_reaction({0!r}, {1},'.format(rxnstr, self.arrheniusHigh.rateStr())]
        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

    def options(self):
        if self.arrheniusHigh.A[0] < 0:
            return ['negative_A']
        else:
            return []


class Falloff(ThirdBody):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k_\\infty(T) \\left[ \\frac{P_\\mathrm{r}}{1 + P_\\mathrm{r}} \\right] F

    where

    .. math::

        P_\\mathrm{r} &= \\frac{k_0(T)}{k_\\infty(T)} [\\ce{M}]

        k_0(T) &= A_0 T^{n_0} \\exp \\left( - \\frac{E_0}{RT} \\right)

        k_\\infty(T) &= A_\\infty T^{n_\\infty} \\exp \\left( - \\frac{E_\\infty}{RT} \\right)

    and :math:`[\\ce{M}] \\approx P/RT` is the concentration of the
    bath gas. The Arrhenius expressions :math:`k_0(T)` and :math:`k_\\infty(T)`
    represent the low-pressure and high-pressure limit kinetics, respectively.
    The former is necessarily one reaction order higher than the latter.
    Several different parameterizations are allowed for the falloff function
    :math:`F(P_r, T)`. A collision efficiency can be used to further correct
    the value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    `F`                                     Falloff function parameterization
    =============== ======================= ====================================
    """
    def __init__(self, arrheniusLow=None, F=None, **kwargs):
        ThirdBody.__init__(self, **kwargs)
        self.arrheniusLow = arrheniusLow
        self.F = F

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + arrow + productstr
        prefix = ' '*(indent + 17)
        lines = ['falloff_reaction({0!r},'.format(rxnstr)]
        lines.append(prefix + 'kf={0},'.format(self.arrheniusHigh.rateStr()))
        lines.append(prefix + 'kf0={0},'.format(self.arrheniusLow.rateStr()))
        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))
        if self.F:
            lines.append(prefix + 'falloff={0},'.format(self.F.to_cti()))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)


class ChemicallyActivated(ThirdBody):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = k_0(T) \\left[ \\frac{1}{1 + P_\\mathrm{r}} \\right] F

    where

    .. math::

        P_\\mathrm{r} &= \\frac{k_0(T)}{k_\\infty(T)} [\\ce{M}]

        k_0(T) &= A_0 T^{n_0} \\exp \\left( - \\frac{E_0}{RT} \\right)

        k_\\infty(T) &= A_\\infty T^{n_\\infty} \\exp \\left( - \\frac{E_\\infty}{RT} \\right)

    and :math:`[\\ce{M}] \\approx P/RT` is the concentration of the bath gas.
    The Arrhenius expressions :math:`k_0(T)` and :math:`k_\\infty(T)`
    represent the low-pressure and high-pressure limit kinetics, respectively.
    The former is necessarily one reaction order higher than the latter. The
    allowable parameterizations for the function *F* are the same as for the
    `Falloff` class. A collision efficiency can be used to further correct the
    value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    `F`                                     Falloff function parameterization
    =============== ======================= ====================================
    """
    def __init__(self, arrheniusLow=None, F=None, **kwargs):
        ThirdBody.__init__(self, **kwargs)
        self.arrheniusLow = arrheniusLow
        self.F = F

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + arrow + productstr
        prefix = ' '*(indent + 30)
        lines = ['chemically_activated_reaction({0!r},'.format(rxnstr)]
        lines.append(prefix + 'kLow={0},'.format(self.arrheniusLow.rateStr()))
        lines.append(prefix + 'kHigh={0},'.format(self.arrheniusHigh.rateStr()))
        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))
        if self.F:
            lines.append(prefix + 'falloff={0},'.format(self.F.to_cti()))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)


class Troe(object):
    """
    For the Troe model the parameter :math:`F` is computed via

    .. math::

        \\log F &= \\left\\{1 + \\left[ \\frac{\\log P_\\mathrm{r} + c}{n - d (\\log P_\\mathrm{r} + c)} \\right]^2 \\right\\}^{-1} \\log F_\\mathrm{cent}

        c &= -0.4 - 0.67 \\log F_\\mathrm{cent}

        n &= 0.75 - 1.27 \\log F_\\mathrm{cent}

        d &= 0.14

        F_\\mathrm{cent} &= (1 - \\alpha) \\exp \\left( -T/T_3 \\right) + \\alpha \\exp \\left( -T/T_1 \\right) + \\exp \\left( -T_2/T \\right)

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `alpha`         :class:`Quantity`       The :math:`\\alpha` parameter
    `T1`            :class:`Quantity`       The :math:`T_1` parameter
    `T2`            :class:`Quantity`       The :math:`T_2` parameter
    `T3`            :class:`Quantity`       The :math:`T_3` parameter
    =============== ======================= ====================================

    """

    def __init__(self, alpha=0.0, T3=0.0, T1=0.0, T2=None):
        self.alpha = alpha
        self.T3 = T3
        self.T1 = T1
        self.T2 = T2

    def to_cti(self):
        if self.T2:
            return 'Troe(A={0.alpha[0]}, T3={0.T3[0]}, T1={0.T1[0]}, T2={0.T2[0]})'.format(self)
        else:
            return 'Troe(A={0.alpha[0]}, T3={0.T3[0]}, T1={0.T1[0]})'.format(self)


class Sri(object):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T, P)` using the
    "SRI" formulation of the blending function :math:`F` using either 3 or
    5 parameters. See `The SRI Falloff Function
    <https://cantera.org/science/reactions.html#sec-sri-falloff>`__.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `A`            ``float``                The :math:`a` parameter
    `B`            ``float``                The :math:`b` parameter
    `C`            ``float``                The :math:`c` parameter
    `D`            ``float``                The :math:`d` parameter
    `E`            ``float``                The :math:`e` parameter
    =============== ======================= ====================================
    """

    def __init__(self, A=0.0, B=0.0, C=0.0, D=1.0, E=0.0):
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def to_cti(self):
        if self.D == 1.0 and self.E == 0.0:
            return 'SRI(A={0.A}, B={0.B}, C={0.C})'.format(self)
        else:
            return 'SRI(A={0.A}, B={0.B}, C={0.C}, D={0.D}, E={0.E})'.format(self)


class TransportData(object):
    geometryFlags = ['atom', 'linear', 'nonlinear']

    def __init__(self, label, geometry, wellDepth, collisionDiameter,
                 dipoleMoment, polarizability, zRot, comment=None):

        try:
            geometry = int(geometry)
        except ValueError:
            raise InputParseError("Bad geometry flag '{0}' for species '{1}', is the flag a float "
                                  "or character? It should be an integer.".format(geometry, label))
        if geometry not in (0, 1, 2):
            raise InputParseError("Bad geometry flag '{0}' for species '{1}'".format(geometry, label))

        self.label = label
        self.geometry = self.geometryFlags[int(geometry)]
        self.wellDepth = float(wellDepth)
        self.collisionDiameter = float(collisionDiameter)
        self.dipoleMoment = float(dipoleMoment)
        self.polarizability = float(polarizability)
        self.zRot = float(zRot)
        self.comment = comment or ''  # @todo: include this in the CTI

    def __repr__(self):
        return ('TransportData({label!r}, {geometry!r}, {wellDepth!r}, '
                '{collisionDiameter!r}, {dipoleMoment!r}, {polarizability!r}, '
                '{zRot!r}, {comment!r})').format(**self.__dict__)

    def to_cti(self, indent=0):
        prefix = ' '*(indent+18)
        lines = ['gas_transport(geom={0!r},'.format(self.geometry),
                 prefix+'diam={0},'.format(self.collisionDiameter),
                 prefix+'well_depth={0},'.format(self.wellDepth)]
        if self.dipoleMoment:
            lines.append(prefix+'dipole={0},'.format(self.dipoleMoment))
        if self.polarizability:
            lines.append(prefix+'polar={0},'.format(self.polarizability))
        if self.zRot:
            lines.append(prefix+'rot_relax={0},'.format(self.zRot))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)


def fortFloat(s):
    """
    Convert a string representation of a floating point value to a float,
    allowing for some of the peculiarities of allowable Fortran representations.
    """
    s = s.strip()
    s = s.replace('D', 'E').replace('d', 'e')
    s = s.replace('E ', 'E+').replace('e ', 'e+')
    return float(s)


def get_index(seq, value):
    """
    Find the first location in *seq* which contains a case-insensitive,
    whitespace-insensitive match for *value*. Returns *None* if no match is
    found.
    """
    if isinstance(seq, string_types):
        seq = seq.split()
    value = value.lower().strip()
    for i, item in enumerate(seq):
        if item.lower() == value:
            return i
    return None


def contains(seq, value):
    if isinstance(seq, string_types):
        return value.lower() in seq.lower()
    else:
        return get_index(seq, value) is not None


class Surface(object):
    def __init__(self, name, density):
        self.name = name
        self.siteDensity = density
        self.speciesList = []
        self.reactions = []


class Parser(object):
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
        self.speciesList = []  # bulk species only
        self.speciesDict = {}  # bulk and surface species
        self.surfaces = []
        self.reactions = []
        self.finalReactionComment = ''
        self.headerLines = []

    def warn(self, message):
        if self.warning_as_error:
            raise InputParseError(message)
        else:
            logging.warning(message)

    def parseComposition(self, elements, nElements, width):
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
                count = int(float(count))
                if count:
                    composition[symbol.capitalize()] = count
            except ValueError:
                pass
        return composition

    def getRateConstantUnits(self, length_dims, length_units, quantity_dims,
                             quantity_units, time_dims=1, time_units='s'):

        units = ''
        if length_dims:
            units += length_units
        if length_dims > 1:
            units += str(length_dims)
        if quantity_dims:
            units += '/' + quantity_units
        if quantity_dims > 1:
            units += str(quantity_dims)
        if time_dims:
            units += '/' + time_units
        if time_dims > 1:
            units += str(time_dims)
        if units.startswith('/'):
            units = '1' + units
        return units

    def addElement(self, element_string):
        if '/' in element_string:
            name, weight, _ = element_string.split('/')
            weight = fortFloat(weight)
            name = name.capitalize()
            self.elements.append(name)
            self.element_weights[name] = weight
        else:
            self.elements.append(element_string.capitalize())

    def readThermoEntry(self, lines, TintDefault):
        """
        Read a thermodynamics entry for one species in a Chemkin-format file
        (consisting of two 7-coefficient NASA polynomials). Returns the label of
        the species, the thermodynamics model as a :class:`MultiNASA` object, the
        elemental composition of the species, and the comment/note associated with
        the thermo entry.
        """
        identifier = lines[0][0:24].split()
        species = identifier[0].strip()

        if len(identifier) > 1:
            note = ''.join(identifier[1:]).strip()
        else:
            note = ''

        # Normal method for specifying the elemental composition
        composition = self.parseComposition(lines[0][24:44], 4, 5)

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
            composition2 = self.parseComposition(elements, len(elements)//10, 10)
            composition.update(composition2)

        if not composition:
            raise InputParseError("Error parsing elemental composition for "
                                  "species '{0}'".format(species))

        # Extract the NASA polynomial coefficients
        # Remember that the high-T polynomial comes first!
        Tmin = fortFloat(lines[0][45:55])
        Tmax = fortFloat(lines[0][55:65])
        try:
            Tint = fortFloat(lines[0][65:75])
        except ValueError:
            Tint = TintDefault

        coeffs_high = [fortFloat(lines[i][j:k])
                       for i,j,k in [(1,0,15), (1,15,30), (1,30,45), (1,45,60),
                                     (1,60,75), (2,0,15), (2,15,30)]]
        coeffs_low = [fortFloat(lines[i][j:k])
                       for i,j,k in [(2,30,45), (2,45,60), (2,60,75), (3,0,15),
                                     (3,15,30), (3,30,45), (3,45,60)]]

        # Duplicate the valid set of coefficients if only one range is provided
        if all(c == 0 for c in coeffs_low) and Tmin == Tint:
            coeffs_low = coeffs_high
        elif all(c == 0 for c in coeffs_high) and Tmax == Tint:
            coeffs_high = coeffs_low

        # Construct and return the thermodynamics model
        thermo = MultiNASA(
            polynomials=[
                NASA(Tmin=(Tmin,"K"), Tmax=(Tint,"K"), coeffs=coeffs_low),
                NASA(Tmin=(Tint,"K"), Tmax=(Tmax,"K"), coeffs=coeffs_high)
            ],
            Tmin=(Tmin,"K"),
            Tmax=(Tmax,"K"),
        )

        return species, thermo, composition, note

    def readNasa9Entry(self, entry):
        """
        Read a thermodynamics `entry` for one species given as one or more
        9-coefficient NASA polynomials, written in the format described in
        Appendix A of NASA Reference Publication 1311 (McBride and Gordon, 1996).
        Returns the label of the species, the thermodynamics model as a
        :class:`MultiNASA` object, the elemental composition of the species, and
        the comment/note associated with the thermo entry.
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

        composition = self.parseComposition(entry[1][10:50], 5, 8)

        polys = []
        totalTmin = 1e100
        totalTmax = -1e100
        try:
            for i in range(N):
                A,B,C = entry[2+3*i:2+3*(i+1)]
                Tmin = fortFloat(A[1:11])
                Tmax = fortFloat(A[11:21])
                coeffs = [fortFloat(B[0:16]), fortFloat(B[16:32]),
                          fortFloat(B[32:48]), fortFloat(B[48:64]),
                          fortFloat(B[64:80]), fortFloat(C[0:16]),
                          fortFloat(C[16:32]), fortFloat(C[48:64]),
                          fortFloat(C[64:80])]
                polys.append(NASA(Tmin=(Tmin,"K"), Tmax=(Tmax,"K"), coeffs=coeffs))
                totalTmin = min(Tmin, totalTmin)
                totalTmax = max(Tmax, totalTmax)
        except (IndexError, ValueError) as err:
            raise InputParseError('Error while reading thermo entry for species {0}:\n{1}'.format(species, err))

        thermo = MultiNASA(polynomials=polys,
                           Tmin=(totalTmin,"K"),
                           Tmax=(totalTmax,"K"))

        return species, thermo, composition, note

    def setupKinetics(self):
        # We look for species including the next permissible character. '\n' is
        # appended to the reaction string to identify the last species in the
        # reaction string. Checking this character is necessary to correctly
        # identify species with names ending in '+' or '='.
        self.species_tokens = set()
        for next_char in ('<', '=', '(', '+', '\n'):
            self.species_tokens.update(k + next_char for k in self.speciesDict)
        self.other_tokens = {'M': 'third-body', 'm': 'third-body',
                             '(+M)': 'falloff3b', '(+m)': 'falloff3b',
                             '<=>': 'equal', '=>': 'equal', '=': 'equal',
                             'HV': 'photon', 'hv': 'photon'}
        self.other_tokens.update(('(+%s)' % k, 'falloff3b: %s' % k) for k in self.speciesDict)
        self.Slen = max(map(len, self.other_tokens))

    def readKineticsEntry(self, entry, surface):
        """
        Read a kinetics `entry` for a single reaction as loaded from a
        Chemkin-format file. Returns a :class:`Reaction` object with the
        reaction and its associated kinetics.
        """

        # Handle non-default units which apply to this entry
        energy_units = self.energy_units
        quantity_units = self.quantity_units
        if 'units' in entry.lower():
            for units in sorted(QUANTITY_UNITS, key=lambda k: -len(k)):
                pattern = re.compile(r'units *\/ *%s *\/' % re.escape(units),
                                     flags=re.IGNORECASE)
                m = pattern.search(entry)
                if m:
                    entry = pattern.sub('', entry)
                    quantity_units = QUANTITY_UNITS[units]
                    break

            for units in sorted(ENERGY_UNITS, key=lambda k: -len(k)):
                pattern = re.compile(r'units *\/ *%s *\/' % re.escape(units),
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
        original_reaction = reaction # for use in error messages

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
                    raise InputParseError('Unexpected token "{0}" in reaction expression "{1}".'.format(token, original_reaction))

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

        if lhs is True:
            raise InputParseError("Failed to find reactant/product delimiter in reaction string.")

        # Create a new Reaction object for this reaction
        reaction = Reaction(reactants=[], products=[], reversible=reversible,
                            parser=self)

        def parseExpression(expression, dest):
            falloff3b = None
            thirdBody = False  # simple third body reaction (non-falloff)
            photon = False
            for stoichiometry, species, kind in expression:
                if kind == 'third-body':
                    thirdBody = True
                elif kind == 'falloff3b':
                    falloff3b = 'M'
                elif kind.startswith('falloff3b:'):
                    falloff3b = kind.split()[1]
                elif kind == 'photon':
                    photon = True
                else:
                    dest.append((stoichiometry, self.speciesDict[species]))

            return falloff3b, thirdBody, photon

        falloff_3b_r, thirdBody, photon_r = parseExpression(reactants, reaction.reactants)
        falloff_3b_p, thirdBody, photon_p = parseExpression(products, reaction.products)

        if falloff_3b_r != falloff_3b_p:
            raise InputParseError('Third bodies do not match: "{0}" and "{1}" in'
                ' reaction entry:\n\n{2}'.format(falloff_3b_r, falloff_3b_p, entry))

        if photon_r:
            raise InputParseError('Reactant photon not supported. '
                                  'Found in reaction:\n{0}'.format(entry.strip()))
        if photon_p and reversible:
            self.warn('Found reversible reaction containing a product photon:'
                '\n{0}\nIf the "--permissive" option was specified, this will '
                'be converted to an irreversible reaction with the photon '
                'removed.'.format(entry.strip()))
            reaction.reversible = False

        reaction.thirdBody = falloff_3b_r

        # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
        # This assumes elementary kinetics for all reactions
        rStoich = sum(r[0] for r in reaction.reactants) + (1 if thirdBody else 0)
        if rStoich < 1:
            raise InputParseError('No reactant species for reaction {1}.'.format(reaction))

        length_dim = 3 * (rStoich - 1)
        quantity_dim = rStoich - 1
        kunits = self.getRateConstantUnits(length_dim, 'cm',
                                           quantity_dim, quantity_units)
        klow_units = self.getRateConstantUnits(length_dim + 3, 'cm',
                                               quantity_dim + 1, quantity_units)

        # The rest of the first line contains Arrhenius parameters
        reaction_type = SurfaceArrhenius if surface else Arrhenius
        arrhenius = reaction_type(
            A=(A, kunits),
            b=b,
            Ea=(Ea, energy_units),
            T0=(1,"K"),
            parser=self
        )

        arrheniusLow = None
        arrheniusHigh = None
        falloff = None
        chebyshev = None
        pdepArrhenius = None
        efficiencies = {}
        chebyshevCoeffs = []
        revReaction = None

        # Note that the subsequent lines could be in any order
        for line in lines[1:]:
            if not line.strip():
                continue  # blank line, or (erased) units field
            tokens = line.split('/')
            parsed = False

            if 'stick' in line.lower():
                parsed = True
                arrhenius.is_sticking = True

            if 'mwon' in line.lower():
                parsed = True
                arrhenius.motz_wise = True

            if 'mwoff' in line.lower():
                parsed = True
                arrhenius.motz_wise = False

            if 'dup' in line.lower():
                # Duplicate reaction
                parsed = True
                reaction.duplicate = True

            if 'low' in line.lower():
                # Low-pressure-limit Arrhenius parameters for "falloff" reaction
                parsed = True
                tokens = tokens[1].split()
                arrheniusLow = Arrhenius(
                    A=(float(tokens[0].strip()), klow_units),
                    b=float(tokens[1].strip()),
                    Ea=(float(tokens[2].strip()), energy_units),
                    T0=(1,"K"),
                    parser=self
                )

            elif 'high' in line.lower():
                # High-pressure-limit Arrhenius parameters for "chemically
                # activated" reaction
                parsed = True
                tokens = tokens[1].split()
                arrheniusHigh = Arrhenius(
                    A=(float(tokens[0].strip()), kunits),
                    b=float(tokens[1].strip()),
                    Ea=(float(tokens[2].strip()), energy_units),
                    T0=(1,"K"),
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
                                           thirdBody=reaction.thirdBody,
                                           reversible=False,
                                           parser=self)

                    revReaction.kinetics = reaction_type(
                        A=(float(tokens[0].strip()), klow_units),
                        b=float(tokens[1].strip()),
                        Ea=(float(tokens[2].strip()), energy_units),
                        T0=(1,"K"),
                        parser=self
                    )
                    if thirdBody:
                        revReaction.kinetics = ThirdBody(
                            arrheniusHigh=revReaction.kinetics,
                            parser=self)

            elif 'ford' in line.lower():
                parsed = True
                tokens = tokens[1].split()
                reaction.fwdOrders[tokens[0].strip()] = tokens[1].strip()

            elif 'troe' in line.lower():
                # Troe falloff parameters
                parsed = True
                tokens = tokens[1].split()
                alpha = float(tokens[0].strip())
                T3 = float(tokens[1].strip())
                T1 = float(tokens[2].strip())
                T2 = float(tokens[3].strip()) if len(tokens) > 3 else None

                falloff = Troe(
                    alpha=(alpha,''),
                    T3=(T3,"K"),
                    T1=(T1,"K"),
                    T2=(T2,"K") if T2 is not None else None,
                )
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
                arrhenius.coverages.append(
                    [C[0], fortFloat(C[1]), fortFloat(C[2]), fortFloat(C[3])])

            elif 'cheb' in line.lower():
                # Chebyshev parameters
                parsed = True
                if chebyshev is None:
                    chebyshev = Chebyshev()
                tokens = [t.strip() for t in tokens]
                if contains(tokens, 'TCHEB'):
                    index = get_index(tokens, 'TCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Tmin = float(tokens2[0].strip())
                    chebyshev.Tmax = float(tokens2[1].strip())
                if contains(tokens, 'PCHEB'):
                    index = get_index(tokens, 'PCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Pmin = (float(tokens2[0].strip()), 'atm')
                    chebyshev.Pmax = (float(tokens2[1].strip()), 'atm')
                if contains(tokens, 'TCHEB') or contains(tokens, 'PCHEB'):
                    pass
                elif chebyshev.degreeT == 0 or chebyshev.degreeP == 0:
                    tokens2 = tokens[1].split()
                    chebyshev.degreeT = int(float(tokens2[0].strip()))
                    chebyshev.degreeP = int(float(tokens2[1].strip()))
                    chebyshev.coeffs = np.zeros((chebyshev.degreeT, chebyshev.degreeP), np.float64)
                    chebyshevCoeffs.extend([float(t.strip()) for t in tokens2[2:]])
                else:
                    tokens2 = tokens[1].split()
                    chebyshevCoeffs.extend([float(t.strip()) for t in tokens2])

            elif 'plog' in line.lower():
                # Pressure-dependent Arrhenius parameters
                parsed = True
                if pdepArrhenius is None:
                    pdepArrhenius = []
                tokens = tokens[1].split()
                pdepArrhenius.append([float(tokens[0].strip()), Arrhenius(
                    A=(float(tokens[1].strip()), kunits),
                    b=float(tokens[2].strip()),
                    Ea=(float(tokens[3].strip()), energy_units),
                    T0=(1,"K"),
                    parser=self
                )])
            elif len(tokens) >= 2:
                # Assume a list of collider efficiencies
                parsed = True
                for collider, efficiency in zip(tokens[0::2], tokens[1::2]):
                    efficiencies[collider.strip()] = float(efficiency.strip())

            if not parsed:
                raise InputParseError('Unparsable line:\n"""\n{0}\n"""'.format(line))


        # Decide which kinetics to keep and store them on the reaction object
        # Only one of these should be true at a time!
        if chebyshev is not None:
            if chebyshev.Tmin is None or chebyshev.Tmax is None:
                raise InputParseError('Missing TCHEB line for reaction {0}'.format(reaction))
            if chebyshev.Pmin is None or chebyshev.Pmax is None:
                raise InputParseError('Missing PCHEB line for reaction {0}'.format(reaction))
            index = 0
            for t in range(chebyshev.degreeT):
                for p in range(chebyshev.degreeP):
                    chebyshev.coeffs[t,p] = chebyshevCoeffs[index]
                    index += 1
            if quantity_units == 'mol' and self.quantity_units == 'molec':
                chebyshev.coeffs[0, 0] -= np.log10(Avogadro)*quantity_dim
            elif quantity_units == 'molec' and self.quantity_units == 'mol':
                chebyshev.coeffs[0, 0] += np.log10(Avogadro)*quantity_dim

            reaction.kinetics = chebyshev
        elif pdepArrhenius is not None:
            reaction.kinetics = PDepArrhenius(
                pressures=([P for P, arrh in pdepArrhenius], "atm"),
                arrhenius=[arrh for P, arrh in pdepArrhenius],
                parser=self
            )
        elif arrheniusLow is not None:
            reaction.kinetics = Falloff(arrheniusHigh=arrhenius,
                                        arrheniusLow=arrheniusLow,
                                        F=falloff,
                                        parser=self,
                                        efficiencies=efficiencies)
        elif arrheniusHigh is not None:
            reaction.kinetics = ChemicallyActivated(arrheniusHigh=arrheniusHigh,
                                                    arrheniusLow=arrhenius,
                                                    F=falloff,
                                                    parser=self,
                                                    efficiencies=efficiencies)
        elif thirdBody:
            reaction.kinetics = ThirdBody(arrheniusHigh=arrhenius,
                                          parser=self,
                                          efficiencies=efficiencies)
        elif reaction.thirdBody:
            raise InputParseError('Reaction equation implies pressure '
                'dependence but no alternate rate parameters (i.e. HIGH or '
                'LOW) were given for reaction {0}'.format(reaction))
        else:
            reaction.kinetics = arrhenius

        if revReaction:
            revReaction.duplicate = reaction.duplicate
            revReaction.kinetics.efficiencies = reaction.kinetics.efficiencies

        return reaction, revReaction

    def loadChemkinFile(self, path, skipUndeclaredSpecies=True, surface=False):
        """
        Load a Chemkin-format input file to `path` on disk.
        """

        transportLines = []
        self.line_number = 0

        with open(path, 'rU') as ck_file:

            def readline():
                self.line_number += 1
                line = strip_nonascii(ck_file.readline())
                if '!' in line:
                    return line.split('!', 1)
                elif line:
                    return line, ''
                else:
                    return None, None

            line, comment = readline()
            advance = True
            inHeader = True
            while line is not None:
                tokens = line.split() or ['']
                if inHeader and not line.strip():
                    self.headerLines.append(comment.rstrip())

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
                        self.addElement(token)

                elif tokens[0].upper().startswith('SPEC'):
                    # List of species identifiers
                    tokens = tokens[1:]
                    inHeader = False
                    while line is not None and get_index(line, 'END') is None:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('REAC', 'REACTIONS', 'TRAN',
                                                  'TRANSPORT', 'THER', 'THERMO'):
                            self.warn('"SPECIES" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
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
                        if token in self.speciesDict:
                            species = self.speciesDict[token]
                            self.warn('Found additional declaration of species {0}'.format(species))
                        else:
                            species = Species(label=token)
                            self.speciesDict[token] = species
                            self.speciesList.append(species)

                elif tokens[0].upper().startswith('SITE'):
                    # List of species identifers for surface species
                    if '/' in tokens[0]:
                        surfname = tokens[0].split('/')[1]
                    else:
                        surfname = 'surface{}'.format(len(self.surfaces)+1)
                    tokens = tokens[1:]
                    siteDensity = None
                    for token in tokens[:]:
                        if token.upper().startswith('SDEN/'):
                            siteDensity = fortFloat(token.split('/')[1])
                            tokens.remove(token)

                    if siteDensity is None:
                        raise InputParseError('SITE section defined with no site density')
                    self.surfaces.append(Surface(name=surfname,
                                                 density=siteDensity))
                    surf = self.surfaces[-1]

                    inHeader = False
                    while line is not None and get_index(line, 'END') is None:
                        # Grudging support for implicit end of section
                        start = line.strip().upper().split()
                        if start and start[0] in ('REAC', 'REACTIONS', 'THER',
                                                  'THERMO'):
                            self.warn('"SITE" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
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
                        if token in self.speciesDict:
                            species = self.speciesDict[token]
                            self.warn('Found additional declaration of species {0}'.format(species))
                        else:
                            species = Species(label=token, sites=sites)
                            self.speciesDict[token] = species
                            surf.speciesList.append(species)

                elif tokens[0].upper().startswith('THER') and contains(line, 'NASA9'):
                    inHeader = False
                    entryLength = None
                    entry = []
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
                            label, thermo, comp, note = self.readNasa9Entry(entry)
                            entry = []
                            if label not in self.speciesDict:
                                if skipUndeclaredSpecies:
                                    logging.info('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
                                    continue
                                else:
                                    # Add a new species entry
                                    species = Species(label=label)
                                    self.speciesDict[label] = species
                                    self.speciesList.append(species)
                            else:
                                species = self.speciesDict[label]

                            # use the first set of thermo data found
                            if species.thermo is not None:
                                self.warn('Found additional thermo entry for species {0}. '
                                          'If --permissive was given, the first entry is used.'.format(label))
                            else:
                                species.thermo = thermo
                                species.composition = comp
                                species.note = note

                elif tokens[0].upper().startswith('THER'):
                    # List of thermodynamics (hopefully one per species!)
                    inHeader = False
                    line, comment = readline()
                    if line is not None and get_index(line, 'END') is None:
                        TintDefault = float(line.split()[1])
                    thermo = []
                    current = []
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
                                    label, thermo, comp, note = self.readThermoEntry(thermo, TintDefault)
                                except Exception as e:
                                    error_line_number = self.line_number - len(current) + 1
                                    error_entry = ''.join(current).rstrip()
                                    logging.info(
                                        'Error while reading thermo entry starting on line {0}:\n'
                                        '"""\n{1}\n"""'.format(error_line_number, error_entry)
                                    )
                                    raise

                                if label not in self.speciesDict:
                                    if skipUndeclaredSpecies:
                                        logging.info('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
                                        thermo = []
                                        line, comment = readline()
                                        current = []
                                        continue
                                    else:
                                        # Add a new species entry
                                        species = Species(label=label)
                                        self.speciesDict[label] = species
                                        self.speciesList.append(species)
                                else:
                                    species = self.speciesDict[label]

                                # use the first set of thermo data found
                                if species.thermo is not None:
                                    self.warn('Found additional thermo entry for species {0}. '
                                              'If --permissive was given, the first entry is used.'.format(label))
                                else:
                                    species.thermo = thermo
                                    species.composition = comp
                                    species.note = note

                                thermo = []
                                current = []
                        elif thermo and thermo[-1].rstrip().endswith('&'):
                            # Include Chemkin-style extended elemental composition
                            thermo.append(line)
                        line, comment = readline()

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
                            raise InputParseError("Unrecognized token on REACTIONS line, {0!r}".format(token))

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
                    # correctly, we associate them with the next reaction (and
                    # keep track of the final trailing comment separately)
                    if kineticsList and kineticsList[0] == '':
                        kineticsList.pop(0)
                        self.finalReactionComment = commentsList.pop()

                    self.setupKinetics()
                    for kinetics, comment, line_number in zip(kineticsList, commentsList, startLines):
                        try:
                            reaction, revReaction = self.readKineticsEntry(kinetics, surface)
                        except Exception as e:
                            self.line_number = line_number
                            logging.info('Error reading reaction starting on '
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
                    raise InputParseError('Section starts with unrecognized keyword.'
                        '\n"""\n{0}\n"""'.format(line.rstrip()))

                if advance:
                    line, comment = readline()
                else:
                    advance = True

        self.checkDuplicateReactions()

        index = 0
        for reaction in self.reactions:
            index += 1
            reaction.index = index

        if transportLines:
            self.parseTransportData(transportLines, path, transport_start_line)

    def checkDuplicateReactions(self):
        """
        Check for marked (and unmarked!) duplicate reactions. Raise exception
        for unmarked duplicate reactions.

        Pressure-independent and pressure-dependent reactions are treated as
        different, so they don't need to be marked as duplicate.
        """
        message = ('Encountered unmarked duplicate reaction {0} '
                   '(See lines {1} and {2} of the input file.).')

        possible_duplicates = defaultdict(list)
        for r in self.reactions:
            k = (tuple(r.reactants), tuple(r.products), r.kinetics.isPressureDependent())
            possible_duplicates[k].append(r)

        for reactions in possible_duplicates.values():
            for r1,r2 in itertools.combinations(reactions, 2):
                if r1.duplicate and r2.duplicate:
                    pass  # marked duplicate reaction
                elif (r1.thirdBody and r1.thirdBody.upper() == 'M' and
                      r1.kinetics.efficiencies.get(r2.thirdBody) == 0):
                    pass  # explicit zero efficiency
                elif (r2.thirdBody and r2.thirdBody.upper() == 'M' and
                      r2.kinetics.efficiencies.get(r1.thirdBody) == 0):
                    pass  # explicit zero efficiency
                elif r1.thirdBody != r2.thirdBody:
                    pass  # distinct third bodies
                else:
                    raise InputParseError(message.format(r1, r1.line_number, r2.line_number))

    def parseTransportData(self, lines, filename, line_offset):
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
                comment = None

            data = line.split()

            speciesName = data[0]
            if speciesName in self.speciesDict:
                if len(data) != 7:
                    raise InputParseError('Unable to parse line {0} of {1}:\n"""\n{2}"""\n'
                        '6 transport parameters expected, but found {3}.'.format(
                            line_offset + i, filename, original_line, len(data)-1))

                if self.speciesDict[speciesName].transport is None:
                    self.speciesDict[speciesName].transport = TransportData(*data, comment=comment)
                else:
                    self.warn('Ignoring duplicate transport data'
                         ' for species "{0}" on line {1} of "{2}".'.format(
                            speciesName, line_offset + i, filename))

    def getSpeciesString(self, speciesList, indent):
        speciesNameLength = 1
        allElements = set(self.elements)
        for s in speciesList:
            if s.composition is None:
                raise InputParseError('No thermo data found for species: {0!r}'.format(s.label))
            missingElements = set(s.composition) - allElements
            if missingElements:
                raise InputParseError("Undefined elements in species '{}':"
                    " {}".format(s.label, ','.join(repr(e) for e in missingElements)))
            speciesNameLength = max(speciesNameLength, len(s.label))

        speciesNames = ['']
        speciesPerLine = max(int((80-indent)/(speciesNameLength + 2)), 1)

        for i,s in enumerate(speciesList):
            if i and not i % speciesPerLine:
                speciesNames.append(' '*indent)
            speciesNames[-1] += '{0:{1}s}'.format(s.label, speciesNameLength+2)

        speciesNames = '\n'.join(line.rstrip() for line in speciesNames)

        return speciesNames

    def writeCTI(self, header=None, name='gas', transportModel='Mix',
                 outName='mech.cti'):

        delimiterLine = '#' + '-'*79
        haveTransport = True
        for s in self.speciesList:
            if not s.transport:
                haveTransport = False

        lines = []
        surface_names = []

        # Assign IDs to reactions if necessary
        nReactingPhases = 0
        if self.reactions:
            nReactingPhases += 1
        for surf in self.surfaces:
            surface_names.append(surf.name)
            if surf.reactions:
                nReactingPhases += 1

        use_reaction_ids = nReactingPhases > 1
        if use_reaction_ids:
            for i,R in enumerate(self.reactions):
                R.ID = '{0}-{1}'.format(name, i+1)
            for surf in self.surfaces:
                for i,R in enumerate(surf.reactions):
                    R.ID = '{0}-{1}'.format(surf.name, i+1)

        # Original header
        if self.headerLines:
            lines.append('"""')
            lines.extend(self.headerLines)
            lines.extend(('"""', ''))

        # Cantera-generated header
        if header:
            lines.extend(header)

        if name is not None:
            speciesNames = self.getSpeciesString(self.speciesList, 21)
            # Write the gas definition
            lines.append("units(length='cm', time='s', quantity={0!r}, act_energy={1!r})".format(self.output_quantity_units, self.output_energy_units))
            lines.append('')
            lines.append('ideal_gas(name={0!r},'.format(name))
            lines.append('          elements="{0}",'.format(' '.join(self.elements)))
            lines.append('          species="""{0}""",'.format(speciesNames))
            if self.reactions:
                if not use_reaction_ids:
                    lines.append("          reactions='all',")
                else:
                    lines.append("          reactions='{0}-*',".format(name))
            if haveTransport:
                lines.append("          transport={0!r},".format(transportModel))
            lines.append('          initial_state=state(temperature=300.0, pressure=OneAtm))')
            lines.append('')

        for surf in self.surfaces:
            # Write definitions for surface phases
            speciesNames = self.getSpeciesString(surf.speciesList, 26)
            lines.append('ideal_interface(name={0!r},'.format(surf.name))
            lines.append('                elements="{0}",'.format(' '.join(self.elements)))
            lines.append('                species="""{0}""",'.format(speciesNames))
            lines.append('                site_density={0},'.format(surf.siteDensity))
            lines.append('                phases="{0}",'.format(name))
            if surf.reactions:
                if not use_reaction_ids:
                    lines.append("          reactions='all',")
                else:
                    lines.append("          reactions='{0}-*',".format(surf.name))
            lines.append('                initial_state=state(temperature=300.0, pressure=OneAtm))')
            lines.append('')

        # Write data on custom elements
        if self.element_weights:
            lines.append(delimiterLine)
            lines.append('# Element data')
            lines.append(delimiterLine)
            lines.append('')
            for name, weight in sorted(self.element_weights.items()):
                lines.append('element(symbol={0!r}, atomic_mass={1})'.format(name, weight))

        # Write the individual species data
        lines.append(delimiterLine)
        lines.append('# Species data')
        lines.append(delimiterLine)
        lines.append('')

        for s in self.speciesList:
            lines.append(s.to_cti())
        for surf in self.surfaces:
            for s in surf.speciesList:
                lines.append(s.to_cti())

        if self.reactions or any(surf.reactions for surf in self.surfaces):
            # Write the reactions
            lines.append(delimiterLine)
            lines.append('# Reaction data')
            lines.append(delimiterLine)

            if self.motz_wise is True:
                lines.append('enable_motz_wise()')
            elif self.motz_wise is False:
                lines.append('disable_motz_wise()')

            for i,r in enumerate(self.reactions):
                lines.extend('# '+c for c in r.comment.split('\n') if c)
                lines.append('\n# Reaction {0}'.format(i+1))
                lines.append(r.to_cti())

            for surf in self.surfaces:
                for i,r in enumerate(surf.reactions):
                    lines.extend('# '+c for c in r.comment.split('\n') if c)
                    lines.append('\n# {0} Reaction {1}'.format(surf.name, i+1))
                    lines.append(r.to_cti())

            # Comment after the last reaction
            lines.extend('# '+c for c in self.finalReactionComment.split('\n') if c)

            lines.append('')

        with open(outName, 'w') as f:
            f.write('\n'.join(lines))

        return surface_names

    @staticmethod
    def showHelp():
        print("""
ck2cti.py: Convert Chemkin-format mechanisms to Cantera input files (.cti)

Usage:
    ck2cti [--input=<filename>]
           [--thermo=<filename>]
           [--transport=<filename>]
           [--surface=<filename>]
           [--id=<phase-id>]
           [--output=<filename>]
           [--permissive]
           [-d | --debug]

Example:
    ck2cti --input=chem.inp --thermo=therm.dat --transport=tran.dat

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.cti'.

An input file containing only species definitions (which can be referenced from
phase definitions in other input files) can be created by specifying only a
thermo file.

For the case of a surface mechanism, the gas phase input file should be
specified as 'input' and the surface phase input file should be specified as
'surface'.

The '--permissive' option allows certain recoverable parsing errors (e.g.
duplicate transport data) to be ignored.

""")

    @staticmethod
    def convertMech(inputFile, thermoFile=None, transportFile=None,
                    surfaceFile=None, phaseName='gas', outName=None,
                    quiet=False, permissive=None):

        parser = Parser()
        if inputFile:
            inputFile = os.path.expanduser(inputFile)
        if thermoFile:
            thermoFile = os.path.expanduser(thermoFile)
        if transportFile:
            transportFile = os.path.expanduser(transportFile)
        if surfaceFile:
            surfaceFile = os.path.expanduser(surfaceFile)
        if outName:
            outName = os.path.expanduser(outName)

        if quiet:
            logging.basicConfig(level=logging.ERROR)
        else:
            logging.basicConfig(level=logging.INFO)

        if permissive is not None:
            parser.warning_as_error = not permissive

        if inputFile:
            if not os.path.exists(inputFile):
                raise IOError('Missing input file: {0!r}'.format(inputFile))
            try:
                # Read input mechanism files
                parser.loadChemkinFile(inputFile)
            except Exception as err:
                logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n{2}\n".format(
                                inputFile, parser.line_number, err))
                raise
        else:
            phaseName = None

        if surfaceFile:
            if not os.path.exists(surfaceFile):
                raise IOError('Missing input file: {0!r}'.format(surfaceFile))
            try:
                # Read input mechanism files
                parser.loadChemkinFile(surfaceFile, surface=True)
            except Exception as err:
                logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n{2}\n".format(
                                surfaceFile, parser.line_number, err))
                raise

        if thermoFile:
            if not os.path.exists(thermoFile):
                raise IOError('Missing thermo file: {0!r}'.format(thermoFile))
            try:
                parser.loadChemkinFile(thermoFile,
                                     skipUndeclaredSpecies=bool(inputFile))
            except Exception:
                logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n".format(
                                thermoFile, parser.line_number))
                raise

        if transportFile:
            if not os.path.exists(transportFile):
                raise IOError('Missing transport file: {0!r}'.format(transportFile))
            with open(transportFile, 'rU') as f:
                lines = [strip_nonascii(line) for line in f]
            parser.parseTransportData(lines, transportFile, 1)

            # Transport validation: make sure all species have transport data
            for s in parser.speciesList:
                if s.transport is None:
                    raise InputParseError("No transport data for species '{0}'.".format(s))

        if not outName:
            outName = os.path.splitext(inputFile)[0] + '.cti'

        # Write output file
        surface_names = parser.writeCTI(name=phaseName, outName=outName)
        if not quiet:
            nReactions = len(parser.reactions) + sum(len(surf.reactions) for surf in parser.surfaces)
            print('Wrote CTI mechanism file to {0!r}.'.format(outName))
            print('Mechanism contains {0} species and {1} reactions.'.format(len(parser.speciesList), nReactions))
        return surface_names


def convertMech(inputFile, thermoFile=None, transportFile=None, surfaceFile=None,
                phaseName='gas', outName=None, quiet=False, permissive=None):
    return Parser.convertMech(inputFile, thermoFile, transportFile, surfaceFile,
                              phaseName, outName, quiet, permissive)

def main(argv):

    longOptions = ['input=', 'thermo=', 'transport=', 'surface=', 'id=',
                   'output=', 'permissive', 'help', 'debug']

    try:
        optlist, args = getopt.getopt(argv, 'dh', longOptions)
        options = dict()
        for o,a in optlist:
            options[o] = a

        if args:
            raise getopt.GetoptError('Unexpected command line option: ' +
                                     repr(' '.join(args)))

    except getopt.GetoptError as e:
        print('ck2cti.py: Error parsing arguments:')
        print(e)
        print('Run "ck2cti.py --help" to see usage help.')
        sys.exit(1)

    if not options or '-h' in options or '--help' in options:
        Parser.showHelp()
        sys.exit(0)

    if '--input' in options:
        inputFile = options['--input']
    else:
        inputFile = None

    thermoFile = options.get('--thermo')

    if '--output' in options:
        outName = options['--output']
        if not outName.endswith('.cti'):
            outName += '.cti'
    elif inputFile:
        outName = os.path.splitext(inputFile)[0] + '.cti'
    else:
        outName = os.path.splitext(thermoFile)[0] + '.cti'

    permissive = '--permissive' in options
    transportFile = options.get('--transport')
    surfaceFile = options.get('--surface')
    phaseName = options.get('--id', 'gas')

    surfaces = Parser.convertMech(inputFile, thermoFile, transportFile,
                                  surfaceFile, phaseName, outName,
                                  permissive=permissive)

    # Do full validation by importing the resulting mechanism
    if not inputFile:
        # Can't validate input file that don't define a phase
        return

    try:
        import cantera as ct
    except ImportError:
        print('WARNING: Unable to import Cantera Python module. Output '
              'mechanism has not been validated')
        sys.exit(0)

    try:
        print('Validating mechanism...', end='')
        gas = ct.Solution(outName)
        for surfname in surfaces:
            phase = ct.Interface(outName, surfname, [gas])
        print('PASSED.')
    except RuntimeError as e:
        print('FAILED.')
        print(e)
        sys.exit(1)


def script_entry_point():
    main(sys.argv[1:])

if __name__ == '__main__':
    main(sys.argv[1:])
