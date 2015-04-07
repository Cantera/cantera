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
"""

from __future__ import print_function

from collections import defaultdict
import logging
import os.path
import numpy as np
import re
import itertools

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


def compatible_quantities(quantity_basis, units):
    if quantity_basis == 'mol':
        return 'molec' not in units
    elif quantity_basis == 'molec':
        return 'molec' in units or 'mol' not in units
    else:
        raise Exception('Unknown quantity basis: "{0}"'.format(quantity_basis))


class InputParseError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin-format
    mechanism files. Pass a string describing the circumstances that caused
    the exceptional behavior.
    """
    pass


class Species(object):
    def __init__(self, label):
        self.label = label
        self.thermo = None
        self.transport = None
        self.note = None
        self.composition = None

    def __str__(self):
        return self.label

    def to_cti(self, indent=0):
        lines = []
        atoms = ' '.join('{0}:{1}'.format(*a)
                         for a in self.composition.items())

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
        vals = ['{0: 15.8E}'.format(i) for i in self.coeffs]
        if len(self.coeffs) == 7:
            lines = ['NASA([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                     prefix+'     [{0}, {1}, {2},'.format(*vals[0:3]),
                     prefix+'      {0}, {1}, {2},'.format(*vals[3:6]),
                     prefix+'      {0}]),'.format(vals[6])]
        else:
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

    def __init__(self, index=-1, reactants=None, products=None, kinetics=None,
                 reversible=True, duplicate=False, fwdOrders=None,
                 thirdBody=None):
        self.index = index
        self.reactants = reactants  # list of (stoichiometry, species) tuples
        self.products = products  # list of (stoichiometry, specis) tuples
        self.kinetics = kinetics
        self.reversible = reversible
        self.duplicate = duplicate
        self.fwdOrders = fwdOrders if fwdOrders is not None else {}
        self.thirdBody = thirdBody

    def _coeff_string(self, coeffs):
        L = []
        for stoichiometry,species in coeffs:
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

        if self.duplicate:
            kinstr = kinstr[:-1] + ",\n{0}options='duplicate')".format(k_indent)

        if self.fwdOrders:
            order = ' '.join('{0}:{1}'.format(k,v)
                             for (k,v) in self.fwdOrders.items())
            kinstr = kinstr[:-1] + ",\n{0}order='{1}')".format(k_indent, order)

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

    def efficiencyString(self):
        return ' '.join('{0}:{1}'.format(mol, eff)
                        for mol,eff in self.efficiencies.items())


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

    .. math:: k(T) = A \\left( \\frac{T}{T_0} \\right)^n \\exp \\left( - \\frac{E_\\mathrm{a}}{RT} \\right)

    where :math:`A`, :math:`n`, :math:`E_\\mathrm{a}`, and :math:`T_0` are the
    parameters to be set, :math:`T` is absolute temperature, and :math:`R` is
    the gas law constant. The attributes are:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `A`             :class:`Quantity`   The preexponential factor in s^-1, m^3/mol*s, etc.
    `T0`            :class:`Quantity`   The reference temperature in K
    `n`             :class:`Quantity`   The temperature exponent
    `Ea`            :class:`Quantity`   The activation energy in J/mol
    =============== =================== ========================================

    """

    def __init__(self, A=0.0, n=0.0, Ea=0.0, T0=1.0, **kwargs):
        KineticsModel.__init__(self, **kwargs)
        self.A = A
        self.T0 = T0
        self.n = n
        self.Ea = Ea

    def isPressureDependent(self):
        """
        Returns ``False`` since Arrhenius kinetics are not pressure-dependent.
        """
        return False

    def rateStr(self):
        if compatible_quantities(self.parser.quantity_units, self.A[1]):
            A = '{0:e}'.format(self.A[0])
        else:
            A = "({0:e}, '{1}')".format(*self.A)

        if self.Ea[1] == self.parser.energy_units:
            Ea = str(self.Ea[0])
        else:
            Ea = "({0}, '{1}')".format(*self.Ea)

        return '[{0}, {1}, {2}]'.format(A, self.n, Ea)

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstring = reactantstr + arrow + productstr
        return 'reaction({0!r}, {1})'.format(rxnstring, self.rateStr())


class PDepArrhenius(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    expression

    .. math:: k(T,P) = A(P) T^{n(P)} \\exp \\left[ \\frac{-E_\\mathrm{a}(P)}{RT} \\right]

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
        template = '[({0}, {1!r}), {2.A[0]:e}, {2.n}, {2.Ea[0]}],'
        for pressure,arrhenius in zip(self.pressures[0], self.arrhenius):
            lines.append(prefix + template.format(pressure,
                                                  self.pressures[1],
                                                  arrhenius))
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

    def __init__(self, coeffs=None, kunits='', **kwargs):
        KineticsModel.__init__(self, **kwargs)
        if coeffs is not None:
            self.coeffs = np.array(coeffs, np.float64)
            self.degreeT = self.coeffs.shape[0]
            self.degreeP = self.coeffs.shape[1]
        else:
            self.coeffs = None
            self.degreeT = 0
            self.degreeP = 0
        self.kunits = kunits

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
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    "SRI" formulation of the blending function :math:`F` using either 3 or
    5 parameters. See :ref:`sec-sri-falloff`.

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

        if int(geometry) not in (0,1,2):
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

def isnumberlike(text):
    """ Returns true if `text` can be interpreted as a floating point number. """
    try:
        float(text)
        return True
    except ValueError:
        return False

def get_index(seq, value):
    """
    Find the first location in *seq* which contains a case-insensitive,
    whitespace-insensitive match for *value*. Returns *None* if no match is
    found.
    """
    if isinstance(seq, str):
        seq = seq.split()
    value = value.lower().strip()
    for i,item in enumerate(seq):
        if item.lower() == value:
            return i
    return None

def contains(seq, value):
    if isinstance(seq, str):
        return value.lower() in seq.lower()
    else:
        return get_index(seq, value) is not None


class Parser(object):
    def __init__(self):
        self.processed_units = False
        self.energy_units = 'cal/mol'
        self.quantity_units = 'mol'
        self.warning_as_error = True

        self.elements = []
        self.speciesList = []
        self.speciesDict = {}
        self.reactions = []

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

        # Extract the NASA polynomial coefficients
        # Remember that the high-T polynomial comes first!
        try:
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

        except (IndexError, ValueError) as err:
            raise InputParseError('Error while reading thermo entry for species {0}:\n{1}'.format(species, err))

        composition = self.parseComposition(lines[0][24:44], 4, 5)

        # Non-standard extended elemental composition data may be located beyond
        # column 80 on the first line of the thermo entry
        if len(lines[0]) > 80:
            elements = lines[0][80:]
            composition2 = self.parseComposition(elements, len(elements)//10, 10)
            composition.update(composition2)

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
        for next_char in ('<','=','(','+','\n'):
            self.species_tokens.update(k + next_char for k in self.speciesDict)
        self.other_tokens = {'M': 'third-body', 'm': 'third-body',
                             '(+M)': 'falloff3b', '(+m)': 'falloff3b',
                             '<=>': 'equal', '=>': 'equal', '=': 'equal'}
        self.other_tokens.update(('(+%s)' % k, 'falloff3b: %s' % k) for k in self.speciesDict)
        self.Slen = max(map(len, self.other_tokens))

    def readKineticsEntry(self, entry):
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
        n = float(tokens[-2])
        Ea = float(tokens[-1])
        reaction = ''.join(tokens[:-3]) + '\n'

        # Identify species tokens in the reaction expression in order of
        # decreasing length
        locs = {}
        for i in range(self.Slen, 0, -1):
            for j in range(len(reaction)-i+1):
                test = reaction[j:j+i]
                if test in self.species_tokens:
                    reaction = reaction[:j] + ' '*(i-1) + reaction[j+i-1:]
                    locs[j] = test[:-1], 'species'

        # Identify other tokens in the reaction expression in order of
        # descending length
        for i in range(self.Slen, 0, -1):
            for j in range(len(reaction)-i+1):
                test = reaction[j:j+i]
                if test in self.other_tokens:
                    reaction = reaction[:j] + ' '*i + reaction[j+i:]
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
                    raise InputParseError('Unexpected token "{0}" in reaction expression "{1}".'.format(token, reaction))

        reactants = []
        products = []
        stoichiometry = 1
        lhs = True
        for token,kind in [v for k,v in sorted(locs.items())]:
            if kind == 'equal':
                reversible = token in ('<=>', '=')
                lhs = False
            elif kind == 'coeff':
                stoichiometry = token
            elif lhs:
                reactants.append((stoichiometry,token,kind))
                stoichiometry = 1
            else:
                products.append((stoichiometry,token,kind))
                stoichiometry = 1

        if lhs is True:
            raise InputParseError("Failed to find reactant/product delimiter in reaction string.")

        # Create a new Reaction object for this reaction
        reaction = Reaction(reactants=[], products=[], reversible=reversible)

        def parseExpression(expression, dest):
            falloff3b = None
            thirdBody = False  # simple third body reaction (non-falloff)
            for stoichiometry,species,kind in expression:
                if kind == 'third-body':
                    thirdBody = True
                elif kind == 'falloff3b':
                    falloff3b = 'M'
                elif kind.startswith('falloff3b:'):
                    falloff3b = kind.split()[1]
                else:
                    dest.append((stoichiometry, self.speciesDict[species]))

            return falloff3b, thirdBody

        falloff_3b_r, thirdBody = parseExpression(reactants, reaction.reactants)
        falloff_3b_p, thirdBody = parseExpression(products, reaction.products)

        if falloff_3b_r != falloff_3b_p:
            raise InputParseError('Third bodies do not match: "{0}" and "{1}" in'
                ' reaction entry:\n\n{2}'.format(falloff_3b_r, falloff_3b_p, entry))

        reaction.thirdBody = falloff_3b_r

        # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
        # This assumes elementary kinetics for all reactions
        rStoich = sum(r[0] for r in reaction.reactants) + (1 if thirdBody else 0)
        if rStoich > 3 or rStoich < 1:
            raise InputParseError('Invalid number of reactant species ({0}) for reaction {1}.'.format(rStoich, reaction))

        length_dim = 3 * (rStoich - 1)
        quantity_dim = rStoich - 1
        kunits = self.getRateConstantUnits(length_dim, 'cm',
                                           quantity_dim, quantity_units)
        klow_units = self.getRateConstantUnits(length_dim + 3, 'cm',
                                               quantity_dim + 1, quantity_units)

        # The rest of the first line contains Arrhenius parameters
        arrhenius = Arrhenius(
            A=(A,kunits),
            n=n,
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
            tokens = line.split('/')
            if 'dup' in line.lower():
                # Duplicate reaction
                reaction.duplicate = True

            elif 'low' in line.lower():
                # Low-pressure-limit Arrhenius parameters for "falloff" reaction
                tokens = tokens[1].split()
                arrheniusLow = Arrhenius(
                    A=(float(tokens[0].strip()),klow_units),
                    n=float(tokens[1].strip()),
                    Ea=(float(tokens[2].strip()),energy_units),
                    T0=(1,"K"),
                    parser=self
                )

            elif 'high' in line.lower():
                # High-pressure-limit Arrhenius parameters for "chemically
                # activated" reaction
                tokens = tokens[1].split()
                arrheniusHigh = Arrhenius(
                    A=(float(tokens[0].strip()),kunits),
                    n=float(tokens[1].strip()),
                    Ea=(float(tokens[2].strip()),energy_units),
                    T0=(1,"K"),
                    parser=self
                )
                # Need to fix units on the base reaction:
                arrhenius.A = (arrhenius.A[0], klow_units)

            elif 'rev' in line.lower():
                reaction.reversible = False

                # Create a reaction proceeding in the opposite direction
                revReaction = Reaction(reactants=reaction.products,
                                       products=reaction.reactants,
                                       thirdBody=reaction.thirdBody,
                                       reversible=False)
                tokens = tokens[1].split()
                revReaction.kinetics = Arrhenius(
                    A=(float(tokens[0].strip()),klow_units),
                    n=float(tokens[1].strip()),
                    Ea=(float(tokens[2].strip()),energy_units),
                    T0=(1,"K"),
                    parser=self
                )
                if thirdBody:
                    revReaction.kinetics = ThirdBody(
                        arrheniusHigh=revReaction.kinetics,
                        parser=self)

            elif 'ford' in line.lower():
                tokens = tokens[1].split()
                reaction.fwdOrders[tokens[0].strip()] = tokens[1].strip()

            elif 'troe' in line.lower():
                # Troe falloff parameters
                tokens = tokens[1].split()
                alpha = float(tokens[0].strip())
                T3 = float(tokens[1].strip())
                T1 = float(tokens[2].strip())
                try:
                    T2 = float(tokens[3].strip())
                except (IndexError, ValueError):
                    T2 = None

                falloff = Troe(
                    alpha=(alpha,''),
                    T3=(T3,"K"),
                    T1=(T1,"K"),
                    T2=(T2,"K") if T2 is not None else None,
                )
            elif 'sri' in line.lower():
                # SRI falloff parameters
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

            elif 'cheb' in line.lower():
                # Chebyshev parameters
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
                    chebyshev.coeffs = np.zeros((chebyshev.degreeT,chebyshev.degreeP), np.float64)
                    chebyshevCoeffs.extend([float(t.strip()) for t in tokens2[2:]])
                else:
                    tokens2 = tokens[1].split()
                    chebyshevCoeffs.extend([float(t.strip()) for t in tokens2])

            elif 'plog' in line.lower():
                # Pressure-dependent Arrhenius parameters
                if pdepArrhenius is None:
                    pdepArrhenius = []
                tokens = tokens[1].split()
                pdepArrhenius.append([float(tokens[0].strip()), Arrhenius(
                    A=(float(tokens[1].strip()),kunits),
                    n=float(tokens[2].strip()),
                    Ea=(float(tokens[3].strip()),energy_units),
                    T0=(1,"K"),
                    parser=self
                )])
            else:
                # Assume a list of collider efficiencies
                for collider, efficiency in zip(tokens[0::2], tokens[1::2]):
                    efficiencies[collider.strip()] = float(efficiency.strip())

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
            reaction.kinetics = chebyshev
        elif pdepArrhenius is not None:
            reaction.kinetics = PDepArrhenius(
                pressures=([P for P, arrh in pdepArrhenius],"atm"),
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
        else:
            reaction.kinetics = arrhenius

        if revReaction:
            revReaction.duplicate = reaction.duplicate
            revReaction.kinetics.efficiencies = reaction.kinetics.efficiencies

        return reaction, revReaction

    def loadChemkinFile(self, path):
        """
        Load a Chemkin-format input file to `path` on disk.
        """

        transportLines = []

        with open(path, 'rU') as ck_file:
            self.line_number = 0

            def readline():
                self.line_number += 1
                line = ck_file.readline()
                if '!' in line:
                    return line.split('!', 1)
                elif line:
                    return line, ''
                else:
                    return None, None

            line, comment = readline()
            advance = True
            while line is not None:
                tokens = line.split() or ['']

                if tokens[0].upper().startswith('ELEM'):
                    tokens = tokens[1:]
                    while line is not None and not contains(line, 'END'):
                        # Grudging support for implicit end of section
                        if contains(line, 'SPEC'):
                            self.warn('"ELEMENTS" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            tokens.pop()
                            break

                        line, comment = readline()
                        tokens.extend(line.split())

                    for token in tokens:
                        if token.upper() == 'END':
                            break
                        self.elements.append(token.capitalize())

                elif tokens[0].upper().startswith('SPEC'):
                    # List of species identifiers
                    tokens = tokens[1:]
                    while line is not None and not contains(line, 'END'):
                        # Grudging support for implicit end of section
                        if (contains(line, 'REAC') or contains(line, 'TRAN') or
                            contains(line, 'THER')):
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

                elif tokens[0].upper().startswith('THER') and contains(line, 'NASA9'):
                    entryPosition = 0
                    entryLength = None
                    entry = []
                    while line is not None and not get_index(line, 'END') == 0:
                        # Grudging support for implicit end of section
                        if (contains(line, 'REAC') or contains(line, 'TRAN')):
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

                        if entryPosition == 0:
                            entry.append(line)
                        elif entryPosition == 1:
                            entryLength = 2 + 3 * int(line.split()[0])
                            entry.append(line)
                        elif entryPosition < entryLength:
                            entry.append(line)

                        if entryPosition == entryLength-1:
                            label, thermo, comp, note = self.readNasa9Entry(entry)
                            try:
                                species = self.speciesDict[label]
                                # use the first set of thermo data found
                                if species.thermo is not None:
                                    self.warn('Found additional thermo entry for species {0}'.format(label))
                                else:
                                    species.thermo = thermo
                                    species.composition = comp
                                    species.note = note
                            except KeyError:
                                logging.info('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))

                            entryPosition = -1
                            entry = []

                        entryPosition += 1

                elif tokens[0].upper().startswith('THER'):
                    # List of thermodynamics (hopefully one per species!)
                    line, comment = readline()
                    if line is not None and not contains(line, 'END'):
                        TintDefault = float(line.split()[1])
                    thermo = []
                    while line is not None and not contains(line, 'END'):
                        # Grudging support for implicit end of section
                        if contains(line, 'REAC') or contains(line, 'TRAN'):
                            self.warn('"THERMO" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            tokens.pop()
                            break

                        if len(line) >= 80 and line[79] in ['1', '2', '3', '4']:
                            thermo.append(line)
                            if line[79] == '4':
                                label, thermo, comp, note = self.readThermoEntry(thermo, TintDefault)
                                try:
                                    species = self.speciesDict[label]
                                    # use the first set of thermo data found
                                    if species.thermo is not None:
                                        self.warn('Found additional thermo entry for species {0}'.format(label))
                                    else:
                                        species.thermo = thermo
                                        species.composition = comp
                                        species.note = note
                                except KeyError:
                                    logging.info('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
                                thermo = []
                        line, comment = readline()

                elif tokens[0].upper().startswith('REAC'):
                    # Reactions section

                    for token in tokens[1:]:
                        units = token.upper()
                        if units in ENERGY_UNITS:
                            if (self.processed_units and
                                self.energy_units != ENERGY_UNITS[units]):
                                raise InputParseError("Multiple REACTIONS sections with "
                                                      "different units are not supported.")
                            self.energy_units = ENERGY_UNITS[units]
                        elif units in QUANTITY_UNITS:
                            if (self.processed_units and
                                self.quantity_units != QUANTITY_UNITS[units]):
                                raise InputParseError("Multiple REACTIONS sections with "
                                                      "different units are not supported.")
                            self.quantity_units = QUANTITY_UNITS[units]
                        else:
                            raise InputParseError("Unrecognized energy or quantity unit, {0!r}".format(units))

                    if len(tokens) > 1:
                        self.processed_units = True

                    kineticsList = []
                    commentsList = []
                    startLines = []
                    kinetics = ''
                    comments = ''

                    line, comment = readline()
                    while line is not None and not contains(line, 'END'):
                        # Grudging support for implicit end of section
                        if contains(line, 'TRAN'):
                            self.warn('"REACTIONS" section implicitly ended by start of '
                                      'next section on line {0}.'.format(self.line_number))
                            advance = False
                            break

                        lineStartsWithComment = not line and comment
                        line = line.strip()
                        comment = comment.strip()

                        if '=' in line and not lineStartsWithComment:
                            # Finish previous record
                            kineticsList.append(kinetics)
                            commentsList.append(comments)
                            startLines.append(self.line_number)
                            kinetics = ''
                            comments = ''

                        if line:
                            kinetics += line + '\n'
                        if comment:
                            comments += comment + '\n'

                        line, comment = readline()

                    # Don't forget the last reaction!
                    if kinetics.strip() != '':
                        kineticsList.append(kinetics)
                        commentsList.append(comments)

                    if kineticsList[0] == '' and commentsList[-1] == '':
                        # True for mechanism files generated from RMG-Py
                        kineticsList.pop(0)
                        commentsList.pop(-1)
                    elif kineticsList[0] == '' and commentsList[0] == '':
                        # True for mechanism files generated from RMG-Java
                        kineticsList.pop(0)
                        commentsList.pop(0)
                    else:
                        # In reality, comments can occur anywhere in the mechanism
                        # file (e.g. either or both of before and after the
                        # reaction equation)
                        # If we can't tell what semantics we are using, then just
                        # throw the comments away
                        # (This is better than failing to load the mechanism file at
                        # all, which would likely occur otherwise)
                        if kineticsList[0] == '':
                            kineticsList.pop(0)
                        if len(kineticsList) != len(commentsList):
                            commentsList = ['' for kinetics in kineticsList]

                    self.setupKinetics()
                    for kinetics, comments, line_number in zip(kineticsList, commentsList, startLines):
                        try:
                            reaction,revReaction = self.readKineticsEntry(kinetics)
                        except Exception as e:
                            logging.error('Error reading reaction entry starting on line {0}:'.format(line_number))
                            raise
                        reaction.line_number = line_number
                        self.reactions.append(reaction)
                        if revReaction is not None:
                            revReaction.line_number = line_number
                            self.reactions.append(revReaction)

                elif tokens[0].upper().startswith('TRAN'):
                    line, comment = readline()
                    while line is not None and not contains(line, 'END'):
                        # Grudging support for implicit end of section
                        if contains(line, 'REAC'):
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
            self.parseTransportData(transportLines)

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
                    pass # marked duplicate reaction
                elif (r1.thirdBody and r1.thirdBody.upper() == 'M' and
                      r1.kinetics.efficiencies.get(r2.thirdBody) == 0):
                    pass # explicit zero efficiency
                elif (r2.thirdBody and r2.thirdBody.upper() == 'M' and
                      r2.kinetics.efficiencies.get(r1.thirdBody) == 0):
                    pass # explicit zero efficiency
                elif r1.thirdBody != r2.thirdBody:
                    pass # distinct third bodies
                else:
                    raise InputParseError(message.format(r1, r1.line_number, r2.line_number))

    def parseTransportData(self, lines):
        """
        Parse the Chemkin-format transport data in ``lines`` (a list of strings)
        and add that transport data to the previously-loaded species.
        """

        for line in lines:
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            if get_index(line, 'END') == 0:
                break

            if '!' in line:
                line, comment = line.split('!', 1)
                data = line.split() + [comment]
            else:
                data = line.split()
            if len(data) < 7:
                raise InputParseError('Unable to parse transport data: not enough parameters')

            speciesName = data[0]
            if speciesName in self.speciesDict:
                if self.speciesDict[speciesName].transport is None:
                    self.speciesDict[speciesName].transport = TransportData(*data)
                else:
                    self.warn('Ignoring duplicate transport data'
                         ' for species "{0}".'.format(speciesName))

    def writeCTI(self, header=None, name='gas', transportModel='Mix',
                 outName='mech.cti'):

        delimiterLine = '#' + '-'*79
        haveTransport = True
        speciesNameLength = 1
        elementsFromSpecies = set()
        for s in self.speciesList:
            if not s.transport:
                haveTransport = False
            if s.composition is None:
                raise InputParseError('No thermo data found for species: {0!r}'.format(s.label))
            elementsFromSpecies.update(s.composition)
            speciesNameLength = max(speciesNameLength, len(s.label))

        # validate list of elements
        missingElements = elementsFromSpecies - set(self.elements)
        if missingElements:
            raise InputParseError('Undefined elements: ' + str(missingElements))

        speciesNames = ['']
        for i,s in enumerate(self.speciesList):
            if i and not i % 5:
                speciesNames.append(' '*21)
            speciesNames[-1] += '{0:{1}s}'.format(s.label, speciesNameLength+2)

        speciesNames = '\n'.join(speciesNames).strip()

        lines = []
        if header:
            lines.extend(header)

        # Write the gas definition
        lines.append("units(length='cm', time='s', quantity={0!r}, act_energy={1!r})".format(self.quantity_units, self.energy_units))
        lines.append('')
        lines.append('ideal_gas(name={0!r},'.format(name))
        lines.append('          elements="{0}",'.format(' '.join(self.elements)))
        lines.append('          species="""{0}""",'.format(speciesNames))
        if self.reactions:
            lines.append("          reactions='all',")
        if haveTransport:
            lines.append("          transport={0!r},".format(transportModel))
        lines.append('          initial_state=state(temperature=300.0, pressure=OneAtm))')
        lines.append('')

        # Write the individual species data
        lines.append(delimiterLine)
        lines.append('# Species data')
        lines.append(delimiterLine)
        lines.append('')

        for s in self.speciesList:
            lines.append(s.to_cti())

        # Write the reactions
        lines.append(delimiterLine)
        lines.append('# Reaction data')
        lines.append(delimiterLine)

        for i,r in enumerate(self.reactions):
            lines.append('\n# Reaction {0}'.format(i+1))
            lines.append(r.to_cti())

        lines.append('')

        f = open(outName, 'w')
        f.write('\n'.join(lines))

    def showHelp(self):
        print("""
ck2cti.py: Convert Chemkin-format mechanisms to Cantera input files (.cti)

Usage:
    ck2cti --input=<filename>
           [--thermo=<filename>]
           [--transport=<filename>]
           [--id=<phase-id>]
           [--output=<filename>]
           [--permissive]
           [-d | --debug]

Example:
    ck2cti --input=chem.inp --thermo=therm.dat --transport=tran.dat

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.cti'.

The '--permissive' option allows certain recoverable parsing errors (e.g.
duplicate transport data) to be ignored.

""")

    def convertMech(self, inputFile, thermoFile=None,
                    transportFile=None, phaseName='gas',
                    outName=None, quiet=False, permissive=None):
        if quiet:
            logging.basicConfig(level=logging.ERROR)
        else:
            logging.basicConfig(level=logging.INFO)

        if permissive is not None:
            self.warning_as_error = not permissive

        if not os.path.exists(inputFile):
            raise IOError('Missing input file: {0!r}'.format(inputFile))
        try:
            # Read input mechanism files
            self.loadChemkinFile(inputFile)
        except Exception:
            logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n".format(
                            inputFile, self.line_number))
            raise

        if thermoFile:
            if not os.path.exists(thermoFile):
                raise IOError('Missing thermo file: {0!r}'.format(thermoFile))
            try:
                self.loadChemkinFile(thermoFile)
            except Exception:
                logging.warning("\nERROR: Unable to parse '{0}' near line {1}:\n".format(
                                thermoFile, self.line_number))
                raise

        if transportFile:
            if not os.path.exists(transportFile):
                raise IOError('Missing transport file: {0!r}'.format(transportFile))
            lines = open(transportFile, 'rU').readlines()
            self.parseTransportData(lines)

            # Transport validation: make sure all species have transport data
            for s in self.speciesList:
                if s.transport is None:
                    raise InputParseError("No transport data for species '{0}'.".format(s))

        if not outName:
            outName = os.path.splitext(inputFile)[0] + '.cti'

        # Write output file
        self.writeCTI(name=phaseName, outName=outName)
        if not quiet:
            print('Wrote CTI mechanism file to {0!r}.'.format(outName))
            print('Mechanism contains {0} species and {1} reactions.'.format(len(self.speciesList), len(self.reactions)))


def main(argv):
    import getopt
    import sys

    longOptions = ['input=', 'thermo=', 'transport=', 'id=', 'output=',
                   'permissive', 'help', 'debug']

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

    parser = Parser()

    if not options or '-h' in options or '--help' in options:
        parser.showHelp()
        sys.exit(0)

    if '--input' in options:
        inputFile = options['--input']
    else:
        print('Error: no mechanism input file specified')
        sys.exit(1)

    if '--output' in options:
        outName = options['--output']
        if not outName.endswith('.cti'):
            outName += '.cti'
    else:
        outName = None

    permissive = '--permissive' in options
    thermoFile = options.get('--thermo')
    transportFile = options.get('--transport')
    phaseName = options.get('--id', 'gas')

    parser.convertMech(inputFile, thermoFile, transportFile, phaseName,
                       outName, permissive=permissive)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
