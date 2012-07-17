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

import logging
import types

import numpy as np

################################################################################

UNIT_OPTIONS = {'CAL/': 'cal/mol',
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
                'KJOULES/MOLE': 'kJ/mol',
                'MOL': 'mol',
                'MOLE': 'mol',
                'MOLES': 'mol',
                'MOLEC': 'molec',
                'MOLECULES': 'molec'}

ENERGY_UNITS = 'cal/mol'
QUANTITY_UNITS = 'mol'

################################################################################

class InputParseError(Exception):
    """
    An exception class for exceptional behavior involving Chemkin-format
    mechanism files. Pass a string describing the circumstances that caused
    the exceptional behavior.
    """
    pass

################################################################################

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
                         for a in self.composition.iteritems())

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

################################################################################

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

################################################################################

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

    The coefficients are stored internally in the nine-coefficient format, even
    when only seven coefficients are provided.
    """

    def __init__(self, coeffs, Tmin=None, Tmax=None, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        coeffs = coeffs or (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        if len(coeffs) == 7:
            self.cm2 = 0.0; self.cm1 = 0.0
            self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
        elif len(coeffs) == 9:
            self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = coeffs
        else:
            raise InputParseError('Invalid number of NASA polynomial coefficients; '
                                  'should be 7 or 9.')

    def to_cti(self, indent=0):
        prefix = ' '*indent
        vals = self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6
        vals = ['{0: 15.8E}'.format(i) for i in vals]
        lines = ['NASA([{0:.2f}, {1:.2f}],'.format(self.Tmin[0], self.Tmax[0]),
                 prefix+'     [{0}, {1}, {2},'.format(*vals[0:3]),
                 prefix+'      {0}, {1}, {2},'.format(*vals[3:6]),
                 prefix+'      {0}]),'.format(vals[6])]

        return '\n'.join(lines)

################################################################################

class MultiNASA(ThermoModel):
    """
    A set of thermodynamic parameters given by NASA polynomials. This class
    stores a list of :class:`NASA` objects in the `polynomials`
    attribute. When evaluating a thermodynamic quantity, a polynomial that
    contains the desired temperature within its valid range will be used.
    """

    def __init__(self, polynomials=None, Tmin=0.0, Tmax=0.0, comment=''):
        ThermoModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
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

################################################################################

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
                 reversible=True, duplicate=False, fwdOrders=None):
        self.index = index
        self.reactants = reactants
        self.products = products
        self.kinetics = kinetics
        self.reversible = reversible
        self.duplicate = duplicate
        self.fwdOrders = fwdOrders if fwdOrders is not None else {}

    def __str__(self):
        """
        Return a string representation of the reaction, in the form 'A + B <=> C + D'.
        """
        arrow = ' <=> '
        if not self.reversible: arrow = ' -> '
        return arrow.join([' + '.join([str(s) for s in self.reactants]),
                           ' + '.join([str(s) for s in self.products])])

    def to_cti(self, indent=0):
        arrow = ' <=> ' if self.reversible else ' => '
        reactantstr = ' + '.join(str(s) for s in self.reactants)
        productstr= ' + '.join(str(s) for s in self.products)

        kinstr = self.kinetics.to_cti(reactantstr, arrow, productstr, indent)

        k_indent = ' ' * (kinstr.find('(') + 1)

        if self.duplicate:
            kinstr = kinstr[:-1] + ",\n{0}options='duplicate')".format(k_indent)

        if self.fwdOrders:
            order = ' '.join('{0}:{1}'.format(k,v)
                             for (k,v) in self.fwdOrders.items())
            kinstr = kinstr[:-1] + ",\n{0}order='{1}')".format(k_indent, order)

        return kinstr

################################################################################

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

    def __init__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.comment = comment

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
        if hasattr(self, 'efficiencies'):
            return ' '.join('{0}:{1}'.format(mol, eff)
                            for mol,eff in self.efficiencies.iteritems())
        else:
            return ''

################################################################################

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

    def __init__(self, Tdata=None, kdata=None, Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.Tdata = Tdata
        self.kdata = kdata

    def isPressureDependent(self):
        """
        Returns ``False`` since KineticsData kinetics are not
        pressure-dependent.
        """
        return False

################################################################################

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

    def __init__(self, A=0.0, n=0.0, Ea=0.0, T0=1.0, Tmin=None, Tmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
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
        return '[{0.A[0]:e}, {0.n}, {0.Ea[0]}]'.format(self)

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstring = reactantstr + arrow + productstr
        return 'reaction({0!r}, {1})'.format(rxnstring, self.rateStr())

################################################################################

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

    def __init__(self, pressures=None, arrhenius=None, highPlimit=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
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

################################################################################

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

    def __init__(self, coeffs=None, kunits='', Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
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
        rxnstr = reactantstr + ' (+ M)' + arrow + productstr + ' (+ M)'
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

################################################################################

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

    def __init__(self, arrheniusHigh=None, efficiencies=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrheniusHigh = arrheniusHigh
        self.efficiencies = {}
        if efficiencies is not None:
            for mol, eff in efficiencies.iteritems():
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

################################################################################

class Lindemann(ThirdBody):
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
    The former is necessarily one reaction order higher than the latter. For
    the Lindemann model, :math:`F = 1`. A collision efficiency can be used to
    further correct the value of :math:`k(T,P)`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusLow=None, arrheniusHigh=None, efficiencies=None,
                 Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        ThirdBody.__init__(self, arrheniusHigh=arrheniusHigh,
                           efficiencies=efficiencies, Tmin=Tmin, Tmax=Tmax,
                           Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.arrheniusLow = arrheniusLow

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' (+ M)' + arrow + productstr + ' (+ M)'
        prefix = ' '*(indent + 17)
        lines = ['falloff_reaction({0!r},'.format(rxnstr)]
        lines.append(prefix + 'kf={0},'.format(self.arrheniusHigh.rateStr()))
        lines.append(prefix + 'kf0={0},'.format(self.arrheniusLow.rateStr()))
        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))

        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

class Troe(Lindemann):
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
    The former is necessarily one reaction order higher than the latter. A
    collision efficiency can be used to further correct the value of
    :math:`k(T,P)`.

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
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    `alpha`         :class:`Quantity`       The :math:`\\alpha` parameter
    `T1`            :class:`Quantity`       The :math:`T_1` parameter
    `T2`            :class:`Quantity`       The :math:`T_2` parameter
    `T3`            :class:`Quantity`       The :math:`T_3` parameter
    =============== ======================= ====================================

    """

    def __init__(self, arrheniusLow=None, arrheniusHigh=None, efficiencies=None,
                 alpha=0.0, T3=0.0, T1=0.0, T2=None, Tmin=None, Tmax=None,
                 Pmin=None, Pmax=None, comment=''):
        Lindemann.__init__(self, arrheniusLow=arrheniusLow,
                           arrheniusHigh=arrheniusHigh,
                           efficiencies=efficiencies, Tmin=Tmin, Tmax=Tmax,
                           Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.alpha = alpha
        self.T3 = T3
        self.T1 = T1
        self.T2 = T2

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' (+ M)' + arrow + productstr + ' (+ M)'
        prefix = ' '*17
        lines = ['falloff_reaction({0!r},'.format(rxnstr),
                 prefix + 'kf={0},'.format(self.arrheniusHigh.rateStr()),
                 prefix + 'kf0={0},'.format(self.arrheniusLow.rateStr())]

        if self.T2:
            troeArgs = 'A={0.alpha[0]}, T3={0.T3[0]}, T1={0.T1[0]}, T2={0.T2[0]}'.format(self)
        else:
            troeArgs = 'A={0.alpha[0]}, T3={0.T3[0]}, T1={0.T1[0]}'.format(self)
        lines.append(prefix + 'falloff=Troe({0}),'.format(troeArgs))

        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))

        # replace trailing comma
        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

class Sri(Lindemann):
    """
    A kinetic model of a phenomenological rate coefficient k(T, P) using the
    "SRI" formulation of the blending function :math:`F` using either 3 or
    5 parameters. See :ref:`sec-sri-falloff`.

    The attributes are:

    =============== ======================= ====================================
    Attribute       Type                    Description
    =============== ======================= ====================================
    `arrheniusLow`  :class:`Arrhenius`      The Arrhenius kinetics at the low-pressure limit
    `arrheniusHigh` :class:`Arrhenius`      The Arrhenius kinetics at the high-pressure limit
    `efficiencies`  ``dict``                A mapping of species to collider efficiencies
    `A`            ``float``                The :math:`a` parameter
    `B`            ``float``                The :math:`b` parameter
    `C`            ``float``                The :math:`c` parameter
    `D`            ``float``                The :math:`d` parameter
    `E`            ``float``                The :math:`e` parameter
    =============== ======================= ====================================
    """

    def __init__(self, arrheniusLow=None, arrheniusHigh=None, efficiencies=None,
                 A=0.0, B=0.0, C=0.0, D=1.0, E=0.0,
                 Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        Lindemann.__init__(self, arrheniusLow=arrheniusLow,
                           arrheniusHigh=arrheniusHigh,
                           efficiencies=efficiencies, Tmin=Tmin, Tmax=Tmax,
                           Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def to_cti(self, reactantstr, arrow, productstr, indent=0):
        rxnstr = reactantstr + ' (+ M)' + arrow + productstr + ' (+ M)'
        prefix = ' '*17
        lines = ['falloff_reaction({0!r},'.format(rxnstr),
                 prefix + 'kf={0},'.format(self.arrheniusHigh.rateStr()),
                 prefix + 'kf0={0},'.format(self.arrheniusLow.rateStr())]

        if self.D == 1.0 and self.E == 0.0:
            sriArgs = 'A={0.A}, B={0.B}, C={0.C}'.format(self)
        else:
            sriArgs = 'A={0.A}, B={0.B}, C={0.C}, D={0.D}, E={0.E}'.format(self)
        lines.append(prefix + 'falloff=SRI({0}),'.format(sriArgs))

        if self.efficiencies:
            lines.append(prefix + 'efficiencies={0!r},'.format(self.efficiencyString()))

        # replace trailing comma
        lines[-1] = lines[-1][:-1] + ')'
        return '\n'.join(lines)

################################################################################

class TransportData(object):
    geometryFlags = ['atom', 'linear', 'nonlinear']

    def __init__(self, label, geometry, wellDepth, collisionDiameter,
                 dipoleMoment, polarizability, zRot, comment=None):

        assert isinstance(label, types.StringTypes)
        assert int(geometry) in (0,1,2)

        self.label = label
        self.geometry = self.geometryFlags[int(geometry)]
        self.wellDepth = float(wellDepth)
        self.collisionDiameter = float(collisionDiameter)
        self.dipoleMoment = float(dipoleMoment)
        self.polarizability = float(polarizability)
        self.zRot = float(zRot)
        self.comment = comment or '' # @todo: include this in the CTI

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

################################################################################

def fortFloat(s):
    """
    Convert a string representation of a floating point value to a float,
    allowing for some of the peculiarities of allowable Fortran representations.
    """
    s = s.strip()
    s = s.replace('D', 'E').replace('d', 'e')
    s = s.replace('E ', 'E+').replace('e ', 'e+')
    return float(s)

################################################################################

def readThermoEntry(entry, TintDefault):
    """
    Read a thermodynamics `entry` for one species in a Chemkin-format file.
    Returns the label of the species, the thermodynamics model as a
    :class:`MultiNASA` object and the elemental composition of the species.
    """
    lines = entry.splitlines()
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

        a0_high = fortFloat(lines[1][0:15])
        a1_high = fortFloat(lines[1][15:30])
        a2_high = fortFloat(lines[1][30:45])
        a3_high = fortFloat(lines[1][45:60])
        a4_high = fortFloat(lines[1][60:75])

        a5_high = fortFloat(lines[2][0:15])
        a6_high = fortFloat(lines[2][15:30])
        a0_low = fortFloat(lines[2][30:45])
        a1_low = fortFloat(lines[2][45:60])
        a2_low = fortFloat(lines[2][60:75])

        a3_low = fortFloat(lines[3][0:15])
        a4_low = fortFloat(lines[3][15:30])
        a5_low = fortFloat(lines[3][30:45])
        a6_low = fortFloat(lines[3][45:60])
    except (IndexError, ValueError) as err:
        raise InputParseError('Error while reading thermo entry for species {0}'.format(species))

    elements = lines[0][24:44]
    composition = {}
    for i in range(4):
        symbol = elements[5*i:5*i+2].strip()
        count = elements[5*i+2:5*i+5].strip()
        if not symbol:
            continue
        try:
            count = int(float(count))
            if count:
                composition[symbol.capitalize()] = count
        except ValueError:
            pass

    # Construct and return the thermodynamics model
    thermo = MultiNASA(
        polynomials = [
            NASA(Tmin=(Tmin,"K"), Tmax=(Tint,"K"), coeffs=[a0_low, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low]),
            NASA(Tmin=(Tint,"K"), Tmax=(Tmax,"K"), coeffs=[a0_high, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high])
        ],
        Tmin = (Tmin,"K"),
        Tmax = (Tmax,"K"),
    )

    return species, thermo, composition, note

################################################################################

def readKineticsEntry(entry, speciesDict, energyUnits, moleculeUnits):
    """
    Read a kinetics `entry` for a single reaction as loaded from a
    Chemkin-format file. The associated mapping of labels to species
    `speciesDict` should also be provided. Returns a :class:`Reaction` object
    with the reaction and its associated kinetics.
    """

    lines = entry.strip().splitlines()

    # The first line contains the reaction equation and a set of modified Arrhenius parameters
    tokens = lines[0].split()
    A = float(tokens[-3])
    n = float(tokens[-2])
    Ea = float(tokens[-1])
    reaction = ''.join(tokens[:-3])
    revReaction = None
    thirdBody = False

    # Split the reaction equation into reactants and products
    if '<=>' in reaction:
        reversible = True
        reactants, products = reaction.split('<=>')
    elif '=>' in reaction:
        reversible = False
        reactants, products = reaction.split('=>')
    elif '=' in reaction:
        reversible = True
        reactants, products = reaction.split('=')
    else:
        raise InputParseError("Failed to find reactant/product delimiter in reaction string.")

    if '(+M)' in reactants: reactants = reactants.replace('(+M)','')
    if '(+m)' in reactants: reactants = reactants.replace('(+m)','')
    if '(+M)' in products:  products = products.replace('(+M)','')
    if '(+m)' in products:  products = products.replace('(+m)','')

    # Create a new Reaction object for this reaction
    reaction = Reaction(reactants=[], products=[], reversible=reversible)

    # Convert the reactants and products to Species objects using the speciesDict
    for reactant in reactants.split('+'):
        reactant = reactant.strip()
        stoichiometry = 1
        if reactant[0].isdigit():
            # This allows for reactions to be of the form 2A=B+C instead of A+A=B+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(reactant[0])
            reactant = reactant[1:]
        if reactant == 'M' or reactant == 'm':
            thirdBody = True
        elif reactant not in speciesDict:
            raise InputParseError('Unexpected reactant "{0}" in reaction {1}.'.format(reactant, reaction))
        else:
            for _ in range(stoichiometry):
                reaction.reactants.append(speciesDict[reactant])
    for product in products.split('+'):
        product = product.strip()
        stoichiometry = 1
        if product[0].isdigit():
            # This allows for reactions to be of the form A+B=2C instead of A+B=C+C
            # The implementation below assumes an integer between 0 and 9, inclusive
            stoichiometry = int(product[0])
            product = product[1:]
        if product.upper() == 'M' or product == 'm':
            pass
        elif product not in speciesDict:
            raise InputParseError('Unexpected product "{0}" in reaction {1}.'.format(product, reaction))
        else:
            for _ in range(stoichiometry):
                reaction.products.append(speciesDict[product])

    # Determine the appropriate units for k(T) and k(T,P) based on the number of reactants
    # This assumes elementary kinetics for all reactions
    if len(reaction.reactants) + (1 if thirdBody else 0) == 3:
        kunits = "cm^6/(mol^2*s)"
        klow_units = "cm^9/(mol^3*s)"
    elif len(reaction.reactants) + (1 if thirdBody else 0) == 2:
        kunits = "cm^3/(mol*s)"
        klow_units = "cm^6/(mol^2*s)"
    elif len(reaction.reactants) + (1 if thirdBody else 0) == 1:
        kunits = "s^-1"
        klow_units = "cm^3/(mol*s)"
    else:
        raise InputParseError('Invalid number of reactant species for reaction {0}.'.format(reaction))

    # The rest of the first line contains the high-P limit Arrhenius parameters (if available)
    #tokens = lines[0][52:].split()
    tokens = lines[0].split()[1:]
    arrheniusHigh = Arrhenius(
        A = (A,kunits),
        n = n,
        Ea = (Ea, energyUnits),
        T0 = (1,"K"),
    )

    if len(lines) == 1:
        # If there's only one line then we know to use the high-P limit kinetics as-is
        reaction.kinetics = arrheniusHigh
    else:
        # There's more kinetics information to be read
        arrheniusLow = None
        troe = None
        sri = None
        chebyshev = None
        pdepArrhenius = None
        efficiencies = {}
        chebyshevCoeffs = []

        # Note that the subsequent lines could be in any order
        for line in lines[1:]:
            tokens = line.split('/')
            if 'DUP' in line or 'dup' in line:
                # Duplicate reaction
                reaction.duplicate = True

            elif 'LOW' in line or 'low' in line:
                # Low-pressure-limit Arrhenius parameters
                tokens = tokens[1].split()
                arrheniusLow = Arrhenius(
                    A = (float(tokens[0].strip()),klow_units),
                    n = float(tokens[1].strip()),
                    Ea = (float(tokens[2].strip()),"kcal/mol"),
                    T0 = (1,"K"),
                )

            elif 'rev' in line.lower():
                reaction.reversible = False

                # Create a reaction proceeding in the opposite direction
                revReaction = Reaction(reactants=reaction.products,
                                       products=reaction.reactants,
                                       reversible=False)
                tokens = tokens[1].split()
                revReaction.kinetics = Arrhenius(
                    A = (float(tokens[0].strip()),klow_units),
                    n = float(tokens[1].strip()),
                    Ea = (float(tokens[2].strip()),"kcal/mol"),
                    T0 = (1,"K"),
                )

            elif 'ford' in line.lower():
                tokens = tokens[1].split()
                reaction.fwdOrders[tokens[0].strip()] = tokens[1].strip()

            elif 'TROE' in line or 'troe' in line:
                # Troe falloff parameters
                tokens = tokens[1].split()
                alpha = float(tokens[0].strip())
                T3 = float(tokens[1].strip())
                T1 = float(tokens[2].strip())
                try:
                    T2 = float(tokens[3].strip())
                except (IndexError, ValueError):
                    T2 = None

                troe = Troe(
                    alpha = (alpha,''),
                    T3 = (T3,"K"),
                    T1 = (T1,"K"),
                    T2 = (T2,"K") if T2 is not None else None,
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
                    sri = Sri(A=A, B=B, C=C)
                else:
                    sri = Sri(A=A, B=B, C=C, D=D, E=E)

            elif 'CHEB' in line or 'cheb' in line:
                # Chebyshev parameters
                if chebyshev is None:
                    chebyshev = Chebyshev()
                tokens = [t.strip() for t in tokens]
                if 'TCHEB' in line:
                    index = tokens.index('TCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Tmin = float(tokens2[0].strip())
                    chebyshev.Tmax = float(tokens2[1].strip())
                if 'PCHEB' in line:
                    index = tokens.index('PCHEB')
                    tokens2 = tokens[index+1].split()
                    chebyshev.Pmin = (float(tokens2[0].strip()), 'atm')
                    chebyshev.Pmax = (float(tokens2[1].strip()), 'atm')
                if 'TCHEB' in line or 'PCHEB' in line:
                    pass
                elif chebyshev.degreeT == 0 or chebyshev.degreeP == 0:
                    tokens2 = tokens[1].split()
                    chebyshev.degreeT = int(float(tokens2[0].strip()))
                    chebyshev.degreeP = int(float(tokens2[1].strip()))
                    chebyshev.coeffs = np.zeros((chebyshev.degreeT,chebyshev.degreeP), np.float64)
                else:
                    tokens2 = tokens[1].split()
                    chebyshevCoeffs.extend([float(t.strip()) for t in tokens2])

            elif 'PLOG' in line or 'plog' in line:
                # Pressure-dependent Arrhenius parameters
                if pdepArrhenius is None:
                    pdepArrhenius = []
                tokens = tokens[1].split()
                pdepArrhenius.append([float(tokens[0].strip()), Arrhenius(
                    A = (float(tokens[1].strip()),kunits),
                    n = float(tokens[2].strip()),
                    Ea = (float(tokens[3].strip()),"kcal/mol"),
                    T0 = (1,"K"),
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
                pressures = ([P for P, arrh in pdepArrhenius],"atm"),
                arrhenius = [arrh for P, arrh in pdepArrhenius],
            )
        elif troe is not None:
            troe.arrheniusHigh = arrheniusHigh
            troe.arrheniusLow = arrheniusLow
            troe.efficiencies = efficiencies
            reaction.kinetics = troe
        elif sri is not None:
            sri.arrheniusHigh = arrheniusHigh
            sri.arrheniusLow = arrheniusLow
            sri.efficiencies = efficiencies
            reaction.kinetics = sri
        elif arrheniusLow is not None:
            reaction.kinetics = Lindemann(arrheniusHigh=arrheniusHigh, arrheniusLow=arrheniusLow)
            reaction.kinetics.efficiencies = efficiencies
        elif thirdBody:
            reaction.kinetics = ThirdBody(arrheniusHigh=arrheniusHigh)
            reaction.kinetics.efficiencies = efficiencies
        else:
            reaction.kinetics = arrheniusHigh

    return reaction, revReaction

################################################################################

def loadChemkinFile(path, speciesList=None):
    """
    Load a Chemkin-format input file to `path` on disk, returning lists of
    the species and reactions in the Chemkin file.
    """
    speciesDict = {}
    if speciesList is None:
        speciesList = []
    else:
        for species in speciesList:
            speciesDict[species.label] = species

    reactionList = []
    transportLines = []

    def removeCommentFromLine(line):
        if '!' in line:
            index = line.index('!')
            comment = line[index+1:-1]
            line = line[0:index] + '\n'
            return line, comment
        else:
            comment = ''
            return line, comment

    with open(path, 'r') as f:
        line = f.readline()
        while line != '':
            line = removeCommentFromLine(line)[0]
            line = line.strip()
            tokens = line.split()

            if 'SPECIES' in line:
                # List of species identifiers
                index = tokens.index('SPECIES')
                tokens = tokens[index+1:]
                while 'END' not in tokens:
                    line = f.readline()
                    line = removeCommentFromLine(line)[0]
                    line = line.strip()
                    tokens.extend(line.split())

                for token in tokens:
                    if token == 'END':
                        break
                    if token in speciesDict:
                        species = speciesDict[token]
                    else:
                        species = Species(label=token)
                        speciesDict[token] = species
                    speciesList.append(species)

            elif 'THERM' in line:
                # List of thermodynamics (hopefully one per species!)
                line = f.readline()
                TintDefault = float(line.split()[1])
                thermo = ''
                while line != '' and 'END' not in line:
                    line = removeCommentFromLine(line)[0]
                    if len(line) >= 80:
                        if line[79] in ['1', '2', '3', '4']:
                            thermo += line
                            if line[79] == '4':
                                label, thermo, comp, note = readThermoEntry(thermo, TintDefault)
                                try:
                                    speciesDict[label].thermo = thermo
                                    speciesDict[label].composition = comp
                                    speciesDict[label].note = note
                                except KeyError:
                                    logging.warning('Skipping unexpected species "{0}" while reading thermodynamics entry.'.format(label))
                                thermo = ''
                    line = f.readline()

            elif 'REACTIONS' in line:
                # Reactions section
                energyUnits = 'CAL/MOL'
                moleculeUnits = 'MOLES'
                try:
                    energyUnits = tokens[1]
                    moleculeUnits = tokens[2]
                except IndexError:
                    pass

                ENERGY_UNITS = UNIT_OPTIONS[energyUnits]
                QUANTITY_UNITS = UNIT_OPTIONS[moleculeUnits]

                kineticsList = []
                commentsList = []
                kinetics = ''
                comments = ''

                line = f.readline()
                while line != '' and 'END' not in line:

                    lineStartsWithComment = line.startswith('!')
                    line, comment = removeCommentFromLine(line)
                    line = line.strip(); comment = comment.strip()

                    if '=' in line and not lineStartsWithComment:
                        # Finish previous record
                        kineticsList.append(kinetics)
                        commentsList.append(comments)
                        kinetics = ''
                        comments = ''

                    if line: kinetics += line + '\n'
                    if comment: comments += comment + '\n'

                    line = f.readline()

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

                for kinetics, comments in zip(kineticsList, commentsList):
                    reaction,revReaction = readKineticsEntry(kinetics, speciesDict, energyUnits, moleculeUnits)
                    reactionList.append(reaction)
                    if revReaction is not None:
                        reactionList.append(revReaction)

            elif 'TRAN' in line:
                line = f.readline()
                while 'END' not in line:
                    transportLines.append(line)

            line = f.readline()

    # Check for marked (and unmarked!) duplicate reactions
    # Raise exception for unmarked duplicate reactions
    for index1 in range(len(reactionList)):
        reaction1 = reactionList[index1]
        for index2 in range(index1+1, len(reactionList)):
            reaction2 = reactionList[index2]
            if reaction1.reactants == reaction2.reactants and reaction1.products == reaction2.products:
                if reaction1.duplicate and reaction2.duplicate:
                    pass
                elif reaction1.kinetics.isPressureDependent() == reaction2.kinetics.isPressureDependent():
                    # If both reactions are pressure-independent or both are pressure-dependent, then they need duplicate tags
                    # pdep and non-pdep reactions are treated as different, so those are okay
                    raise InputParseError('Encountered unmarked duplicate reaction {0}.'.format(reaction1))

    index = 0
    for reaction in reactionList:
        index += 1
        reaction.index = index

    if transportLines:
        parseTransportData(transportLines, speciesList)

    return speciesList, reactionList

################################################################################

def parseTransportData(lines, speciesList):
    """
    Parse the Chemkin-format transport data in ``lines`` (a list of strings)
    and add that transport data to the species in ``speciesList``.
    """
    speciesDict = dict((species.label, species) for species in speciesList)
    for line in lines:
        line = line.strip()
        if not line or line.startswith('!'):
            continue
        if line.startswith('END'):
            break

        data = line.split()
        if len(data) < 7:
            raise InputParseError('Unable to parse transport data: not enough parameters')
        if len(data) >= 8:
            # comment may contain spaces. Rejoin into a single field.
            comment = ''.join(data[7:]).lstrip('!')
            data = data[:7] + [comment]

        speciesName = data[0]
        if speciesName in speciesDict:
            speciesDict[speciesName].transport = TransportData(*data)

################################################################################

def writeCTI(species,
             reactions=None,
             header=None,
             name='gas',
             transportModel='Mix',
             outName='mech.cti'):

    delimiterLine = '#' + '-'*79
    haveTransport = True
    speciesNameLength = 1
    elements = set()
    for s in species:
        if not s.transport:
            haveTransport = False
        if s.composition is None:
            raise InputParseError('No thermo data found for species: {0!r}'.format(s.label))
        elements.update(s.composition)
        speciesNameLength = max(speciesNameLength, len(s.label))

    speciesNames = ['']
    for i,s in enumerate(species):
        if i and not i % 5:
            speciesNames.append(' '*21)
        speciesNames[-1] += '{0:{1}s}'.format(s.label, speciesNameLength+2)

    speciesNames = '\n'.join(speciesNames).strip()

    lines = []
    if header:
        lines.extend(header)

    # Write the gas definition
    lines.append("units(length='cm', time='s', quantity={0!r}, act_energy={1!r})".format(QUANTITY_UNITS, ENERGY_UNITS))
    lines.append('')
    lines.append('ideal_gas(name={0!r},'.format(name))
    lines.append('          elements="{0}",'.format(' '.join(elements)))
    lines.append('          species="""{0}""",'.format(speciesNames))
    if reactions:
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

    for s in species:
        lines.append(s.to_cti())

    # Write the reactions
    lines.append(delimiterLine)
    lines.append('# Reaction data')
    lines.append(delimiterLine)

    for i,r in enumerate(reactions):
        lines.append('\n# Reaction {0}'.format(i+1))
        lines.append(r.to_cti())

    lines.append('')

    f = open(outName, 'w')
    f.write('\n'.join(lines))

################################################################################

def showHelp():
    print """
ck2cti.py: Convert Chemkin-format mechanisms to Cantera input files (.cti)

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.cti'.

Usage:
    ck2cti --input=<filename>
           [--thermo=<filename>]
           [--transport=<filename>]
           [--id=<phase-id>]
           [--output=<filename>]
           [-d | --debug]

Example:
    ck2cti --input=chem.inp --thermo=therm.dat --transport=tran.dat

"""

################################################################################

def convertMech(inputFile, thermoFile=None,
                transportFile=None, phaseName='gas',
                outName=None):
    # Read input mechanism files
    species, reactions = loadChemkinFile(inputFile)

    if thermoFile:
        species, _ = loadChemkinFile(thermoFile, species)

    if transportFile:
        lines = open(transportFile).readlines()
        parseTransportData(lines, species)

    if not outName:
        outName = os.path.splitext(inputFile)[0] + '.cti'

    # Write output file
    writeCTI(species, reactions, name=phaseName, outName=outName)
    print 'Wrote CTI mechanism file to {0!r}.'.format(outName)
    print 'Mechanism contains {0} species and {1} reactions.'.format(len(species), len(reactions))

################################################################################

if __name__ == '__main__':
    import getopt
    import sys
    import os.path

    longOptions = ['input=', 'thermo=', 'transport=', 'id=', 'output=',
                   'help', 'debug']

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'dh', longOptions)
        options = dict()
        for o,a in optlist:
            options[o] = a

        if args:
            raise getopt.GetoptError('Unexpected command line option: ' +
                                     repr(' '.join(args)))

    except getopt.GetoptError as e:
        print 'ck2cti.py: Error parsing arguments:'
        print e
        print 'Run "ck2cti.py --help" to see usage help.'
        sys.exit(1)

    if not options or '-h' in options or '--help' in options:
        showHelp()
        sys.exit(0)

    if '--input' in options:
        inputFile = options['--input']
    else:
        print 'Error: no mechanism input file specified'
        sys.exit(1)

    if '--output' in options:
        outName = options['--output']
        if not outName.endswith('.cti'):
            outName += '.cti'
    else:
        outName = None

    thermoFile = options.get('--thermo')
    transportFile = options.get('--transport')

    phaseName = options.get('--id', 'gas')

    convertMech(inputFile, thermoFile, transportFile, phaseName, outName)
